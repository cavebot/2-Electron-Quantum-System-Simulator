C#
C# 
C#             CPC/ SUBMITTED VERSION  
C#             BSPCI2E PACKAGE : 1.0    / F77
C#
C*             PROGRAM NUMBER 8          CS1PH
C#
C#LAAN2002
C#######################################################################
      PROGRAM CS1PH
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      INCLUDE "parameter.2e.inc"

      PARAMETER ( DSQ2     = 1.414 213 562 4 D+00 )
      PARAMETER ( PI       = 3.141 592 653 589 793 D+00 )
      PARAMETER ( ENAU     = 27.211 396 181  D+00 )
      PARAMETER ( T0       = 2.418 884 326 555 53 10D-017)
      PARAMETER ( A0       = 0.529 177 249 24D-08)
      PARAMETER ( ALPHA    = 1.0D+00/137.035 989 561 D+00 )
      PARAMETER ( CS1_MB   = 28.00286 D+00 )
      PARAMETER ( NPHOTONS = 1 )

      INTEGER LINITIAL, LFINAL, LOI, LOF

      DIMENSION ENI(NWF2E),ENF(NWF2E)
      DIMENSION EHF(NS-2), DMX(NWF2E), DDMX(NWF2E)
      INTEGER NINP, NOUT, ND2E
      integer   n_state_initial, n_state_final
      CHARACTER*6 GAUGE, dip_type
      CHARACTER*100 argv

      
C..............................................................


      NINP = 9
      NOUT = 16
      ND2E = 17

C..............................................................

      CALL GETARG(1, ARGV)
      READ(ARGV,*) linitial
      CALL GETARG(2, ARGV)
      READ(ARGV,*) lfinal
      CALL GETARG(3, ARGV)
      READ(ARGV,*) gauge
      CALL GETARG(4, ARGV)
      READ(ARGV,*) dip_type     !'cs', 'os'

C
      n_state_initial = 1        ! from first state of the initial symmetry


      INPMODE = 1
      IF(gauge.EQ.'v')          INPMODE = 0     
      
C.............................................................

      OPEN(NINP, FILE='dat/en1e.dat')
      READ(NINP,*) LMIN, LMAX

      IF((LMIN-1).NE.0) THEN
         
         WRITE(*,*) ' ION ENERGIES NOT FOR L = 0. EXITING..'
         STOP
      ENDIF

         DO  IE = 1, NS-2

            EHF(IE) = 0.0D+00
         ENDDO

         READ(NINP,*) NMX, NVALENCE

         DO  IE = 1, NMX

            READ(NINP,*)  EHF(IE)
         ENDDO

      CLOSE(NINP)


C======================================================

C#    1ST ( N = 1 ) IONIZATION THRESHOLD

      ETHR1 = EHF(NVALENCE) 

      WRITE(*,*) 'FIRST IONIZATION THRESHOLD : ETHR1  = ', ETHR1


C#
C#     READ  INFORMATION FROM DME FILE
C#


      IF(linitial.LT.lfinal ) THEN

         CALL dmxfile(nd2e, "dat", "he","bin", gauge, linitial, lfinal)
 
C         CALL DMX2EFILE(ND2E, LIN, INPMODE)

         READ(ND2E)      MODE
         READ(ND2E)      LOI, LOF, NSI, NSF

         write(*,*) nsi, nsf, nwf2e

         READ(ND2E)    ( ENI(K), K = 1, NSI)
         READ(ND2E)    ( ENF(K), K = 1, NSF)

         n_state_final   = nsf

         DO I = 1, n_state_initial
           READ(ND2E) ( DMX(NSR), NSR = 1, NSF)
         END DO

      ELSE

!     dmx file!
         CALL dmxfile(nd2e, "dat", "he","bin", gauge, lfinal, linitial)
!         CALL DMX2EFILE(ND2E, LFIN, INPMODE)

         READ(ND2E)      MODE
         READ(ND2E)      LOF, LOI, NSF, NSI
         READ(ND2E)    ( ENF(K), K = 1, NSF)
         READ(ND2E)    ( ENI(K), K = 1, NSI)

         n_state_final   = nsf
         DO I = 1, n_state_final
            READ(ND2E) ( DDMX(NSR), NSR = 1, NSI)
            DMX(I) = DDMX(n_state_initial)
         END DO



      ENDIF


      
      IF(INPMODE.NE.MODE) THEN 
         WRITE(*,*)'# cs2e1ph_b:: incosistent gauge'
         WRITE(*,*)'# cs2e1ph_b::          mode = ', mode
         WRITE(*,*)'# cs2e1ph_b::      inp_mode = ', inpmode
         WRITE(*,*)'# cs2e1ph_b aborts            '
         STOP
      ENDIF




C........................................................

     

C#        RENORMALIZATION OF THE INITIAL STATE ( RHO = SQRT( 2 /DE(I) )

         RHO_I = 1.0D+00
         IF( ENI(n_state_initial).GT.ETHR1) THEN
            RHO_I = 1.0D+00/SQRT( ABS( ENI(n_state_initial + 1) 
     &                               - ENI(n_state_final   - 1))/4.0)
         ENDIF


         CALL CS1PHFILE(16, L, MODE) 

         sf_n2 = 0.0D+00
         sf_n1 = 0.0D+00
         sf_0 = 0.0D+00
         sf_1 = 0.0D+00
         sf_2 = 0.0D+00
      
         
         DO J = 1,  n_state_final


C#   NORMALIZE CONTINUUM FINAL STATES ( RHO = SQRT( 2 /DE(J) )

         DE_AU  = ABS( ENF(J) - ENI(n_state_initial) )/2.0

         RHO_J = 1.0D+00 / SQRT( ABS( ENF(J + 1) - ENF(J - 1))/4 )


         WPH = DE_AU 

         IF(MODE.EQ.0) THEN

            WPH = 1.0D+00 / DE_AU 

         ENDIF

C#
C#        OS IN A.U. F_AB
C#

         F_BB = ( 2.0D0/3.0D0) * ( 2*LOF + 1) * WPH * (DMX(J) * DMX(J))  


         sf_n2 = sf_n2 + f_bb / ( 2.0D+00 * de_au )**2
         sf_n1 = sf_n1 + f_bb / ( 2.0D+00 * de_au )
         sf_0 =  sf_0  + f_bb
         sf_1 =  sf_1  + f_bb * ( 2.0D+00 * de_au )
         sf_2 =  sf_2  + f_bb * ( 2.0D+00 * de_au )**2

C#  CS ---> CROSS SECTIONS

         IF(dip_type.EQ.'cs') THEN


            IF(ENF(J).GT.ETHR1) THEN

C....................................................................

C#
C#  CS_N ( IN A.U.) =  2 * PI * ( 2 * PI * ALPHA )^N
C#                       * WPH^N * | D_L^(N) ( IN A.U.) |^2
C#
C#  CS_N ( SI )     =  2 * PI * ( 2 * PI * ALPHA )^N
C#                       * TO**(N-1) * A0**(2N) 
C#                       * WPH^N * | D_L^(N) ( IN A.U.) |^2



               CS_N_AU  =  2 * PI * ( 2 * PI * ALPHA )**NPHOTONS
               CS_N_SI  =  CS_N_AU * T0**(NPHOTONS-1) * A0**(2*NPHOTONS) 


               CS_1_AU = CS_N_AU * WPH**NPHOTONS * ( DMX(J) * RHO_J )**2
               CS_1_SI = CS_N_SI * WPH**NPHOTONS * ( DMX(J) * RHO_J )**2


C               WRITE(*,*) ' NPHOTONS = ', NPHOTONS
C               WRITE(*,*) ' CS_N_AU  = ', CS_N_AU
Cq               WRITE(*,*) ' CS_N_SI  = ', CS_N_SI


C#....................................

               CS_1 =  PI * ( 2 * PI * ALPHA ) * F_BB * RHO_J * RHO_J
C#
 
               WRITE(16,5) DE_AU * ENAU, CS_1 * CS1_MB, CS_1,
     &                     CS_1 * CS1_MB * 1.0D-18, DMX(J) * RHO_J


               WRITE(*,*) DE_AU * ENAU, CS_1_SI *1.0D+18, CS_1_AU



C               WRITE(16,5) DE_AU * ENAU, DMX(J)**2



            ENDIF


         ELSE IF(dip_type.eq.'os') THEN


C#        OS 

            IF( ENF(J).LT.ETHR1 ) THEN

               WRITE(*,'(I3,1x,5(G12.5))') 
     &          j, sf_n2, sf_n1, sf_0, sf_1, sf_2

               WRITE(16,5)  ENF(J) * ENAU * 0.5D0, F_BB * CS1_MB, F_BB,
     &                      F_BB * CS1_MB * 1.0D-18 
 
            ENDIF
            
         ELSE
            WRITE(*,*) ' cs2e1ph_b: no calculation has selected '
            WRITE(*,*) '     choice = cs (cross sections)       '
            WRITE(*,*) '            = os (oscillator strengths) '
            
         ENDIF


      END DO

C         WRITE(*,*)"# cs1ph::   sum of oscillators  S_f = ", sum_f

C      WRITE(*,*)"   S_2         S_1         S_0         S_1         S_2"
C      WRITE(*,'(5(G12.5))') sf_n2, sf_n1, sf_0, sf_1, sf_2

      CLOSE(16) 

 5    FORMAT(2X,1P8E15.7,2X,1P8E15.7,2X,1P8E15.7,2X)
      END
C#######################################################################
C#EOF
