C#
C# 
C#             CPC/ SUBMITTED VERSION  
C#             BSPCI2E PACKAGE : 1.0
C#
C#             PROGRAM NUMBER 9          CS2PH
C#
C#
C#######################################################################
      PROGRAM CS2PH

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

      INCLUDE "parameter.2e.inc"

      PARAMETER ( DSQ2     = 1.414 213 562 4 D+00 )
      PARAMETER ( PI       = 3.141 592 653 589 793 D+00 )
      PARAMETER ( ENAU     = 27.211 396 181  D+00 )
      PARAMETER ( T0       = 2.418 884 326 555 53 10D-017)
      PARAMETER ( A0       = 0.529 177 249 24D-08)
      PARAMETER ( ALPHA    = 1.0D+00/137.035 989 561 D+00 )
      PARAMETER ( CS2_SI   = 2.50547D-052 )
      PARAMETER ( NPHOTONS = 2 )


      DIMENSION DMX(NWF2E,NWF2E,2),ENI(NWF2E)
      DIMENSION EHF(NS-2),ENF(NWF2E),ENM(NWF2E)
      INTEGER   NINP, NOUT, ND2E
      INTEGER   linitial , loi, n_state_initial
      INTEGER   lfinal   , lof, n_state_final
      INTEGER   lintermed, lom, n_state_intermed
C
      CHARACTER*6  gauge
C      CHARACTER*6  dip_type
      CHARACTER*100 argv

C..............................................................


      NINP = 9
      NOUT = 16
      ND2E = 17


C..............................................................

C      OPEN(NINP,FILE='inp/cs2ph.inp',STATUS='OLD')
C      READ(NINP, *) LIN, LFIN, GAUGE
C      READ(NINP, *) NIN, NFIN
C      READ(NINP, '(A14)') CHOICE

      CLOSE(NINP) 

      CALL GETARG(1, ARGV)
      READ(ARGV,*) linitial
      CALL GETARG(2, ARGV)
      READ(ARGV,*) lfinal
      CALL GETARG(3, ARGV)
      READ(ARGV,*) gauge
C      CALL GETARG(4, ARGV)
C      READ(ARGV,*) dip_type


      n_state_initial = 1 

C..............................................................
C#
C#                READ 1-E ENERGY DATA  AND GET E+
C#

      OPEN(NINP,FILE='dat/en1e.dat')
      
      READ(NINP, *) LMIN, LMAX

      IF((LMIN-1).NE.0) THEN
         
         WRITE(*, *) ' ION ENERGIES NOT FOR L = 0. EXITING..'
         STOP
      ENDIF

         DO  IE = 1, NS - 2

            EHF(IE) = 0.0D+00
         ENDDO

         READ(NINP,*) NMX, NVALENCE

         DO  IE = 1, NMX

            READ(NINP, *)  EHF(IE)
         ENDDO

      CLOSE(NINP)


C#  FIRST IONIZATION THRESHOLD   E+


      ETHR1 = EHF(NVALENCE) 


C#
C#  I.E. FOR HELIUM  SHOULD BE     E+(HE) = - 4.0 RYDBERG
C.........................................................


C#
C#         GET TRANSITION DATA FROM 2-E DMX
C#

      INPMODE = 1
      IF(GAUGE.EQ.'v')          INPMODE = 0 


C#                FIRST TRANSITION
C#
C#                       HW  
C#   |INITIAL STATES > -----> |INTERMEDIATE STATES >
C#


      lintermed = linitial + 1


C#           DMX2E_1( N L , M L + 1) 

C      CALL DMX2EFILE(ND2E, linitial, INPMODE)

      CALL dmxfile(nd2e, "dat", "he","bin", gauge, linitial, lintermed)

      READ(ND2E)  MODE  
      READ(ND2E)  LOI, LOM, NSI, NSM
      READ(ND2E) (ENI(NSC), NSC = 1, NSI)
      READ(ND2E) (ENM(NSR), NSR = 1, NSM)


      WRITE(*,*) ' READ ENERGIES'
      WRITE(*,*) ' NSI = ', NSI
      WRITE(*,*) ' NSM = ', NSM

      DO NSR = 1, NSI
         READ(ND2E) ( DMX(NSC, NSR, 1), NSC = 1, NSM)
      ENDDO


      WRITE(*,*)'# cs2e2ph_b:: first transition : ', loi, lom

      CLOSE(ND2E)




      IF(INPMODE.NE.MODE) THEN 
         WRITE(*,*)'# cs2e2ph_b:: incosistent gauge'
         WRITE(*,*)'# cs2e2ph_b::          mode = ', mode
         WRITE(*,*)'# cs2e2ph_b::      inp_mode = ', inpmode
         WRITE(*,*)'# cs2e2ph_b aborts            '
         STOP
      ENDIF




C#                        SECOND NUMBER 
C#                                HW
C#   | INTERMEDIATE STATES > -----------> | FINAL STATES >




      IF( linitial.NE.lfinal ) THEN


C#       DMX_2( M L + 1, K L + 2 ) 

C#         CALL DMX2EFILE(ND2E, linitial + 1, INPMODE)

         CALL dmxfile(nd2e, "dat", "he","bin", gauge, lintermed, lfinal)


         READ(ND2E)  MODE
         READ(ND2E)  LOM, LOF, NSM, NSF
         READ(ND2E) (ENM(NSC), NSC = 1, NSM)
         READ(ND2E) (ENF(NSR), NSR = 1, NSF)


         DO NSC = 1, NSM
            READ(ND2E) ( DMX(NSC, NSR, 2), NSR = 1, NSF)
         ENDDO

         WRITE(*,*)'# cs2e2ph_b:: second transition : ', lom, lof

      ELSE


C#       DMX_2( K L + 2, M L + 1) 


         CALL dmxfile(nd2e, "dat", "he","bin", gauge,linitial,lintermed)


         READ(ND2E)  MODE
         READ(ND2E)  LOF, LOM, NSF, NSM
         READ(ND2E) (ENF(NSC), NSC = 1, NSF)
         READ(ND2E) (ENM(NSR), NSR = 1, NSM)

         
         DO NSR = 1, NSF
            READ(ND2E) ( DMX( NSC, NSR, 2), NSC = 1, NSM)            
         ENDDO

         WRITE(*,*)'# cs2e2ph_b:: second transition : ', lom, lof
      ENDIF
      
      CLOSE(ND2E)


      IF(INPMODE.NE.MODE) THEN 
         WRITE(*,*)'# cs2e2ph_b:: incosistent gauge'
         WRITE(*,*)'# cs2e2ph_b::          mode = ', mode
         WRITE(*,*)'# cs2e2ph_b::      inp_mode = ', inpmode
         WRITE(*,*)'# cs2e2ph_b aborts            '
         STOP
      ENDIF


C................................................................
C#
C#            TWO-PHOTON GENERALIZED CROSS SECTIONS 
C#

C      n_state_final = nsf

      n_state_final  =  300

C      ENI(N_STATE_INITIAL) = -5.806D+00           

      WRITE(*,*) '        E+     = ',     ETHR1 * 0.5D+00 * ENAU
      WRITE(*,*) ' E_INITIAL(IN) = ', ENI(N_STATE_INITIAL)*0.5D+00*ENAU
      WRITE(*,*) ' E_FINAL(IN)   = ',    ENF(1) * 0.5D+00 * ENAU
      WRITE(*,*) ' E_FINAL(FIN)  = ', ENF(N_STATE_FINAL)*0.5D+00* ENAU



C#  OUTPOUT FILE FOR THE 2-PHOTON CROSS SECTIONS DATA


      CALL CS2PHFILE(16, linitial, lfinal, mode) 


C.................................

        CS_N_AU  = 2 * PI * ( 2 * PI * ALPHA )**NPHOTONS
        CS_N_SI  = CS_N_AU * T0**(NPHOTONS-1) * A0**(2*NPHOTONS) 


        WRITE(*,*) ' NPHOTONS = ', NPHOTONS
        WRITE(*,*) ' CS_N_AU  = ', CS_N_AU
        WRITE(*,*) ' CS_N_SI  = ', CS_N_SI

C................................

C#  LOOP OVER THE FINAL STATES



      DO J = 1, N_STATE_FINAL

C#  PHOTON ENERGY IN A.U.
 
        WPH  =  0.25D+00 * ABS( ENF(J) -  ENI(N_STATE_INITIAL) )

C#  AVOID THE POLE INTO THE CONTINUUM (ATI CASE)


        IF(WPH.GT.ABS(0.5D+00 *(ETHR1-ENI(N_STATE_INITIAL)))) GOTO 200
  
C#...........................................................
C#
C#           INTERMEDIATE STATES SUMMATION
C#

        DMX_2     = 0.0D+00
        DMX_2_POS = 0.0D+00
        DMX_2_NEG = 0.0D+00
        
        DO IM = 1, NSM

C#      DETUNING IN A.U.

           DE = ENI(N_STATE_INITIAL) * 0.5D+00 + WPH - ENM(IM)* 0.5D+00
           DE_IM = 1.0D+00
           DE_MF = 1.0D+00
           
           IF(MODE.EQ.0) THEN
              DE_IM = ( ENI(N_STATE_INITIAL) - ENM(IM) ) * 0.5D+00
              DE_MF = ( ENM(IM)   -  ENF(J) ) * 0.5D+00
           ENDIF
           
           DE_IF = DE 

C* DE_IM * DE_MF
C#
           DMX_2  = DMX(IM, N_STATE_INITIAL,1) * DMX(IM, J, 2) / DE_IF

C#      FOR NUMERICAL REASONS NEXT SEGMENT
 

          
            IF(DMX_2.GT.0D+00) THEN
               DMX_2_POS = DMX_2_POS + DMX_2
            ELSE            
               DMX_2_NEG = DMX_2_NEG + DMX_2
            ENDIF

        ENDDO

C............................................................
C#
C#             NORMALIZATION FACTOR  IN A.U.


C# BOUND STATES

        DENS_J = 1.0D+00

C# CONTINUUM STATES

        IF(ENF(J).GT.ETHR1) THEN

           DENS_J = 1.0D+00/SQRT( ABS( ENF(J+1) - ENF(J-1) ) * 0.25D+00)

        ENDIF

C        DMX_2 = DENS_J * DENS_J * ( DMX_2_POS + DMX_2_NEG )

          DMX_2 = DENS_J * ( DMX_2_POS + DMX_2_NEG )



C....................................................................
C#
C#  CS_N ( IN A.U.) =  2 * PI * ( 2 * PI * ALPHA )^N
C#                       * WPH^N * | D_L^(N) ( IN A.U.) |^2
C#
C#  CS_N ( SI )     =  2 * PI * ( 2 * PI * ALPHA )^N
C#                       * TO**(N-1) * A0**(2N) 
C#                       * WPH^N * | D_L^(N) ( IN A.U.) |^2
C#
C#
C#  CS_2( CM^4)      = 2.50547 x 10^{-52} 
C#                   * WPH(AU)**2  * | D^(2)(AU) |^2
C#
C#


C#    LENGTH GAUGE
          DWPH = WPH**(2*mode-1)
C#  VELOCITY GAUGE
!          IF(MODE.EQ.0) THEN
!             DWPH = 1.0D+00 / WPH
!          ENDIF
          

        CS_2_AU = CS_N_AU * DWPH**2  * DMX_2**2
        CS_2_SI = CS_N_SI * DWPH**2  * DMX_2**2


        IF(ENF(J).GT.ETHR1) THEN

C           WRITE(16,5)  WPH * ENAU, DMX_2**2, CS_2_SI, CS_2_AU
           WRITE(16,5)  WPH * ENAU, CS_2_SI, CS_2_AU
        ENDIF

      ENDDO

 200  CONTINUE

      CLOSE(16) 

C#                    FORMAT STATEMENTS
C........................................................

    5 FORMAT(2X,1P8E15.7,2X,1P8E15.7,2X,1P8E15.7)

      END
C#######################################################################
C#EOF
