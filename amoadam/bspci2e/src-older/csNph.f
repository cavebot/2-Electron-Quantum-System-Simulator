C#                                                                       
C#            GENERALISED N-PHOTON CROSS SECTIONS
C#
C#25082004LAANKYOTO::  changes to the output have been made 
C#                     units include has been added
C#
C#
C#25082004LAANKYOTO: 
C#                  (1) input has been changed such that 
C#                      icsNph.f and csNph.f to use the
C#                      same input file inp/csNph.inp
C#
C#
C#29082004LAANKYOTO: 
C#                  (1) in output added the real and the im part
C#                      of the multiphoton dmx
C#                  (2) a subroutine for i/o added
C#                  (3)                     
C#
C#
C#                                                                       
      PROGRAM GCSNPHOTON

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               

      PARAMETER ( MAXPHOT = 4, MAXPATH =10)
      PARAMETER ( MAXDIM = 2000, MAXL = 4, MAXNBPTS = 300 )

      INCLUDE "units.inc"

      COMPLEX*16 DMX, D_N, SIGMA
      REAL*8 S_N_SI,S_N_AU
 
      LOGICAL ANGULPARTINCLUDED 
      INTEGER CHEMIN,ARBRELM,POLAR,FINALSTATE
      INTEGER ATLEVEL
      INTEGER MODE, INP_MODE
      INTEGER NINP, NCSNPH
      CHARACTER*14 MONITORFILE,DIPFILE,ENERFILE,OUTPUTFILE
      CHARACTER*14 ATI
      DIMENSION FINALSTATE(MAXPHOT,2)
      DIMENSION DIPR(MAXDIM,MAXDIM,MAXL)
      DIMENSION DMX( MAXNBPTS )
      DIMENSION EPSILON( MAXL ), EPS( MAXNBPTS ) 
      DIMENSION DELTA( MAXL ), EPSILONMIN( MAXL ), EPSILONMAX( MAXL )
      DIMENSION ENERGYSPACING( MAXL )
      DIMENSION DMX_R( MAXNBPTS )
      DIMENSION DMX_I( MAXNBPTS )
      DIMENSION ENERG(MAXDIM,0:MAXL),RNORME(MAXDIM,0:MAXL) 
      DIMENSION PHASE(MAXDIM,0:MAXL),AN(MAXDIM,0:MAXL)
      INTEGER NDIMENSION
C...........

      CHARACTER*100 ARGV
      CHARACTER*1   SMODE

C..........

      COMMON /PATHS/ ARBRELM(MAXPATH,MAXPHOT,2),ANGULAR(MAXPATH)

C...................  input data section.
c


c............ start here

      WRITE(*,*) '##################################################'
      WRITE(*,*) '#'
      WRITE(*,*) '#'
      WRITE(*,*) '#   CSNPH PACKAGE V.0: csNph_l:   '
      WRITE(*,*) '#'

      NINP = 9
      NCSNPH = 58


      OUTPUTFILE = "dat/csNph.dat"


C......................................................

C GET THE PARAMETERS FROM THE COMMON LINE ARGUMENT

      CALL GETARG(1, ARGV)
      READ(ARGV,*) NBPHOT
      CALL GETARG(2, ARGV)
      READ(ARGV,*) SMODE


      IF(SMODE.EQ."v") THEN
         INP_MODE = 0
      ELSE IF(SMODE.EQ."l") THEN
      INP_MODE = 1           
      ELSE
         INP_MODE = 1

         WRITE(*,*) '# csNph_l::   only the inputs l or v are available'
         WRITE(*,*) '# csNph_l::        set default mode, length gauge '

      ENDIF
      

      WRITE(*,*) '# csNph_l::        inp_mode =  ', inp_mode



      OPEN(NINP,FILE='inp/csNph.inp') 
      READ(NINP,*) ATI
      READ(NINP,*) NINIT
      READ(NINP,*) LINIT
      READ(NINP,*) MINIT
      READ(NINP,*) POLAR
      READ(NINP,*) WMIN
      READ(NINP,*) WMAX
      READ(NINP,*) EPSFACTMIN
      READ(NINP,*) EPSFACTMAX
      READ(NINP,*) NBPTS
      CLOSE(NINP)

C..... get the ion's first ionization  threshold



      CALL READ_ION_THRESHOLD(E_IONIC)


C......                N-photon cross section formulas
C#
C#
C#  CS_N ( IN A.U.) =  2 * PI * ( 2 * PI * ALPHA )^N
C#                       * WPH^N * | D_L^(N) ( IN A.U.) |^2
C#
C#  CS_N ( SI )     =  2 * PI * ( 2 * PI * ALPHA )^N
C#                       * TO**(N-1) * A0**(2N) 
C#                       * WPH^N * | D_L^(N) ( IN A.U.) |^2


        UNIT_N = 1.0D+00         
      IF(NBPHOT.EQ.1) THEN         
         UNIT_N = 1.0D+18
      ENDIF

      CS_N_AU  =  2 * PI * ( 2 * PI * ALPHA )**NBPHOT
      CS_N_SI  =  CS_N_AU * T0**(NBPHOT-1) * A0**(2*NBPHOT) 

c...............


      WRITE(*,*) "#CSNPH_L::   read data files: en.dat dmx.dat"


      CALL READDIPO(DIPFILE,ENERFILE,DIPR,ENERG,RNORME,PHASE,AN,
     &   NDIMENSION,LMIN,LMAX,MAXDIM,MAXL,ANGULPARTINCLUDED,MODE)

C#      ENERG( NINIT,   LINIT ) = -0.832235d+00
C#      ENERG( NINIT, LINIT+1 ) = -0.671735

      IF(INP_MODE.NE.MODE) THEN

         WRITE(*,*) "#CSNPH_P::"
         WRITE(*,*) "# GAUGE MODES DIFFER IN csNph.inp and indmxNph.dat"
         WRITE(*,*) "#CSNPH_P::   dmxNph.dat:       MODE = ", MODE
         WRITE(*,*) "#CSNPH_P::    csNph.inp:   INP_MODE = ", INP_MODE
         STOP
      ENDIF

      CALL CHECK(LMIN,LMAX,NDIMENSION,NINIT,NBPHOT,LINIT,MINIT)
      CALL PATHLM(NBPHOT,LINIT,MINIT,POLAR,NBPATH,ARBRELM)
      CALL COUNTFINALSTATE(NBPHOT,NBPATH,ARBRELM,FINALSTATE,NUMBER,
     S                     MAXPATH,MAXPHOT)


      MONITORFILE = 'out/csNph.out'

      OPEN(14,FILE=MONITORFILE,STATUS='UNKNOWN')
      OPEN(44,FILE='dat/csNph.1.dat',STATUS='UNKNOWN')
      OPEN(54,FILE='dat/csNph.2.dat',STATUS='UNKNOWN')


      WRITE(*,*) '#'
      WRITE(*,*) '#               NBPHOT  = ', NBPHOT
      WRITE(*,*) '#       GAUGE      MODE = ', MODE 
      WRITE(*,*) '#'

      IF(MODE.EQ.0)          WRITE(*,*) '#              VELOCITY'
      IF(MODE.EQ.1)          WRITE(*,*) '#                LENGTH'

      WRITE(*,*) '#'
      WRITE(*,*) '#    CONVERSION FACTORS AU   TO  SI:'
      WRITE(*,*) '#'
      WRITE(*,*) '#              CS_N_AU  = ', CS_N_AU
      WRITE(*,*) '#              CS_N_SI  = ', CS_N_SI
      WRITE(*,*) '#'


CC LES CALCULS SERONT FAITS POUR DES PHOTONS D'ENERGIE ALLANT DE 
c WMIN AC WMAX.


C#    start calculation of cs for the partial wave (lend)

      DO NBFINAL = 1, NUMBER

	 LEND = FINALSTATE(NBFINAL,1)
	 MEND = FINALSTATE(NBFINAL,2)

	 CALL PATHLM(NBPHOT,LINIT,MINIT,POLAR,NBPATH,ARBRELM)
	 CALL SELECTPATH(NBPHOT,LEND,MEND,NBPATH,ARBRELM)
         CALL ANGULCOEFF(ANGULAR,POLAR,NBPHOT,NBPATH,ARBRELM,
     S                   MAXPATH,MAXPHOT,ANGULPARTINCLUDED)

	 WRITE(14,*) '#  csNph::  FINAL SYMMETRY  : ', LEND, MEND
	 WRITE(14,*) '#  csNph::  NUMBER OF PATHS :',  NBPATH

	 DO N = 1, NBPATH
	    WRITE(14,FMT='(8(1X,I3))')(ARBRELM(N,J,1),J=1,NBPHOT+1)
	    WRITE(14,FMT='(8(1X,I3))')(ARBRELM(N,J,2),J=1,NBPHOT+1)
	 ENDDO
         
C....................

C#  NO - ATI

      WRITE(*,*) '#      n_init = ', NINIT
      WRITE(*,*) '#      l_init = ', LINIT
      WRITE(*,*) '#         e_0 = ', ENERG( NINIT,  LINIT ) * enau
      WRITE(*,*) '#         e_1 = ', ENERG( NINIT, LINIT+1) * enau
      WRITE(*,*) '#          E+ = ',  E_IONIC * enau
      WRITE(*,*) '#        Wmin = ',  WMIN * enau
      WRITE(*,*) '#        Wmax = ',  WMAX * enau 
      WRITE(*,*) '#         Nph = ',  nbphot

C      de_3s3p = abs(ENERG( NINIT, LINIT ) - ENERG( NINIT, LINIT+1))
C      WRITE(*,*) '#                   de_3s3p = ',  de_3s3p * enau


      IF(ATI.EQ.'NO'.or.ati.eq.'no') THEN

         WRITE(*,*) '#                 NO ATI CROSS SECTIONS  ' 
         WRITE(*,*) '#                WMIN, WMAX RECALCULATED '
         
         WMIN = ABS( ENERG( NINIT, LINIT ) - E_IONIC )/DBLE(NBPHOT) 

         IF(NBPHOT.NE.1) THEN
         
            WMAX = ABS( ENERG(NINIT,LINIT) - E_IONIC )/DBLE(NBPHOT-1) 
         ENDIF
      ENDIF

      WRITE(*,*) '#                    e_init = ', 
     &                                ENERG( NINIT, LINIT ) * enau
      WRITE(*,*) '#             new     w_min = ', WMIN * enau 
      WRITE(*,*) '#             new     w_max = ', WMAX * enau


C e_fin_min

      EN = ENERG( NINIT, LINIT ) + DBLE( NBPHOT ) * WMIN

      CALL FINDFINALSTATE(LEND,NFINMIN,EN,ENERG,NDIMENSION,MAXDIM,MAXL)

      WRITE(*,*) '#                 e_fin_min = ', en
      WRITE(*,*) '#                 n_fin_min = ', nfinmin

C e_fin_max

      EN = ENERG( NINIT, LINIT ) + DBLE( NBPHOT ) * WMAX


      CALL FINDFINALSTATE(LEND,NFINMAX,EN,ENERG,NDIMENSION,MAXDIM,MAXL)

      WRITE(*,*) '#                 e_fin_max = ', en
      WRITE(*,*) '#                 n_fin_max = ', nfinmax

C...................

C         NTOTAL = NFINMAX-NFINMIN


C#      calculate cs for the final state en(nend,lend)


	 DO 9 NEND = NFINMIN, NFINMAX-1


	    D_N = 0.0D+00

C#     photon

            W =( ENERG(NEND, LEND) - ENERG(NINIT, LINIT))/DBLE(NBPHOT)

C#        ABOVE THRESHOLD LEVEL

            ATLEVEL = IDINT(( E_IONIC - ENERG( NINIT, LINIT ) )/ W ) + 1

	    WRITE(14,*)'#                   PHOTON  : ', W
	    WRITE(14,*)'#  FINAL STATE    NEND( W ) : ', NEND
	    WRITE(14,*)'#               ATI PHOTONS : ', NBPHOT - ATLEVEL
C
	    DO 8 CHEMIN = 1, NBPATH
C
           WRITE(44,FMT='(8(1X,I2))')(ARBRELM(CHEMIN,J,1),J=1,NBPHOT+1)
	   WRITE(14,FMT='(8(1X,I3))')(ARBRELM(CHEMIN,J,1),J=1,NBPHOT+1)

C#         ATI case

            DO I = 1, ATLEVEL - 1

               EPSILON( I ) = 0.0D+00
            ENDDO

C#      if    NBPHOTON - ATLEVEL = 0 THIS PART IS NOT EXECUTED

            DO I = ATLEVEL, NBPHOT - 1

               EINT = ENERG( NINIT, LINIT ) + I * W
               LINT = ARBRELM( CHEMIN , I + 1, 1)

               CALL FINDFINALSTATE(LINT,NINT,EINT,ENERG,NDIMENSION,
     S                            MAXDIM, MAXL)

               ENERGYSPACING( I ) = 1.0D+00 / RNORME( NINT, LINT )**2

c               ENERGYSPACING( I ) = 1.0D+00/AN( NINT, LINT )**2

               EPSILONMIN( I ) = EPSFACTMIN * ENERGYSPACING( I )
               EPSILONMAX( I ) = EPSFACTMAX * ENERGYSPACING( I )

               DELTA( I ) = ( EPSILONMAX( I ) - EPSILONMIN( I ) ) / 
     S                      ( NBPTS - 1 )

            ENDDO   
C#
            DO K = NBPTS, 1, -1

               DO I = ATLEVEL, NBPHOT - 1

                  EPSILON( I ) = EPSILONMIN( I ) + ( K - 1 )* DELTA( I )
               ENDDO

               EPS( K ) = EPSILON( ATLEVEL )
               DMX( K )= SIGMA(CHEMIN,W,DIPR,NINIT,NBPHOT,NEND,EPSILON,
     &                    ENERG,RNORME,PHASE,NDIMENSION)
               
C#	       WRITE(44,FMT='(8(1X,I2))')( ARBRELM(N,J,1),J=1,NBPHOT+1)
C#	       WRITE(44,FMT='(8(1X,I2))')( ARBRELM(N,J,2),J=1,NBPHOT+1)

               WRITE(44,FMT='(1X,I4,4(1X,E15.8))') NEND,W,EPS(K),DMX(K)

C#               WRITE(34,FMT='(6(1X,E15.8))') EPSILON(ATLEVEL) / 
C#     S              ENERGYSPACING(ATLEVEL),DMX(K),-DAP(K)

            ENDDO
C#
            DO K = 1, NBPTS

               
               DMX_R( K ) = DREAL( DMX( K ) )
               DMX_I( K ) = DIMAG( DMX( K ) )

            ENDDO

C# INTERPOLATION BY RATIONAL FRACTION

            CALL RATINT(EPS, DMX_R, NBPTS, 0.D0, AMPLIREER, DELTAREER )
            CALL RATINT(EPS, DMX_I, NBPTS, 0.D0, AMPLIIMER, DELTAIMER )

            WRITE(14,FMT='(1X,6(1X,E15.8))') W, AMPLIREER, DELTAIMER 
            WRITE(14,FMT='(1X,6(1X,E15.8))') W, AMPLIIMER, DELTAIMER
C
	    D_N = D_N + DCMPLX( AMPLIREER, AMPLIIMER )

C            WRITE(64,FMT='(1X,6(1X,E15.8))')
C     S            AMPLIREER, AMPLIIMER, AMPLIREER, AMPLIIMER
            
 8       CONTINUE

         WN = W**(NBPHOT)

         IF(MODE.EQ.0) THEN                
            WN = W**(-NBPHOT) 
         ENDIF

C         wn = w**(2*mode-nbphot)

         S_N_AU = CS_N_AU * wn * ABS(D_N)**2
         S_N_SI = CS_N_SI * wn * ABS(D_N)**2

C FOR N = 1 move to Mb

         S_N_SI = S_N_SI * UNIT_N

C #...............................

         CALL CSPHFILE(NCSNPH, LEND, MODE)

C           WRITE(NCSNPH,FMT='(1x,i3,8(1X,1PE15.8))') 
C     &          NEND, 

           WRITE(NCSNPH,FMT='(8(1X,1PE15.8))')          
     &          W * ENAU, 
     &          S_N_SI, 
     &          S_N_AU, 
     &          ABS(D_N)**2, 
     &          DREAL(D_N), 
     &          DIMAG(D_N),  
     &          PHASE(NEND,LEND), 
     &          AN(NEND,LEND)
C#
C#
C#
         WRITE(54,FMT='(3(1X,1PE15.8),2(1X,I4))')
     &                   W, 
     &                   S_N_SI, S_N_AU, 
     &                   NEND, LEND 

 9       CONTINUE

C#  end of state (nend) of the partial wave (lend)
         
         CLOSE(28) 

      ENDDO

C#  end of partial waves (lend)
C#                       
      WRITE(*,*) '#'
      WRITE(*,*) '#'
      WRITE(*,*) '#    DATA FILES STORED IN dat DIRECTORY '
      WRITE(*,*) '#    LOG  FILES STORED IN out DIRECTORY '
      WRITE(*,*) '#'
      WRITE(*,*) '#    END OF CALCULATION'
      WRITE(*,*) '#'
      WRITE(*,*) '###########################################'

      END
Cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   
C#eof
