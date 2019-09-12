!###################################################################
      SUBROUTINE CXFIN(NCFG, NHF,LHF,LL,NMIN,NMX,NOL,IS,
     &                       L, LS, NCSMX, NCMX, NDI, IDCS)
      PARAMETER(NSX=130) 
      DIMENSION NHF(NSX),LHF(NSX),LL(NSX),NMIN(NSX),NMX(NSX),NOL(NSX)
      DIMENSION IS(NSX),IDCS(NSX),NDI(NSX)
!..................................................................
 3    FORMAT(6I5)

      CALL HCFGFILE(NCFG, L)
      READ( NCFG, 3 ) NSYMMETRIES
      READ(NCFG, 3) L, LS
      READ(NCFG, 3) NCSMX
      NCMX = 0
      DO   K = 1, NCSMX
        READ(NCFG, 3) NHF(K),LHF(K),LL(K),NMIN(K),NMX(K),IDCS(K)
        NOL(K) = NMX(K) - NMIN(K) + 1
        IS(K)  = NCMX + 1 
        NCMX   = NCMX + NOL(K)         
      ENDDO

      CLOSE(NCFG)
      
C#check
C#
C#     NOL :  TOTAL NUMBER OF B-SPLINES
C#    IS   :  GLOBAL POSITION OF THE VECTOR OF BASIS CHANNELS
C#    NCMX :  NUMBER OF B - SPLINES FOR CHANNEL   
C#                                              | NHF, LHF, LL >
C#
C      WRITE(*,*) NSYMMETRIES, L, LS, NCSMX
C      DO   K = 1, NCSMX
C     WRITE(*,*) NHF(K), LHF(K), LL(K), NMIN(K), NMX(K), IDCS(K)
C      enddo
C      READ(NCFG, 3) ITRY
C      READ(NCFG, 7) ( NDI(K), K = 1, NCSMX)

      DO K = 1, NCSMX
         NDI(K) = NOL(K) 
         IF(NDI(K).NE.NOL(K)) THEN 
            WRITE(*,*) '# cxfin::        error in the cfg file:'
            WRITE(*,*) '# cxfin::                          K = ', K 
            WRITE(*,*) '# cxfin::                        NDI = ', NDI(K)
            WRITE(*,*) '# cxfin::                        NOL = ', NOL(K)
         ENDIF

CWRITE(*,*) K, NDI(K), NOL(K) 

       ENDDO

      WRITE(*,*)  '# cxfin::                           L  = ', L
      WRITE(*,*)  '# cxfin::                           S  = ', LS
      WRITE(*,*)  '# cxfin::                          NCS = ', NCSMX
      WRITE(*,*)  '# cxfin::                          NCF = ', NCMX

C      READ(NCFG, *) (ND(K), K = 1, NCSMX)

      IF(LS.EQ.0.OR.LS.EQ.1) THEN 
         LS = 2*LS + 1 
      ELSE
         LS = LS
      ENDIF

      RETURN
      END
!!!###################################################################
      SUBROUTINE GETNUMBEROFCHANNELS(NCFG, L, NDMAX, NBSP)
      IMPLICIT NONE
      INTEGER NCFG,NSYMMETRIES,L,LS,NDMAX
      INTEGER NBSP
      INTEGER NTMP
!......................................................
 3    FORMAT(6I5)

      IF(NBSP.EQ.-1) THEN
         CALL CFGFILE(NCFG, L)
      ELSE
         CALL HCFGFILE(NCFG, L)
      ENDIF

      READ( NCFG, 3 ) NSYMMETRIES, ntmp
      write(*,*) nsymmetries,ntmp
      READ( NCFG, 3 ) L,LS 
      write(*,*) l,ls
      READ( NCFG, 3 ) NDMAX
      write(*,*) ndmax
      READ( NCFG, 3 ) NTMP, NTMP, NTMP, NTMP, NBSP, NTMP

      
      CLOSE(NCFG)

C.......................................................
      RETURN
      END
!!!###################################################################
!!!#  EOF
