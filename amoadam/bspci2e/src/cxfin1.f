!###################################################################
      SUBROUTINE CXFIN(NCFG, NHF,LHF,LL,NMIN,NMX,NOL,IS,
     &                       L, LS, NCSMX, NCMX, NDI, IDCS)
      PARAMETER(NSX=130) 
      DIMENSION NHF(NSX),LHF(NSX),LL(NSX),NMIN(NSX),NMX(NSX),NOL(NSX)
      DIMENSION IS(NSX),IDCS(NSX),NDI(NSX)
      !
 3    FORMAT(6I5)


      CALL cfile(ncfg,"inp","CFG",l) 

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
C      READ(NCFG, 3) ITRY
C      READ(NCFG, 7) ( NDI(K), K = 1, NCSMX)

      DO K = 1, NCSMX
         NDI(K) = NOL(K) 
         IF(NDI(K).NE.NOL(K)) THEN 
            WRITE(*,'(a40,i10)') '& cxfin::     error in the cfg file:'
            WRITE(*,'(a40,i10)') '& cxfin::               K = ', K 
            WRITE(*,'(a40,i10)') '& cxfin::             NDI = ', NDI(K)
            WRITE(*,'(a40,i10)') '& cxfin::             NOL = ', NOL(K)
         ENDIF
       ENDDO

       WRITE(*,'(a40,i10)') '& cxfiin::  nof channels ncs = ',ncsmx
       WRITE(*,'(a40,i10)') '# cxfin::           nof cfgs = ', ncmx


       RETURN
       END
!!!###################################################################
      SUBROUTINE GETNUMBEROFCHANNELS(NCFG, L, NDMAX, NBSP)
      IMPLICIT NONE
      INTEGER:: NCFG, NSYMMETRIES, L, LS, NDMAX
      INTEGER:: NBSP
      INTEGER:: NTMP
!......................................................
 3    FORMAT(5I5)
 4    FORMAT(6I5)
!

      IF(NBSP.EQ.-1) THEN
         CALL cfile(ncfg,"inp","cfg",l)      !         CALL CFGFILE(NCFG, L)
      ELSE
         CALL cfile(ncfg,"inp","CFG",l)      !         CALL HCFGFILE(NCFG, L)
      ENDIF

      READ( NCFG, 3 ) NSYMMETRIES
      READ( NCFG, 3 ) L,LS 
      READ( NCFG, 3 ) NDMAX

      if(NBSP.EQ.-1) THEN
         READ( NCFG, 3 ) NTMP, NTMP, NTMP, NTMP, NBSP
      ELSE
         READ( NCFG, 4 ) NTMP, NTMP, NTMP, NTMP, NBSP, NTMP
      ENDIF
      
      WRITE(*,'(a22,10x,a28,1X,i10)') '# getnumberofchannels::',
     &' nof channels ncs = ',ndmax
      WRITE(*,'(a22,10x,a28,1X,i10)') '# getnumberofchannels::',
     &' nof basis functions nb = ', nbsp



      
!      write(*,*) nsymmetries,ntmp
!      write(*,*) l,ls
!      write(*,*) ndmax,nbsp

      close(ncfg)

      RETURN
      END
!!!###################################################################
!!!#  EOF
