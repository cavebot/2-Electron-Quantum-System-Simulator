C*  PROGRAM TO PRINT THE SELECTED ENERGY EIGENVALUE AND EIGENVECTOR
C*  OF SPECIFIC EIGENSTATES STORED IN UNFORMATTED DATA FILES.
C*  PROGRAM TO TRANSFER EIGENVALUES AND EIGENVECTORS TO UNFORMATTED
C*  FILE FOR FUTURE CALCULATION

      PROGRAM WF2E
      IMPLICIT REAL*8(B-H,O-Z)
      INCLUDE "parameter.2e.inc"
C     !
      DIMENSION U(nhx)
      DIMENSION LL(ncs),  LHF(ncs) 
      DIMENSION NHF(ncs), NMIN(ncs), NMX(ncs)
      DIMENSION ND(ncs)
      CHARACTER*100 ARGV
Cxxxxxxxxxxxxxx


      CALL GETARG(1, ARGV)
      READ(ARGV,*) L
      CALL GETARG(2, ARGV)
      READ(ARGV,*) ndmin
      CALL GETARG(3, ARGV)
      READ(ARGV,*) ndmax
C............................................

C     #  inpout file

      CALL D2EFILE(9,l)

      READ(9)  LO, LS
      READ(9)  NCSMX
      READ(9)  (NHF(K), K = 1, NCSMX)
      READ(9)  (LHF(K), K = 1, NCSMX)
      READ(9)  (LL(K),  K = 1, NCSMX)
      READ(9)  (NMIN(K),K = 1, NCSMX)
      READ(9)  (NMX(K), K = 1, NCSMX)
      READ(9)  ITRY  
      READ(9)  ( ND(K), K = 1, NCSMX)
      READ(9)  NHMX
      READ(9)  NDTOT
C......

C     #  outpout file

      CALL WF2EFILE(1,l)

      WRITE(1)   LO, LS
      WRITE(1)   NCSMX
      WRITE(1) ( NHF(K),  K = 1, NCSMX )
      WRITE(1) ( LHF(K),  K = 1, NCSMX )
      WRITE(1) ( LL(K),   K = 1, NCSMX )
      WRITE(1) ( NMIN(K), K = 1, NCSMX )
      WRITE(1) ( NMX(K),  K = 1, NCSMX )
      WRITE(1) ( ND(K),   K = 1, NCSMX)
      WRITE(1)   NHMX
      WRITE(1)   NDMAX - NDMIN + 1

C...


      DO 30 NDD = 1, NDTOT

      READ(9) EN
      READ(9) (U(K), K = 1, NHMX)

      IF(NDMIN.EQ.0.AND.NDMAX.EQ.0) GO TO 30
      IF(NDD.LT.NDMIN.OR.NDD.GT.NDMAX) GO TO 30

      WRITE(1)   EN
      WRITE(1) ( U(K), K = 1, NHMX )


 30   CONTINUE

      CLOSE(9)
      CLOSE(1)


C.................................................


    1 FORMAT(2X,1P9D14.6)
    3 FORMAT(3I5)
    4 FORMAT(/2X,'DATA SET # =',I3)
    5 FORMAT(/2X,'EIGEN-VECTOR --'/)
    6 FORMAT(/2X,'L AND LS=',2I3/)
    8 FORMAT(/2X,'EIGENVALUE=',1PE15.6)
   10 FORMAT(2X,'TOTAL # OF DATA SETS IN THIS FILE=',I3/)
   12 FORMAT(/2X,'IF EIGENVALUE & EIGENVECTOR ARE TO BE PRINTED, ENTER'
     1  /2X,'"1" IN any FORMAT, ENTER OTHER INTEGER IF OTHERWISE'/)
   13 FORMAT(2X,'TOTAL # OF EIGENSTATES IN THIS DATA SET=',I3//2X,
     1  'ENTER THE RANGE OF EIGENSTATES STORE IN UNFORMATTED FILE'/2X,
     2  '  NDMIN & NDMAX=? -- any format'/)
   14 FORMAT(2X,'ENTER THE RANGE OF EIGENSTATES TO BE PRINTED'/2X,
     1  'NPMIN & NPMAX=? -- any format'/)
   15 FORMAT(2X,'TO SKIP THIS STEP, ENTER "0" AND "0" IN any format'/)

C.....

      END

C####################################################################
C**   DATA ON "DIAGxxxx.DAT" IS ARRANGED IN THE FOLLOWING ORDER
C**   L AND LS -- ORBITAL AND SPIN ANGULAR MOMENTA
C**   NCSMX -- # OF CONFIFURATION SERIES
C**   NHF,LHF,LL,NMIN,NMX -- ALL FROM 1 TO NCSMX
C**   ITRY - # OF DATA SETS IN THE FILE
C**   ND(K) - K FROM 1 TO NCSMX; # OF CONFIG. IN EACH SERIES
C**   NHMX - TOTAL # OF CONFIGURATIONS IN THE EIGENVECTORS
C**   NDTOT - # OF STATES WITH OUTPUT IN THE FILE, REPEATING 
C**      ENERGY EIGENVALUE
C**      ENERGY EIGENVECTOR  (1 TO NHMX)
C*********************************************************************
C--   DATA ON "EGVxxxx.DAT" IS STORED IN THE FOLLOWING ORDER
C--   NCSMX
C--   NHF,LHF,LL,NMIN,NMX
C--   ND
C--   NHMX
C--   NS -- # OF EIGENVECTORS
C--      REPEATING  NS TIMES
C--         ENERGY EIGENVALUE
C--         ENERGY EIGENVECTOR (TOTAL OF NHMX NUMBERS EACH)
C#####################################################################


      SUBROUTINE EGVPRNT
      IMPLICIT REAL*8(B-H,O-Z)
      DIMENSION U(2100),NHF(152),LHF(152),LL(152),NMIN(152),NMX(152),
     1    ND(152)

C............................................

      READ(9)  LO, LS
      READ(9)  NCSMX
      READ(9)  (NHF(K), K = 1, NCSMX)
      READ(9)  (LHF(K), K = 1, NCSMX)
      READ(9)  (LL(K),  K = 1, NCSMX)
      READ(9)  (NMIN(K),K = 1, NCSMX)
      READ(9)  (NMX(K), K = 1, NCSMX)
      READ(9) ITRY


C............................................

      WRITE(1) LO,LS
      WRITE(1) NCSMX
      WRITE(1) (NHF(K), K=1,NCSMX)
      WRITE(1) (LHF(K), K=1,NCSMX)
      WRITE(1) (LL(K), K=1,NCSMX)
      WRITE(1) (NMIN(K), K=1,NCSMX)
      WRITE(1) (NMX(K), K=1,NCSMX)
      WRITE(*,6) LO, LS
      WRITE(*,10) ITRY
C............................................

      DO 100 IT = 1, ITRY

      WRITE(*,4) IT

  
      READ(9) (ND(K), K=1,NCSMX)
C#
C#    print out the configuration 
C#


C*    put 1 or !=1 from keyboard     

      WRITE(*,12)

      READ(9) NHMX
      READ(9) NDTOT

      READ(*,*) IYES

      print*, 'IYES',iyes

      IF(IYES .NE. 1) GO TO 25

      WRITE(*,13) NDTOT
      WRITE(*,15)

      READ(*,*) NDMIN,NDMAX

      WRITE(*,14)
      WRITE(*,15)

C*    Give the range of eigenstates  NPMIN  up to NPMAX to be printed

      READ(*,*) NPMIN,NPMAX
      WRITE(6,*) (ND(K), K=1,NCSMX)

      WRITE(1) (ND(K), K=1,NCSMX)
      WRITE(1) NHMX
      NS = NDMAX - NDMIN + 1
      WRITE(1) NS

 25   CONTINUE


      DO 30 NDD = 1, NDTOT

      READ(9) EN
      READ(9) (U(K), K = 1,NHMX)

      IF(IYES .NE. 1) GO TO 30

      IF(NDMIN .EQ. 0 .AND. NDMAX .EQ. 0) GO TO 35

      IF(NDD .LT. NDMIN .OR. NDD .GT. NDMAX) GO TO 35

      WRITE(1) EN
      WRITE(1) (U(K), K = 1, NHMX)

   35 IF(NPMIN .EQ. 0 .AND. NPMAX .EQ. 0) GO TO 30
      IF(NDD .LT. NPMIN .OR. NDD .GT. NPMAX) GO TO 30

      WRITE(*,*)  EN
      WRITE(*,5)
      WRITE(*,1) (U(K),  K = 1, NHMX)


   30 CONTINUE

  100 CONTINUE


C---------------------------------------------------------------------
C**   DATA ON "DIAGxxxx.DAT" IS ARRANGED IN THE FOLLOWING ORDER
C**   L AND LS -- ORBITAL AND SPIN ANGULAR MOMENTA
C**   NCSMX -- # OF CONFIFURATION SERIES
C**   NHF,LHF,LL,NMIN,NMX -- ALL FROM 1 TO NCSMX
C**   ITRY - # OF DATA SETS IN THE FILE
C**   ND(K) - K FROM 1 TO NCSMX; # OF CONFIG. IN EACH SERIES
C**   NHMX - TOTAL # OF CONFIGURATIONS IN THE EIGENVECTORS
C**   NDTOT - # OF STATES WITH OUTPUT IN THE FILE, REPEATING 
C**      ENERGY EIGENVALUE
C**      ENERGY EIGENVECTOR  (1 TO NHMX)
C*********************************************************************
C--   DATA ON "EGVxxxx.DAT" IS STORED IN THE FOLLOWING ORDER
C--   NCSMX
C--   NHF,LHF,LL,NMIN,NMX
C--   ND
C--   NHMX
C--   NS -- # OF EIGENVECTORS
C--      REPEATING  NS TIMES
C--         ENERGY EIGENVALUE
C--         ENERGY EIGENVECTOR (TOTAL OF NHMX NUMBERS EACH)
C.......................................................................

    1 FORMAT(2X,1P9D14.6)
    3 FORMAT(3I5)
    4 FORMAT(/2X,'DATA SET # =',I3)
    5 FORMAT(/2X,'EIGEN-VECTOR --'/)
    6 FORMAT(/2X,'L AND LS=',2I3/)
    8 FORMAT(/2X,'EIGENVALUE=',1PE15.6)
   10 FORMAT(2X,'TOTAL # OF DATA SETS IN THIS FILE=',I3/)
   12 FORMAT(/2X,'IF EIGENVALUE & EIGENVECTOR ARE TO BE PRINTED, ENTER'
     1  /2X,'"1" IN any FORMAT, ENTER OTHER INTEGER IF OTHERWISE'/)
   13 FORMAT(2X,'TOTAL # OF EIGENSTATES IN THIS DATA SET=',I3//2X,
     1  'ENTER THE RANGE OF EIGENSTATES STORE IN UNFORMATTED FILE'/2X,
     2  '  NDMIN & NDMAX=? -- any format'/)
   14 FORMAT(2X,'ENTER THE RANGE OF EIGENSTATES TO BE PRINTED'/2X,
     1  'NPMIN & NPMAX=? -- any format'/)
   15 FORMAT(2X,'TO SKIP THIS STEP, ENTER "0" AND "0" IN any format'/)
C ---------------------------------------------------------------

      RETURN
      END

C#######################################################################
      SUBROUTINE DIAGS(NHF,LHF,LL,NMIN,NMX,ND,NCSMX)
      DIMENSION NHF(1),LHF(1),LL(1),NMIN(1),NMX(1),ND(1)
    1 FORMAT(/2X,'# OF CONF. SERIES=',I4/2X,
     1  'NHF,LHF,LL,NMIN,NMAX,NDIAG --'/)
    2 FORMAT(25I3)
      WRITE(*,*)  NCSMX
      WRITE(*,2) ( NHF(K),  K = 1, NCSMX)
      WRITE(*,2) ( LHF(K),  K = 1, NCSMX)
      WRITE(*,2) ( LL(K),   K = 1, NCSMX)
      WRITE(*,2) ( NMIN(K), K = 1, NCSMX)
      WRITE(*,2) ( NMX(K),  K = 1, NCSMX)
      WRITE(*,2) ( ND(K),   K = 1, NCSMX)
      RETURN
      END

C#######################################################################
      SUBROUTINE DIAGW(NHF,LHF,LL,NMIN,NMX,ND,NCSMX)
      DIMENSION NHF(1),LHF(1),LL(1),NMIN(1),NMX(1),ND(1)
    1 FORMAT(/2X,'# OF CONF. SERIES=',I4/2X,
     1  'NHF,LHF,LL,NMIN,NMAX,NDIAG --'/)
    2 FORMAT(25I3)
      WRITE(16,*) NCSMX
      WRITE(16,2) ( NHF(K),  K = 1, NCSMX)
      WRITE(16,2) ( LHF(K),  K = 1, NCSMX)
      WRITE(16,2) ( LL(K),   K = 1, NCSMX)
      WRITE(16,2) ( NMIN(K), K = 1, NCSMX)
      WRITE(16,2) ( NMX(K),  K = 1, NCSMX)
      WRITE(16,2) ( ND(K),   K = 1, NCSMX)
      RETURN
      END

C#######################################################################
C#    SUBROUTINE TO OPEN THREE FILES
C#    DIAGxx.DAT INPUT UNFORMATTED FILE
C#    DIAGxx.EGV OUTPUT FORMATTED FILE INCLUDES EIGENVALUE & VECTORS
C#    EGVxx.DAT  OUTPUT UNFORMATTED FILES.
C#######################################################################
      SUBROUTINE OPENFL

c........................................................................
    1 FORMAT(2X,'ENTER THE TOL'/2X,
     1 ' FOR SINGLET-S, TRIPLET-S, SINGLET-P, TRIPLET-P, SINGLET-D,'
     2 /2X,'     TRIPLET-D, SINGLET-F, AND TRIPLET-F RESPECTIVELY'/)
 3    format(5A16)

C.....................................

      open(7,file='inp/wf2e.inp',status='old')

      read(7,*) l

      close(7)
C.....................................

      
      CALL D2EFILE(9,l)  
      CALL WF2EFILE(1,l)



      RETURN
      END
C#######################################################################
C#EOF
