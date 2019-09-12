C#######################################################################
      SUBROUTINE TRNSMX(HMX,XM,IHR,IHC,IR,IC)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      INCLUDE "parameter.2e.inc"
C      PARAMETER(NHX=2500,NS=802)
      DIMENSION HMX(NHX,NHX),XM(NS-2,NS-2)
      IRF=IHR+IR-1
      ICF=IHC+IC-1
      DO 10 K=IHR,IRF
      DO 10 KK=IHC,ICF
      KR=K-IHR+1
      KC=KK-IHC+1
      HMX(K,KK)=XM(KR,KC)
  10   CONTINUE
      RETURN
      END
C#######################################################################
      SUBROUTINE ANGFC(LLE,LE,LLI,LI,LEX,LIN,ID,ANG)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION ANG(4),ID(4)

!
! l1= lle=ljp
! l3= lli=lip
!
      CALL ANGF(LLE,LE,LLI,LI,LEX,LIN,ID(1),ANG(1)) ! < j'j |D| i' i > (1)
      CALL ANGF(LLE,LE,LI,LLI,LEX,LIN,ID(2),ANG(2)) ! < j'j |D| ii' >  (2)
      CALL ANGF(LE,LLE,LLI,LI,LEX,LIN,ID(3),ANG(3)) ! < jj' |D| i'i >  (3)
      CALL ANGF(LE,LLE,LI,LLI,LEX,LIN,ID(4),ANG(4)) ! < jj' |D| ii' >  (4)

      RETURN
      END
C#######################################################################
      SUBROUTINE ANGF(L1,L2,L3,L4,LE,LI,ID,ANG)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)


C*     LS-coupled states  ---> |JM_J LS >

      ANG = ANGLS(L1, L2, L3, L4, LE, LI)

C*     LS-uncoupled states --->  |LS M_L M_S >


c      ANG = ANOS(L1, L2, L3, L4, LE, LI)

      ID = 0
      IF(ANG.EQ.0.0D0) RETURN

      CALL IDMX(L3, L1, ID)         ! 

      RETURN
      END
C#######################################################################
      SUBROUTINE IDMX(LI,LE,ID)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

      IT = IABS(LI-LE)
      IF(IT .NE. 1) STOP
      IF(LI .EQ. 0 .AND. LE .EQ. 1) ID = 1         ! <li=0|D|le=1> = <lip|t|ljp>
      IF(LI .EQ. 1 .AND. LE .EQ. 0) ID = 2         ! <li=1|D|le=0>
      IF(LI .EQ. 1 .AND. LE .EQ. 2) ID = 3         ! <li=1|D|le=2>
      IF(LI .EQ. 2 .AND. LE .EQ. 1) ID = 4         ! <li=2|D|le=1>
      IF(LI .EQ. 2 .AND. LE .EQ. 3) ID = 5
      IF(LI .EQ. 3 .AND. LE .EQ. 2) ID = 6
      IF(LI .EQ. 3 .AND. LE .EQ. 4) ID = 7   
      IF(LI .EQ. 4 .AND. LE .EQ. 3) ID = 8
      IF(LI .EQ. 4 .AND. LE .EQ. 5) ID = 9
      IF(LI .EQ. 5 .AND. LE .EQ. 4) ID = 10
      IF(LI .EQ. 5 .AND. LE .EQ. 6) ID = 11
      IF(LI .EQ. 6 .AND. LE .EQ. 5) ID = 12
      IF(LI .EQ. 6 .AND. LE .EQ. 7) ID = 13
      IF(LI .EQ. 7 .AND. LE .EQ. 6) ID = 14
      IF(LI .EQ. 7 .AND. LE .EQ. 8) ID = 15
      IF(LI .EQ. 8 .AND. LE .EQ. 7) ID = 16
      IF(LI .EQ. 8 .AND. LE .EQ. 9) ID = 17
      IF(LI .EQ. 9 .AND. LE .EQ. 8) ID = 18
      IF(LI .EQ. 9 .AND. LE .EQ. 10)ID = 19
      IF(LI .EQ. 10 .AND. LE .EQ. 9)ID = 20
      IF(LI .EQ. 10 .AND. LE .EQ. 11)ID =21
      RETURN
      END
C#######################################################################
c ls coupling (l,s,j|r|l',s,j')
      FUNCTION ANGLS(L11,L12,L21,L22,L1,L2)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      COMMON/ANG/J1,J2,LS
C*
C*    check first if triangular rule of 3j/6j symbols is fulfiled
C*            
C*    if   not  then     angls = 0
C*    else go on to the calculation  
C*

      IF(L12.EQ.L22) GO TO 10
      ANGLS = 0.D0
      RETURN

   10 IT = IABS(J1 - J2)
      IF(IT .EQ. 1) GO TO 20
      ANGLS = 0.D0
      RETURN

   20 IT = IABS(L11 - L21)
      IF(IT .EQ. 1) GO TO 30
      ANGLS = 0.D0
      RETURN

   30 IT = J1 + J2 + L1 + L2 + LS + L22

      A = ( 2.D0*J1 + 1.D0) * (2.D0*J2 +  1.D0)
     1                      * (2.D0*L1 +  1.D0) * (2.D0*L2  + 1.D0)
     1                      * (2.D0*L11 + 1.D0) * (2.D0*L21 + 1.D0)

      ANGLS = (-1)**IT * DSQRT(A)
     1                 * F3J( 2*J1, 2, 2*J2, 0, 0, 0) 
     1                 * F3J( 2*L11,2, 2*L21, 0, 0, 0) 
     1                 * F6J( 2*L1, 2*J1, 2*LS, 2*J2, 2*L2, 2) 
     1                 * F6J( 2*L11,2*L1, 2*L12,2*L2, 2*L21, 2)

      RETURN
      END
C#######################################################################
      FUNCTION ANOS(LLE,LE,LLI,LI,LEX,LIN)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)

C*
C*    check first if triangular rule of 3j/6j symbols is fulfiled
C*            
C*    if   not  then     angls = 0
C*    else go on to the calculation  
C*

      IF(LE .EQ. LI) GO TO 10
      ANOS = 0.D0
      RETURN

   10 IT = IABS(LLE-LLI)

      IF(IT .EQ. 1) GO TO 20
      ANOS = 0.D0
      RETURN

   20 JJE = 2 * LLE

      JE  = 2 * LE
      JJI = 2 * LLI
      JI  = 2 * LI
      JEX = 2 * LEX
      JIN = 2 * LIN

      A = DFLOAT(JJE + 1) * DFLOAT(JJI + 1)

      ANOS = (-1)**IABS(LLE) * DSQRT(A) 
     1                       * F3J(JJE,2,JJI,0,0,0) 
     1                       * F6J(JIN,2,JEX,JJE,JI,JJI)

      RETURN
      END
C#######################################################################
      FUNCTION F6J(JD1,JD2,JD3,LD1,LD2,LD3)
C  F6J FUNCTION CALLS S6J  FORTRAN IV
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION MED(12)
      J1=JD1
      J2=JD2
      J3=JD3
      L1=LD1
      L2=LD2
      L3=LD3
C  ANGULAR MOMENTUM COUPLING TESTS FOR 6J COEFICINETS
      I=-J1+J2+J3
      I1=I/2
      IF (I-2*I1)1000,1010,1000
 1000 F6J=0.0D0
      GOTO 100
 1010 MED(1)=I1
      I=J1-J2+J3
      I1=I/2
      IF (I-2*I1) 1000,1020,1000
 1020 MED(2)=I1
      I=+J1+J2-J3
      I1=I/2
      IF (I-2*I1) 1000,1030,1000
 1030 MED(3)=I1
      I=-J1+L2+L3
      I1=I/2
      IF (I-2*I1) 1000,1040,1000
 1040 MED(4)=I1
      I=+J1-L2+L3
      I1=I/2
      IF (I-2*I1) 1000,1050,1000
 1050 MED(5)=I1
      I=+J1+L2-L3
      I1=I/2
      IF (I-2*I1) 1000,1060,1000
 1060 MED(6)=I1
      I=-L1+J2+L3
      I1=I/2
      IF (I-2*I1) 1000,1070,1000
 1070 MED(7)=I1
      I=+L1-J2+L3
      I1=I/2
      IF (I-2*I1) 1000,1080,1000
 1080 MED(8)=I1
      I=+L1+J2-L3
      I1=I/2
      IF (I-2*I1) 1000,1090,1000
 1090 MED(9)=I1
      I=-L1+L2+J3
      I1=I/2
      IF (I-2*I1) 1000,1100,1000
 1100 MED(10)=I1
      I=+L1-L2+J3
      I1=I/2
      IF (I-2*I1) 1000,1110,1000
 1110 MED(11)=I1
      I=+L1+L2-J3
      I1=I/2
      IF (I-2*I1) 1000,1120,1000
 1120 MED(12)=I1
      DO 10 N=1,12
      IF (MED(N)) 1000,10,10
   10 CONTINUE
      F6J=S6J(J1,J2,J3,L1,L2,L3)
  100 RETURN
      END
C#######################################################################
      DOUBLE PRECISION FUNCTION S6J(JD1,JD2,JD3,LD1,LD2,LD3)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION FL(322),MA(4),MB(3),MED(12)
      DFLOAT(I)=I
      J1=JD1
      J2=JD2
      J3=JD3
      L1=LD1
      L2=LD2
      L3=LD3
      FL(1)=0.0D0
      FL(2)=0.0D0
      DO 50 N=3,322
      FN=DFLOAT(N-1)
   50 FL(N)=FL(N-1)+DLOG(FN)
   15 MED(1)=(-J1+J2+J3)/2
      MED(2)=(+J1-J2+J3)/2
      MED(3)=(+J1+J2-J3)/2
      MED(4)=(-J1+L2+L3)/2
      MED(5)=(+J1-L2+L3)/2
      MED(6)=(+J1+L2-L3)/2
      MED(7)=(-L1+J2+L3)/2
      MED(8)=(+L1-J2+L3)/2
      MED(9)=(+L1+J2-L3)/2
      MED(10)=(-L1+L2+J3)/2
      MED(11)=(+L1-L2+J3)/2
      MED(12)=(+L1+L2-J3)/2
      MA(1)=MED(1)+MED(2)+MED(3)
      MA(2)=MED(4)+MED(5)+MED(6)
      MA(3)=MED(7)+MED(8)+MED(9)
      MA(4)=MED(10)+MED(11)+MED(12)
      MB(1)=MA(1)+MED(12)
      MB(2)=MA(1)+MED(4)
      MB(3)=MA(1)+MED(8)
C  DETERMINE MAXMUM OF (J1+J2+J3), (J1+L2+L3),(L1+J2+L3),(L1+L2+J3)
      MAX=MA(1)
      DO 30 N=2,4
      IF (MAX-MA(N)) 20,30,30
   20 MAX=MA(N)
   30 CONTINUE
C  DETERMINE MINIMUM OF (J1+J2+L1+L2), (J2+J3+L2+L3),(J3+J1+L3+L1)
      MIN=MB(1)
      DO 51 N=2,3
      IF (MIN-MB(N)) 51,51,40
   40 MIN=MB(N)
   51 CONTINUE
      KMAX=MIN-MAX
      MINP1=MIN+1
      MINI=MINP1-MA(1)
      MIN2=MINP1-MA(2)
      MIN3=MINP1-MA(3)
      MIN4=MINP1-MA(4)
      MIN5=MINP1+1
      MIN6=MB(1)-MIN
      MIN7=MB(2)-MIN
      MIN8=MB(3)-MIN
C  SUM SERIES IN DOUBLE PRECISION
      UK=1.D-15
      S=1.D-15
      IF (KMAX)65,65,55
   55 DO 60 K=1,KMAX
      UK=-UK*DFLOAT(MINI-K)*DFLOAT(MIN2-K)*DFLOAT(MIN3-K)*DFLOAT(MIN4-K
     1 )/(DFLOAT(MIN5-K)*DFLOAT(MIN6+K)*DFLOAT(MIN7+K)*DFLOAT(MIN8+K))
C  CUT OFF SERIES AT 1.0D-25
      IF (DABS(UK)-1.D-25) 65,65,60
   60 S=S+UK
   65 S=S*1.0D+15
C  CALCULATE DELTA FUNCTIONS
      DELOG=0.D0
      DO 70 N=1,12
      NUM=MED(N)
   70 DELOG=DELOG+FL(NUM+1)
      NUM1=MA(1)+2
      NUM2=MA(2)+2
      NUM3=MA(3)+2
      NUM4=MA(4)+2
      DELOG=DELOG-FL(NUM1)-FL(NUM2)-FL(NUM3)-FL(NUM4)
      DELOG=0.5D0*DELOG
      ULOG=FL(MIN5)-FL(MINI)-FL(MIN2)-FL(MIN3)-FL(MIN4)-FL(MIN6+1)-
     * FL(MIN7+1)-FL(MIN8+1)
      PLOG=DELOG+ULOG
      IF (SNGL(PLOG)+64.0D0) 72,75,75
   72 Q=PLOG+64.D0
      Q=DEXP(Q)
      S6J=Q*S
      IF (DABS(S6J)-1.D0) 73,73,74
   73 S6J=0.D0
      GOTO 90
   74 S6J=S6J*DEXP(-64.D0)
      GOTO 78
   75 P=DEXP(PLOG)
      S6J=P*S
   78 MIN2=MIN/2
      IF (MIN-2*MIN2) 80,90,80
   80 S6J=-S6J
   90 CONTINUE
      RETURN 
      END
C*********************************************************************
      FUNCTION F3J(JD1,JD2,JD3,MD1,MD2,MD3)
C  F3J VERSION II
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION FL(322),MTRI(9)
      DFLOAT(I)=I
      J1=JD1
      J2=JD2
      J3=JD3
      M1=MD1
      M2=MD2
      M3=MD3
      FL(1)=0.D0
      FL(2)=0.D0
      DO 50 N=3,322
      FN=N-1
   50 FL(N)=FL(N-1)+DLOG(FN)
   15 I=J1+J2-J3
      I1=I/2
      IF (I-2*I1)1000,1010,1000
 1010 MTRI(1)=I1
      I=J1-J2+J3
      I1=I/2
      IF (I-2*I1) 1000,1020,1000
 1020 MTRI(2)=I1
      I=-J1+J2+J3
      I1=I/2
      IF (I-2*I1) 1000,1030,1000
 1030 MTRI(3)=I1
      IF (M1+M2+M3)1000,1040,1000
 1040 I=J1+M1
      I1=I/2
      IF (I-2*I1) 1000,1050,1000
 1050 MTRI(4)=I1
      MTRI(5)=(J1-M1)/2
      I=J2+M2
      I1=I/2
      IF (I-2*I1) 1000,1060,1000
 1060 MTRI(6)=I1
      MTRI(7)=(J2-M2)/2
      I=J3+M3
      I1=I/2
      IF (I-2*I1) 1000,1070,1000
 1070 MTRI(8)=I1
      MTRI(9)=(J3-M3)/2
      DO 30 N=1,9
      IF (MTRI(N)) 1000,30,30
   30 CONTINUE
      IF (J3-J2+M1) 40,45,45
   40 KMIN=-J3+J2-M1
      GOTO 60
   45 KMIN=0
   60 IF (-J3+J1+M2-KMIN) 80,80,70
   70 KMIN=-J3+J1+M2
   80 KMIN=KMIN/2
      IF(J2-J3+M1) 90,100,100
   90 KMAX=J1+J2-J3
      GOTO 110
  100 KMAX=J1-M1
  110 IF (J2+M2-KMAX) 120,130,130
  120 KMAX=J2+M2
  130 KMAX=KMAX/2
      MINI=MTRI(1)-KMIN+1
      MIN2=MTRI(5)-KMIN+1
      MIN3=MTRI(6)-KMIN+1
      MIN4=(J3-J2+M1)/2+KMIN
      MIN5=(J3-J1-M2)/2+KMIN
C   SUM SERIES IN DOUBLE PRECISION
      UK=1.D-10
      S=1.D-10
      NCUT=0
      KMAX=KMAX-KMIN
      IF (KMAX)165,165,155
  155 DO 160 K=1,KMAX
      UK=-UK*DFLOAT(MINI-K)*DFLOAT(MIN2-K)*DFLOAT(MIN3-K)/(DFLOAT(KMIN+
     * K)*DFLOAT(MIN4+K)*DFLOAT(MIN5+K))
      IF (DABS(UK)-1.D30) 158,157,157
  157 UK=1.D-10*UK
      S=1.D-10*S
      NCUT=NCUT+1
  158 IF (DABS(UK)-1.D-20) 165,160,160
  160 S=S+UK
C   CALCULATE DELTA FUNCTIONS
  165 DELOG=0.D0
      DO 170 N=1,9
      NUM=MTRI(N)
  170 DELOG=DELOG+FL(NUM+1)
      NUM=(J1+J2+J3)/2+2
      DELOG=0.5D0*(DELOG-FL(NUM))
      ULOG=-FL(KMIN+1)-FL(MINI)-FL(MIN2)-FL(MIN3)-FL(MIN4+1)-FL(MIN5+1)
      PLOG=DELOG+ULOG
      IF (SNGL(PLOG)+80.0D0) 172,171,171
  171 IF (NCUT)175,175,172
  172 SIG=DSIGN(1.D0,S)
      S=DABS(S)
      SLOG=DLOG(S)+DFLOAT(NCUT+1)*DLOG(1.D+10)
      F3J=SIG*DEXP(SLOG+PLOG)
      GOTO 178
  175 S=S*1.D+10
      P=DEXP(PLOG)
      F3J=P*S
  178 NUM=KMIN+(J1-J2-M3)/2
      IF (MOD(NUM,2)) 180,190,180
  180 F3J=-F3J
  190 CONTINUE
      GOTO 2000
 1000 F3J=0.0D0
 2000 RETURN
      END
C#######################################################################
C#   EOF
