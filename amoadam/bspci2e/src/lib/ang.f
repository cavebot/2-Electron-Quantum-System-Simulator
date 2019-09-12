C*******************************************************************
C#   THIS SUB. CALC ANGULAR PART OF COULOUMB MATRIX
C#
C#   KNO                     :: NUMBER OF THE ALLOWED K VALUE; 
C#   KA(I), I = 1, KNO
C#   FK(I), I = 1, KNO       ::  ANGULAR FACTOR
C#
C#
      SUBROUTINE ANGFK( L1, L2, L3, L4, L, F, KA, KNO)
      IMPLICIT REAL*8(A-H),REAL*8(O-Z)
      DIMENSION F(1),KA(1)

      LS1 = L1 + L3
      LS3 = L2 + L4
      MAX = MIN0(LS1,LS3) + 1

      LS2 = IABS(L1 - L3)
      LS4 = IABS(L2 - L4)
      MIN = MAX0(LS2, LS4) + 1

      C   = FLOAT( (2*L1+1) * (2*L2+1) * (2*L3+1) * (2*L4+1) )

      J1  = 2 * L1
      J2  = 2 * L2
      J3  = 2 * L3
      J4  = 2 * L4
      J   = 2 * L

      KNO = 0
      DO 10 K = MIN, MAX, 2

      K1  = 2*( K - 1)

      A1  = F3J(J1, J3, K1, 0, 0, 0)          
      A2  = F3J(J2, J4, K1, 0, 0, 0)
      B   = F6J(J1, J2, J, J4, J3, K1)

      KNO = KNO + 1
      IP  = IABS( L1 - L3 + L )

      F(KNO) = (-1)**IP * DSQRT(C) * A1 * A2 * B

   10 CONTINUE

      MI = MIN - 1

      DO 20 I = 1, KNO
 20      KA(I) = MI + (I-1)*2

      RETURN
      END
C*******************************************************************
C***  FUNCTION TO CALCULATE THE 6-J COEFICIENTS
C***  F6J FUNCTION CALLS S6J  FORTRAN IV
C  ANGULAR MOMENTUM COUPLING TESTS FOR 6J COEFFICIENTS

      FUNCTION F6J(J1,J2,J3,L1,L2,L3)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ANGF/FL(322),MA(4),MB(3),MED(12),MTRI(9)
C
      I=-J1+J2+J3
      I1=I/2
      IF (I-2*I1)1000,1010,1000
 1000 F6J=0.0
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
C  *****************************************************************
      DOUBLE PRECISION FUNCTION S6J(J1,J2,J3,L1,L2,L3)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ANGF/FL(322),MA(4),MB(3),MED(12),MTRI(9)
      DFLOAT(I)=I
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
      IF (SNGL(PLOG)+64.0) 72,75,75
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
C#######################################################################
C***  FUNCTION TO CALCULATE THE 3-J COEFICIENTS
C***  F3J VERSION II
      FUNCTION F3J(J1,J2,J3,M1,M2,M3)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ANGF/FL(322),MA(4),MB(3),MED(12),MTRI(9)
      DFLOAT(I)=I
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
      IF (SNGL(PLOG)+80.0) 172,171,171
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
 1000 F3J=0.0
 2000 RETURN
      END
C#######################################################################
C#EOF
