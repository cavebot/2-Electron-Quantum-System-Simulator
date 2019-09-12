      SUBROUTINE YKFCTB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NO,IDR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(1),PP(1),YK(1),ZK(1),WY(1),WZ(1),R(1),DR(1),
     1  DRKR(1),RKK(1)

      IF(IDR.EQ.0) CALL YKSINB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NO)

      IF(IDR.EQ.1) CALL YKEXPB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NO)

      IF(IDR.GT.1) THEN

         WRITE(*,*) ' IDR > 1. EXECUTION WILL TERMINATE. IDR = 0 OR 1  ' 

         STOP
   
      ENDIF

      RETURN
      END
C -----------------------------------------------------------------
C --- YK FUNCTION IN  r = ri * exp( (x+xi)**(1/nrp) )  SCALE
C --- x = (n-1) * H
C --- DRKR = 1.0/R**KR
C --- RKK  = R**KK
C -
      SUBROUTINE YKEXPB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NO)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(1),PP(1),YK(1),ZK(1),WY(1),WZ(1),R(1),DR(1),
     1  DRKR(1),RKK(1)
      N=1
      KR=K+N-1
      KK=K+1
      IF(KR .EQ. 0) GO TO 30
      DO 10 J=1,NO
      WZ(J)=P(J)*PP(J)*DRKR(J)/DR(J)
   10 WY(J)=P(J)*PP(J)*RKK(J)/DR(J)
      CALL YINTB(WY,WZ,YK,ZK,NO,H)
      DO 20 J=1,NO
   20 YK(J)=YK(J)*DRKR(J)+ZK(J)*RKK(J)
      RETURN
   30 DO 40 J=1,NO
      WZ(J)=P(J)*PP(J)/DR(J)
   40 WY(J)=P(J)*PP(J)*RKK(J)/DR(J)
      CALL YINTB(WY,WZ,YK,ZK,NO,H)
      DO 50 J=1,NO
   50 YK(J)=YK(J)+ZK(J)*RKK(J)
      RETURN
      END
c ------------------------------------------------------------------
 
C --- DRKR = 1.0/R**KR
C --- RKK  = R**KK
C -
      SUBROUTINE YKSINB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NO)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(1),PP(1),YK(1),ZK(1),WY(1),WZ(1),R(1),DR(1),
     1  DRKR(1),RKK(1)
      N=1
      KR=K+N-1
      KK=K+1
      WZ(1)=0.0D0
      WY(1)=0.0D0
      IF(KR .EQ. 0) GO TO 30
      DO 10 J=2,NO
      WZ(J)=P(J)*PP(J)*DRKR(J)/DR(J)
   10 WY(J)=P(J)*PP(J)*RKK(J)/DR(J)
      CALL YINTB(WY,WZ,YK,ZK,NO,H)
      YK(1)=0.0D0
      DO 20 J=2,NO
   20 YK(J)=YK(J)*DRKR(J)+ZK(J)*RKK(J)
      RETURN
   30 DO 40 J=2,NO
      WZ(J)=P(J)*PP(J)/DR(J)
   40 WY(J)=P(J)*PP(J)*RKK(J)/DR(J)
      CALL YINTB(WY,WZ,YK,ZK,NO,H)
      YK(1)=0.0D0
      DO 50 J=2,NO
   50 YK(J)=YK(J)+ZK(J)*RKK(J)
      RETURN
      END
c ------------------------------------------------------------------
      subroutine yintb (v,w,y,z,m,h)
c
c  this program calculates the indefinite integrals y and z using the
c  lagrange integration formula
c
c  y(r) = integral of v from 0 to r
c  z(r) = integral of w from r to infinity
c  m is the maximum tabulation point of v and w (virtual infinity)
c  h is the step size of the radial grid
c
      IMPLICIT REAL*8(A-H,O-Z)
      dimension v(1),w(1),y(1),z(1)
c ------------------------------------------------------------------
c                lagrange 10 point integration formula
c               --------------------------------------
      dimension aa(5,10),a(10,5),b(5)
c     equivalence (b(1),a(6,5))
      data ia/5/, ja/10/, da/7257600.D0/
      data aa/2082753.D0, -57281.D0, 10625.D0, -3969.D0, 2497.D0,
     &   9449717.D0, 2655563.D0, -163531.D0,   50315.D0, -28939.D0,
     & -11271304.D0, 6872072.D0, 3133688.D0, -342136.D0, 162680.D0,
     &  16002320.D0,-4397584.D0, 5597072.D0, 3609968.D0,-641776.D0,
     & -17283646.D0, 3973310.D0,-2166334.D0, 4763582.D0,4134338.D0,
     &  13510082.D0,-2848834.D0, 1295810.D0,-1166146.D0,4134338.D0,
     &  -7394032.D0, 1481072.D0, -617584.D0,  462320.D0,-641776.D0,
     &   2687864.D0, -520312.D0,  206072.D0, -141304.D0, 162680.D0,
     &   -583435.D0,  110219.D0,  -42187.D0,   27467.D0, -28939.D0,
     &     57281.D0,  -10625.D0,    3969.D0,   -2497.D0,   2497.D0/
      data b/ 4134338.D0, -641776.D0, 162680.D0, -28939.D0, 2497.D0/
      data h0/0.0D0/
c
c  note that a different even order method can be used by replacing the
c  dimension and data statements in this block (also in do loop 20)
c ------------------------------------------------------------------
      if(h.eq.h0) go to 5
      hd=h/da
      do 2 i=1,ia
         do 1 j=1,ja
            a(j,i)=aa(i,j)*hd
 1       continue
         b(i) = b(i)*hd
 2    continue
      h0=h
 5    y(1)=0.0D0
      z(m)=0.0D0
      do 10 i=2,ia
      k=m-i+1
      y(i)=y(i-1)
      z(k)=z(k+1)
      ii=i-1
      do 10 j=1,ja
      y(i)=y(i)+a(j,ii)*v(j)
 10   z(k)=z(k)+a(j,ii)*w(m-j+1)
      im=ia+1
      in=m-ia+1
      do 20 i=im,in
      k=m-i+1
      y(i)=y(i-1)+b(1)*(v(i  )+v(i-1))+b(2)*(v(i+1)+v(i-2))
     &           +b(3)*(v(i+2)+v(i-3))+b(4)*(v(i+3)+v(i-4))
     &           +b(5)*(v(i+4)+v(i-5))
 20   z(k)=z(k+1)+b(1)*(w(k  )+w(k+1))+b(2)*(w(k-1)+w(k+2))
     &           +b(3)*(w(k-2)+w(k+3))+b(4)*(w(k-3)+w(k+4))
     &           +b(5)*(w(k-4)+w(k+5))
      in=in+1
      do 30 i=in,m
      k=m-i+1
      y(i)=y(i-1)
      z(k)=z(k+1)
      do 30 j=1,ja
      y(i)=y(i)+a(j,k)*v(m-j+1)
 30   z(k)=z(k)+a(j,k)*w(j)
      return
      end
C########################################




