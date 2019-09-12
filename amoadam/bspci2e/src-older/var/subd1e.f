c**  geteig   === End of Compilation 1 ===
c**  corewf   === End of Compilation 2 ===
c**  vhfbsp   === End of Compilation 3 ===
c**  rkprep   === End of Compilation 4 ===
c**  directb  === End of Compilation 5 ===
c**  xchb     === End of Compilation 6 ===
c**  ykfctb   === End of Compilation 7 ===
c**  ykexpb   === End of Compilation 8 ===
c**  yksinb   === End of Compilation 9 ===
c**  soleig   === End of Compilation 10 ===
c**  setmat   === End of Compilation 11 ===
c**  q        === End of Compilation 12 ===
c**  hwf      === End of Compilation 13 ===
C.......................................................................

      subroutine geteig(lang,n,k,dr0,rmax,no,nos,nop,ntc,nm,h,
     1 t,r,dr,p,vd,ff,yk,zk,wy,wz,br,bc,amat,bmat,cmat,ad,ae,er,
     2 redmass, fv1, fv2, work, ogp, idr)
C..........
      implicit doubleprecision(a-h,o-z)
      INCLUDE "parameter.1e.inc"
      dimension t(1),r(1),dr(1),p(nl, np),vd(1),ff(1),yk(1),zk(1),
     1 wy(1),wz(1),br(1),bc(1),amat(nm,nm),bmat(nm,nm),cmat(nm,nm),
     1 ad(nm,nm),ae(nm,nm),er(1),fv1(1),fv2(1),work(1),ogp(1)

C.........  determine  knots              t(i)  i=1,n+k


      call mkgrid(n,k,rmax,dr0,t)
      call vhfbsp(t,r,dr,p,vd,ff,yk,zk,wy,wz,br,bc,ad,ae,h,n,k,no,lang,
     1     nos,nop,ntc,idr)
      call setmat(nm,n,k,t,amat,bmat,ad,ae,redmass)
      call soleig(nm,n-2,amat,bmat,cmat,er,fv1,fv2,work,ogp,iflag)

      if(iflag.ne.0) go to 999

      return
 999  stop
      end
c*********************************************************
c** print a matrix
      subroutine printmat(bm,nm,n)
      double precision bm(nm, nm)
      integer i
 8    format(10f12.3)
      write(*,*) ' print matrix' 
      do  i = 1, n
         write(*,8) ( bm(j,i), j =1,n)
      enddo
      return
      end
c***************************************************************
      subroutine corewf(c,p,f,r,dr,coe,t,n,k,ni,nf,nx,no,ntc,h)
      implicit doubleprecision(a-h,o-z)
      INCLUDE "parameter.1e.inc" 
      dimension c(nx,nx),p(nl,np),f(np),r(np),dr(np)
      dimension coe(nx),t(1)
C........
      in = n - 1
      coe(1) = 0.0D+00
      coe(n) = 0.0D+00

      do 20 npp = ni, nf

        in = in - 1

        do 10 i = 1, n - 2
 10        coe( i + 1 ) = c(i, in)

        p(npp, 1 ) = 0.0D+00

        do 12 j = 2, no
 12        p(npp,j) = bvalue(t, coe, n, k, r(j), 0 )

C..... ?
        if(p(npp, 6).gt.0.0D+00) goto 14

        do 13 j = 1, no
 13        p(npp,j) = - p(npp,j)

 14     continue

        f(1) = 0.0D+00

        do 15 j = 2, no
 15        f(j) = p(npp,j) * p(npp,j) * r(j) / dr(j)

           oth = rint(f,1,no,14,h)

c        write(*,*)  ' core wf # ', np, '  norm =', oth
c        write(16,*) ' core wf # ', np, '  norm =', oth

        sq2 = sqrt(oth)

        do 16 j = 1, no

 16        p(npp,j) = p(npp,j) / sq2

 20   continue

      return
      end
C#################################################################
      SUBROUTINE vhfbsp(t,r,dr,p,vd,ff,yk,zk,wy,wz,br,bc,ad,ae,h,n,k,
     1     no,lang,nos,nop,ntc,idr)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      INCLUDE "parameter.1e.inc" 
C      PARAMETER(NS=802,NP=2000, NL=15)
      DIMENSION T(1),R(1),DR(1),p(nl,np),VD(1),FF(1),C(NS),
     1 YK(1),zk(1),wy(1),wz(1),BR(1),BC(1),AD(NS,NS),AE(NS,NS)
      COMMON/BSPF/BF(NS-2,NP)

    4 FORMAT(1X,1P9E14.6)

C........................

      call directb(p,vd,yk,zk,wy,wz,dr,r,ff,h,no,nos,nop,
     1     ntc,idr)

C.......................


      DO 50 NC = 1, N - 2


         DO  NN = 1, N

            C(NN)   = 0.0D+00
         ENDDO


            C( NC + 1 ) = 1.0D0

            DO  J = 1, NO

               BF(NC, J) = BVALUE( T, C, N, K, R(J), 0)

            ENDDO

 50         CONTINUE

C......................

      DO 100 KR = 1, N - 2

         DO  J = 1, NO

            BR(J) = BF( KR, J)
         ENDDO

            call xchb( br,p,ff,yk,zk,wy,wz,lang,no,dr,nos,nop,ntc,h,r,
     1           idr)


 10   DO 101 KC = 1, KR

         DO  J = 1, NO
            BC(J) = BF(KC, J )
         ENDDO

         if(dr(1).gt.1.0D-15) FF(1) = -BR(1) * BC(1) * VD(1) / DR(1)

            DO 61 J = 2, NO
 61            FF(J) = - BR(J) * BC(J) * VD(J) / dr(j)

      AD(KR,KC) = RINT(FF, 1, NO, 14, H)

C      AE(KR,KC) = 0.0D+00

      AD(KC, KR) = AD(KR,KC)

      DO  J = 2, NO
         FF(J) = YK(J) * BC(J) / DR(J)
C         write(*,'(1X,I4,1x,E20.14)') j, ff(j)
      ENDDO

      AE(KR,KC) = RINT( FF, 1, NO, 14, H)

C      do j = 1, 10 
C         write(*,'(1X,2E20.14)') AD(j,j+1), AE(j,j+1) 
C      enddo

      AE(KC,KR) = AE(KR,KC)

 101  CONTINUE
 100  CONTINUE

      RETURN
      END
c-------------------------------------------------------------------
** This routine prepares r**k and 1/r**(k+1) for use in routine YKFCTB.
      subroutine rkprep(r,drkr,rkk,k,no)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      INCLUDE "parameter.1e.inc"
C      PARAMETER(NP=2000, NL=15) 

C.........................................
      COMMON/RARR/RPW(NL+3,NP)
      DIMENSION DRKR(1),RKK(1),r(1)

    
      KR = K 
      KK = K + 1

      IF(KR.EQ.0) GO TO 30

      DRKR(1) = 0.0D+00

      IF(R(1).NE.0.0D0) DRKR(1) = 1.0D+00 / RPW(KR, 1 )

      DO 20 J = 2, NO

   20 DRKR(J) = 1.0D0 / RPW(KR,J)

   30 CONTINUE

      DO 21 J = 1, NO
   21  RKK(J)   =  RPW(KK,J)

      RETURN

      END
C........................................
C..
C..      Direct Potential   Y(r) ~ 2*Y/r
C..
C........................................
      subroutine directb(p,ft,yk,zk,wy,wz,dr,r,f,h,no,nos,nop,ntc,
     1 idr)
      implicit doubleprecision(a-h,o-z)
      INCLUDE "parameter.1e.inc"
C      PARAMETER(np=2000, nl = 15)
      dimension p(nl,np),ft(1),yk(1),zk(1),wy(1),wz(1),dr(1),r(1),
     1 f(1),DRKR(NP),RKK(NP)

c*    (1s^2)(2s^2)(2p^6)(3s^2)(3p^6)(3d^10)

      qS = 2.0d0
      qP = 6.0d0
      qD = 10.0d0

c*
c*   calculate   r^(k+1) , 1/r^k at grid points
c*  

      CALL RKPREP(R, DRKR, RKK, 0, NO)


c********************************************************
c**    no core  He, H- 
c** 
c********************************************************

      do 10 j  = 1, no
         f(j)  = p(1,j)
 10      ft(j) = 0.0D+00

      if(nos.eq.0) return

c********************************************************
c*
c*    ( 1s^2 )  core 
c*
c********************************************************

c*
c*      s = 0   core states
c*

      call ykfctb(f,f,yk,zk,wy,wz,r,dr,DRKR,RKK,h,0,no,idr)

      do 11 j  = 1, no
 11      ft(j) = ft(j) + 2.d0 * qS * yk(j)

      if(nos.eq.1) return
c********************************************************
c*
c*    ( 1s^2 2s^2 2p^6)  core 
c*
c********************************************************

c*
c*      s = 0   core states
c*

      do 12 j = 1, no
 12      f(j) = p(2,j)

      call ykfctb(f,f,yk,zk,wy,wz,r,dr,DRKR,RKK,h,0,no,idr)

      do 13 j = 1, no
         f(j)  = p( nos + 1, j)
 13      ft(j) = ft(j) + 2.d0 * qS * yk(j)


c*
c*      p = 1   core states
c*

      call ykfctb(f,f,yk,zk,wy,wz,r,dr,DRKR,RKK,h,0,no,idr)

      do 14 j = 1, no
 14      ft(j)= ft(j) + 2.d0 * qP *yk(j)

      if(nos.eq.2) return

c********************************************************
c*
c*    ( 1s^2 2s^2 2p^6 3s^2 3p^6 3d^10)  core 
c*
c********************************************************
      do 15 j = 1,no
 15      f(j) = p(3,j)

c*
c*      s = 0   core states
c*

      call ykfctb(f,f,yk,zk,wy,wz,r,dr,DRKR,RKK,h,0,no,idr)

      do 16 j = 1, no
         f(j) =  p(nos+2,j)
 16      ft(j) = ft(j) + 2.d0 * qS *yk(j)


c*
c*      p = 1   core states
c*

      call ykfctb(f,f,yk,zk,wy,wz,r,dr,DRKR,RKK,h,0,no,idr)

      do 17 j = 1, no
 17      ft(j) = ft(j) + 2.d0 * qP * yk(j)

      if(nos.eq.3) return

c*
c*      s = 0   core states
c*


      do 18 j = 1,no
 18      f(j) = p(4,j)

      call ykfctb(f,f,yk,zk,wy,wz,r,dr,DRKR,RKK,h,0,no,idr)

      do 19 j = 1, no
         f(j)  = p(nos+3,j)
 19      ft(j) = ft(j) + 2.d0* qS * yk(j)

c*
c*      p = 1   core states
c*

      call ykfctb(f,f,yk,zk,wy,wz,r,dr,DRKR,RKK,h,0,no,idr)


      do 20 j=1,no
         f(j) = p(nos+nop+1,j)
 20      ft(j) = ft(j) + 2.d0 * qP * yk(j)


c*
c*      d = 2   core states
c*

      call ykfctb(f,f,yk,zk,wy,wz,r,dr,DRKR,RKK,h,0,no,idr)

      do 21 j = 1, no
 21      ft(j) = ft(j) + 2.d0* qD *yk(j)

      if(nos.eq.4) return

      end

C#############################################################
C#      Exchange potential   X(r)
C#
C#############################################################
      subroutine xchb(pr,p,f,yk,zk,wy,wz,l,no,dr,nos,nop,ntc,h,r,
     1     idr)
c---
      implicit doubleprecision(a-h,o-z)
      INCLUDE "parameter.1e.inc"
C      PARAMETER(np=2000, nl = 15)
      dimension pr(1),p(nl,np),f(1),dr(1),r(1),yk(1),y1(np),
     1     ZK(1),WY(1),WZ(1),DRKR(NP),RKK(NP)
c----
      FS   = 2.0D0 / DFLOAT(2*L+1)
      FPM  = 3.d0  * FS  * DFLOAT(L)/ DFLOAT(2*L-1)
      FPP  = 3.d0  * FS  * DFLOAT(L+1)/ DFLOAT(2*L+3)
      FDM  = 2.5d0 * FPM * DFLOAT(L-1)/ DFLOAT(2*L-3)
      FD   = 5.d0  * FPM * DFLOAT(L+1)/ DFLOAT(2*L+3)/3.d0
      FDP  = 2.5d0 * FPP * DFLOAT(L+2)/ DFLOAT(2*L+5)
      LPM  = L - 1
      LPP  = L + 1
      LDM  = L - 2
      LDP  = L + 2

      if(nos.eq.0) return

      do  j = 1, no
         f(j) = p(1,j)
         yk(j) = 0.d0
      enddo

c*
c*   calculate   r^(k+1) , 1/r^k at grid points
c*  

      CALL RKPREP(R,DRKR,RKK,L,no)

      
      CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,L,NO,IDR)


      if(nos.eq.0)  return

      do  j = 1, no
         yk(j) = yk(j) + fs * y1(j) * f(j)
      enddo

      if(nos.eq.1) return

      do  j = 1, no
         f(j) = p(2,j)
      enddo

      CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,L,NO,IDR)

      do  j = 1, no
         f(j) = p(nos+1,j)
         yk(j) = yk(j) + fs * y1(j) * p(2,j)
      enddo

      CALL RKPREP(R,DRKR,RKK,LPP,no)

      CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPP,NO,IDR)

      do  j = 1, no
         yk(j) = yk(j) + fpp * y1(j) * f(j)
      enddo

       if(lpm.ge.0) then

          CALL RKPREP(R,DRKR,RKK,LPM,no)

          CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPM,NO,IDR)

          do  j = 1,no
             yk(j) = yk(j) + fpm * y1(j) * f(j)
          enddo

       endif

      if(nos.eq.2) return


      do  j = 1, no
         f(j) = p(3,j)
      enddo

      CALL RKPREP(R,DRKR,RKK,L,no)

      CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,L,NO,IDR)

      do  j = 1, no
         f(j) = p(nos + 2,j)
         yk(j) = yk(j) + fs * y1(j) * p(3,j)
      enddo

      CALL RKPREP(R,DRKR,RKK,LPP,no)
         
      CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPP,NO,IDR)

      do  j = 1, no
         yk(j) = yk(j) + fpp * y1(j) * f(j)
      enddo

       if(lpm.ge.0) then

          CALL RKPREP(R,DRKR,RKK,LPM,no)

          CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPM,NO,IDR)

          do j = 1, no
             yk(j) = yk(j) + fpm * y1(j) * f(j)
          enddo

       endif

      if(nos.eq.3) return

      do  j = 1, no
         f(j) = p(4,j)
      enddo

      CALL RKPREP(R,DRKR,RKK,L,no)

      CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,L,NO,IDR)

      do j = 1, no
         f(j) = p(nos+3,j)
         yk(j) = yk(j) + fs * y1(j) * p(4,j)
      enddo
      
      CALL RKPREP(R,DRKR,RKK,LPP,no)

      CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPP,NO,IDR)

      do  j = 1,no
         yk(j) = yk(j) + fpp * y1(j) * f(j)
      enddo

      if(lpm.ge.0) then

         CALL RKPREP(R,DRKR,RKK,LPM,no)

         CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LPM,NO,IDR)

          do  j = 1, no
             yk(j) = yk(j) + fpm * y1(j) * f(j)
          enddo

       endif

       do  j = 1,no
          f(j) = p(nos + nop + 1, j)
       enddo
 
       CALL RKPREP(R,DRKR,RKK,LDP,no)

       CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LDP,NO,IDR)

       do  j = 1, no
          yk(j) = yk(j) + fdp * y1(j) * f(j)
       enddo
       
       CALL RKPREP(R,DRKR,RKK,L,no)
       
       CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,L,NO,IDR)

       do  j = 1, no
          yk(j) = yk(j) + fd * y1(j) * f(j)
       enddo

       if(ldm.ge.0) then

          CALL RKPREP(R,DRKR,RKK,LDM,no)

          CALL YKFCTB(Pr,f,Y1,ZK,WY,WZ,R,DR,DRKR,RKK,H,LDM,NO,IDR)

          do  j = 1, no
             yk(j) = yk(j) + fdm * y1(j) * f(j)
          enddo

       endif

       return
       end
c*****************************************************************

      SUBROUTINE YKFCTB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NO,
     1     IDR)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION P(1),PP(1),YK(1),ZK(1),WY(1),WZ(1),R(1),DR(1),
     1  DRKR(1),RKK(1)

      IF(IDR.EQ.0) CALL YKSINB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NO)

      IF(IDR.EQ.1) CALL YKEXPB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NO)

      IF(IDR.GT.1) STOP

      RETURN
      END

c -----------------------------------------------------------------
c --- YK FUNCTION IN  r = ri * exp( (x+xi)**(1/nrp) )  SCALE
c --- x = (n-1) * H
c --- DRKR = 1.0/R**KR
c --- RKK  = R**KK
c ---
      SUBROUTINE YKEXPB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NO)

      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION P(1),PP(1),YK(1),ZK(1),WY(1),WZ(1),R(1),DR(1),
     1  DRKR(1),RKK(1)

      N = 1

      KR = K + N - 1

      KK = K + 1

      IF(KR.EQ.0) GO TO 30

      DO 10 J = 2, NO

      WZ(J)=P(J)*PP(J)*DRKR(J)/DR(J)
   10 WY(J)=P(J)*PP(J)*RKK(J)/DR(J)

      CALL YINT(WY,WZ,YK,ZK,NO,H)

      DO 20 J=2,NO
   20 YK(J)=YK(J)*DRKR(J)+ZK(J)*RKK(J)
      RETURN

   30 DO 40 J=2,NO
      WZ(J)=P(J)*PP(J)/DR(J)
   40 WY(J)=P(J)*PP(J)*RKK(J)/DR(J)

      CALL YINT(WY,WZ,YK,ZK,NO,H)

      DO 50 J=2,NO
   50 YK(J)=YK(J)+ZK(J)*RKK(J)

      RETURN

      END
c*********************************************************************
c*
c*   Initial version :
c*   This routine performs the Y^k(r ;i,j) integration :
c*
c*   relation (4.5), pg 48, ch. 4 of M. Amusia & L. Chernysheva book 
c*   'Computation of Atomic Processes'
c*                        
c* 
c*                                                             |r
c*   Y^k(r; i,j) == ( 1 / r^k)  *  < Pi(r') | (r')^k | Pj(r') >|
c*                                                             |0
c*
c*                                                             |oo
c*               + r^(k+1) * < Pi(r') | (1/r')^(k+1) | Pj(r') >|
c*                                                             |r
c*      
c*   Y^k(r;i,j)  =  Y_1 / r^k   + Y_2 * r^(k+1)
c*   drkr = 1.0/r^k   
c*   rkk  = r^(k+1)
c*   dr   = r * dx/dR   (see subroutine rin for definition of dx,dR
c*
c*                                     l.aa.n  9/07/2000
c --- DRKR = 1.0/R**KR
C --- RKK  = R**KK
c*********************************************************************
      SUBROUTINE YKSINB(P,PP,YK,ZK,WY,WZ,R,DR,DRKR,RKK,H,K,NO)
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION P(1),PP(1),YK(1),ZK(1),WY(1),WZ(1),R(1),DR(1),
     1  DRKR(1),RKK(1)
c*
c*    N is not used. Why not removed ?
c*    kr = k always
c*

      KR    = K 
      KK    = K + 1
      WZ(1) = 0.0D0
      WY(1) = 0.0D0

c*
c*    for k = 0  no integration for Y_1 is needed
c*

      IF(KR .EQ. 0) GO TO 30

c*
c*   prepare integral argument Y_1, Y_2
c*       
c*     Y_1  ~    Pi(r) * Pj(r) * r^k        
c*       
c*     Y_2  ~    Pi(r) * Pj(r) / r^(k+1)           
c*

      DO  J = 2, NO

         WY(J) = P(J) * PP(J) * RKK(J) / DR(J)
         WZ(J) = P(J) * PP(J) * DRKR(J) / DR(J)

      ENDDO

c*
c*    integrate now 
c*

      CALL YINT(WY, WZ, YK, ZK, NO, H)


      YK(1)   =  0.0D+00

         DO J = 2, NO

           YK(J) = YK(J) * DRKR(J) + ZK(J) * RKK(J)  

         ENDDO

      RETURN
c*
c*    if k = 0  r^(k+1) = r    and 1/r^k = 1    
c*

   30  DO  J = 2, NO

          WY(J)  = P(J) * PP(J) * RKK(J) / DR(J)
          WZ(J) =  P(J) * PP(J) / DR(J)

      ENDDO
c*    integrate now

      CALL YINT(WY,WZ,YK,ZK,NO,H)



c   old version
c      DO 50 J=2,NO
c   50 YK(J)= YK(J) + ZK(J) * RKK(J)


      YK(1) = 0.0D0
         DO J = 2,NO

           YK(J) = YK(J) + ZK(J) * RKK(J) 

         ENDDO


      RETURN
      END
***********************************************************
*
*     call eispack to solve general eigenvalue problem
*
*     A c = eps B c
*
***********************************************************
      subroutine soleig(nm,n,a,b,c,er,fv1,fv2,work,ogp,iflag)
      implicit doubleprecision(a-h,o-z)
      dimension a(nm,1),b(nm,1),c(nm,1),er(nm),work(nm,nm),
     * fv1(1),fv2(1),ogp(1)
      data matz /1/

    1 FORMAT(2X,'en(',I2,') & norm =',1P2E19.10)
*  preserve b-matrix


      do j = 1, n
         do i = 1, n
            work(i,j) = b(i,j)
        enddo
      enddo

c      print*,'a'
c      call printmat(a,nm,n) 
c      print*,'b'
c      call printmat(b,nm,n) 

C*     now diagonalize the B-splines hamiltonian matrix


      call rsg(nm,n,a,work,er,matz,c,fv1,fv2,ierr)


C      do  iew = n, 1, -1 
C         write(*,*) n + 1 - iew, er(iew), (c(j,iew), j = 1, 3)  
C      enddo

   
c     call rgg(nm,n,a,work,er,ei,beta,matz,c,ierr)
  
C      DO NB = 1, n

C         OGP(NB)  = 0.D+00

C         DO  IR = 1, n
C            DO  IC = 1, n
C               OGP(NB) = OGP(NB) + C(IR,NB) * B(IR,IC) * C(IC,NB)
C            ENDDO
C         ENDDO
C      ENDDO

      if(ierr.ne.0) go to 999

      do  iew = 1, n
         er(iew) = - er(iew)
      enddo

      iflag = ierr
      return

 901  write(6,*) 'complex eigenvalue encountered'

      iflag = 999
      return

 999  write(6,*) 'matrix eigenvalue routine fails'

      iflag = ierr
      return

      end

**********************************************
*
*      b (i,j) = Int[ Bi(x) Bj(x) dx ]
*
*      b(j,i) = b(i,j)
*      
*      for i<=j  :
*                       t(i+k)
*               b(i,j) =  Int [Bi(x) Bj(x) dx]  j <i+k
*                         t(j)
*               b(i,j) =           0            j>=i+k
*
*      therefore
*                        l=i+k-1  t(l+1)
*               b(i,j)  =  Sum   { Int [ Bi(x) Bj(x) dx ] }
*                          l=j     t(l)
*
*       where t(l), l = 1,n+k = knot sequence
*
*
*     t(l+1)                             k
*      Int [ f(x) dx ]  = [t(l+1)-t(l)] Sum [ w(m) f(x(m)) ]
*     t(l)                              m=1
*
*      x(m) = [t(l+1)-t(l)] z(m) + t(l) 
*      z(m) = gaussian k-point coordinates for [0,1]
*      w(m) = gaussian k-point weights     for [0,1]
*      
******************************************************************

      subroutine setmat(nm,n,k,t,a,b,ad,ae,redmass)
      implicit doubleprecision(a-h,o-z)
      INCLUDE "parameter.1e.inc"
      dimension t(1),a(nm,nm),b(nm,nm),ad(nm,nm),ae(nm,nm)
      dimension db(NK,NK)
      dimension xg(NK),wg(NK)
      data nderiv/2/

*  initialize  a and b
      do  j = 1, n-2
         do   i = 1, n-2
           if(i.eq.j.and.i.le.2) then
              write(*,*) '  V_d(', i, ',', j, ' )       = ', ad(i,j) 
              write(*,*) '  V_e(', i, ',', j, ' )       = ', ae(i,j)                
           endif
C     a(i,j) = ad(i,j) 
           a(i,j) = ad(i,j) + ae(i,j)
           b(i,j) = 0.0D+00

          enddo
        enddo

      call gauss(k,xg,wg)

*  set up non zero loops  i = 2 ... n-1  ,j = i,i+k-1
      do 500 i = 2, n - 1
         jhi = min0(i+k-1,n-1)
         do 400 j = i,jhi
            low = max0(k,j)
            lhi = min0(i+k-1,n)
*        sum over intermediate segments
            do 300 l = low,lhi
               seg = 0.D0
               seh = 0.D0
               dl = (t(l+1)-t(l))
               x0 = t(l)
c           sum over gaussian weights
               do 200 m = 1, k

                  xm = dl*xg(m) + x0

                  call bsplvd(t,k,xm,l,db,nderiv)

                  fm = db(i-l+k,1)*db(j-l+k,1)

                  gm = q(xm)*fm

                  hm = db(i-l+k,2)*db(j-l+k,2)

                  seg = seg + wg(m)*fm

C*                 redmass = 1./mass_reduced  (1==hydrogen)/(2==positronium)

                  seh = seh + wg(m)* redmass * ( gm - hm )

 200           continue

               a(i-1,j-1) = a(i-1,j-1) + dl*seh
               b(i-1,j-1) = b(i-1,j-1) + dl*seg

 300        continue
 400     continue
 500  continue

C* get the remainder
      do  i = 1, n - 2
         do  j = 1,i           
           b(i,j) = b(j,i)
           a(i,j) = a(j,i)
         enddo
       enddo
       
       return
       end
c*****************************************************************
c***  coulomb + centrifugal + core polarization terms

      function q(x)

C..............................
      INCLUDE "parameter.1e.inc"
      implicit doubleprecision(a-h,o-z)
      common/params/alp1(nk),r01(nk),redmass,znuc,lang

C..............................

      QL = DFLOAT(lang*(lang+1))

      IF(alp1(lang+1).EQ.0.0D+00) THEN

         q  = (2.0D+00 * ( znuc/redmass )  - QL/x)/x

      ELSE

         IF(r01(lang+1).EQ.0.0D+00) THEN

            WRITE(*,*)'CUT-OFF PARAMETER EQUALS ZERO R01 = ',r01(lang+1)  
         ENDIF

         q  = (2.0D+00 * ( znuc/redmass )  - QL/x)/x
         q  = q + alp1(lang+1)*(1.0D+00-dexp(-(x/r01(lang+1))**6))/x**4

      ENDIF


      return
      end
c**  coulomb + core polarization terms
C#####################################################################
      FUNCTION HWF(N,L,Z,R,E)
      DOUBLE PRECISION HWF,Z,R,E,FN,RK,RM,RLL,FK,FM,FLL,A,P,X
      M=N+L
      K=N-L-1
      LL=2*L+1
      FN=N
      RK=K
      RM=M
      RLL=LL
      FK=1.D0
      FM=1.D0
      FLL=1.D0
      P=1.D0
      A=1.D0
      X=-2.D0*Z*R/FN
      IF (DABS(X) .GT. 160.D0) GO TO 20
      DO 15 I=1,M
      FM=FM*RM
   15 RM=RM-1.D0
      RM=M
      DO 16 I=1,LL
      FLL=FLL*RLL
   16 RLL=RLL-1.D0
      IF (K) 1,2,3
    3 DO 4 I=1,K
      P=1.D0+A/RK*P/RM*X
      FK=FK*RK
      A=A+1.D0
      RK=RK-1.D0
    4 RM=RM-1.D0

    2 HWF = DSQRT( Z * FM / FK) * P * DEXP(X/2.D0) *(-X)**(L+1)/FLL/FN

      E = - Z**2 / FN**2

      RETURN
   20 HWF=0.d0
      RETURN
    1 WRITE(16,10) N,L,Z,R
c      WRITE(*,10) N,L,Z,R
   10 FORMAT(52H  FORBIDDEN COMBINATION OF N AND L IN HWF SUBPROGRAM/
     17H     N=,I3,5H   L=,I3,5H   Z=,F6.3)
      STOP
      END
C#####################################################################
*********  set up matrix for differential eqn *****************
*
*         y"(x) + q(x) y(x) = eps y(x)
*         y(a) = 0
*         y(b) = 0
*
*                 n
*         y(x) = Sum [ Bj(x) * cj ]
*                j=1
*                  
*         require x=a to be a k-fold knot ==> B1(a) = 1
*                 x=b to be a k-fold knot ==> Bn(b) = 1
*
*         therefore must set c1 = cn = 0 to satisy b.c.
*
*         write d.e. as 
*
*         n-1                                    n-1
*         Sum [ <Bi|Bj"> + <Bi|q|Bj> ] cj  = eps Sum <Bi|Bj> cj
*         j=2                                    j=2
*
*        a(i,j) = <Bi|Bj"> + <Bi|q|Bj>               ;   i = 2,n-1,j=2,n-1
*
*        n.b.     <Bi|Bj"> = - <Bi'|Bj'> = <Bi"|Bj>  ;  i,j != 1,n
*
*        b(i,j) = <Bi|Bj>
*
*        A * c = eps B * c  symmetric (band) matrix eigenvalue problem
*                           of dimension n-2 for vector c
*
*        <y|y> =  c~ B c
*        
*        B is a metrix in c space
*
********************************************************************
C#EOF





