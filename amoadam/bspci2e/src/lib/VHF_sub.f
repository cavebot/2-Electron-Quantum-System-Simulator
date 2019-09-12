c------------------------------------------------------------------ 
      subroutine vhfbsp(t,r,dr,p,ad,ae,h,n,k,no,lang,nos,nop,ntc,idr)
      implicit real*8(a-h,o-z)
      include "parameter.inc"
C      parameter(ns=62, np=2000, nl = 15, kx = 15)
      dimension t(ns+kx),r(np),dr(np),p(nl,np),c(ns),vd(np), yk(np)
     1 ,br(np),bc(np),ad(ns,ns),ae(ns,ns), ff(np),bf(ns-2,np)
    4 format(1x,1p9e14.6)

C.....................................

C....        Direct potential VD


      call directb( p, vd, dr, r, h, no, nos, nop, ntc, idr)


      write(*,*) '  Vd  = ',vd(no)


C..............................
C....   Set all splines functions on a grid r(j)  


      do 50 nc = 1, n - 2

         do 51 nn = 1, n

            c(nn) = 0.0D+00

 51      continue

         c( nc + 1) = 1.0D+00

         do 52 j = 1, no

            bf(nc,j) = bvalue( t, c, n, k, r(j), 0)

 52      continue
 50   continue

      write(*,*) '  Vd1  = ',vd(no)
C................................

      do 100 kr = 1, n - 2
         do j = 1, no

            br(j) = bf(kr,j)

         end do

C................................
C........ Exchange potential 


         call xchb(br,p,yk,lang,no,dr,nos,nop,ntc,h,r,idr)


         do 101 kc = 1, kr

            do j = 1, no

               bc(j) = bf(kc,j)

            end do
            
c-----------------  Calculate direct integral -----------------

            if(dr(1).gt.1.0D-15) ff(1) = - br(1) * bc(1) * vd(1) / dr(1)
            
            do j = 2, no

               ff(j) = -br(j) * bc(j) * vd(j) / dr(j)
            end do

            ad(kr,kc) = rint(ff, 1, no, 14, h)
            ad(kc,kr) = ad(kr, kc)

c-----------------  Calculate exchange integral -----------------
            do j = 2, no

               ff(j) = yk(j) * bc(j) / dr(j)

            end do

c            write(*,*) 'ff='
c            write(*,*) (ff(j),j=2,no,50)

            ae(kr,kc) = rint(ff, 1, no, 14, h)
            ae(kc,kr) = ae(kr,kc)

 101     continue
 100  continue


c      write(*,*) 'ad_sub'
c      write(*,'(6e10.3)') ((ad(i,j),i=1,n-2,10),j=1,n-2,10)
c      write(*,*) 'ae_sub'
c      write(*,'(6e10.3)') ((ae(i,j),i=1,n-2,10),j=1,n-2,10)

      return
      end
c---------------------------------------------------------------------
      subroutine directb(p,ft,dr,r,h,no,nos,nop,ntc,idr)
      implicit real*8(a-h,o-z)
      include "parameter.inc"
C      parameter(np=2000, nl = 15)
      dimension p(nl,np),ft(1),yk(np),f(np)
      dimension drkr(np),rkk(np),dr(np),r(np)
C..............................................

      call rkprep(r, drkr, rkk, 0, no)

      do 10 j=1,no
         f(j)=p(1,j)
 10      ft(j)=0.d0
      if(nos.eq.0) return

      call ykfctb(f,f,yk,r,dr,drkr,rkk,h,0,no,idr)

      do 11 j=1,no
 11      ft(j)=ft(j)+2.d0*2.d0*yk(j)

      if(nos.eq.1)return

      do 12 j=1,no
 12      f(j)=p(2,j)

      call ykfctb(f,f,yk,r,dr,drkr,rkk,h,0,no,idr)

      do 13 j=1,no
         f(j)=p(nos+1,j)
 13      ft(j)=ft(j)+2.d0*2.d0*yk(j)

      call ykfctb(f,f,yk,r,dr,drkr,rkk,h,0,no,idr)

      do 14 j=1,no
 14      ft(j)=ft(j)+2.d0*6.d0*yk(j)
      if(nos.eq.2)return

      do 15 j=1,no
 15      f(j)=p(3,j)

      call ykfctb(f,f,yk,r,dr,drkr,rkk,h,0,no,idr)

      do 16 j=1,no
         f(j)=p(nos+2,j)
 16      ft(j)=ft(j)+2.d0*2.d0*yk(j)

      call ykfctb(f,f,yk,r,dr,drkr,rkk,h,0,no,idr)

      do 17 j=1,no
 17      ft(j)=ft(j)+2.d0*6.d0*yk(j)

      if(nos.eq.3)return

      do 18 j=1,no
 18      f(j)=p(4,j)

      call ykfctb(f,f,yk,r,dr,drkr,rkk,h,0,no,idr)

      do 19 j=1,no
         f(j)=p(nos+3,j)
 19      ft(j)=ft(j)+2.d0*2.d0*yk(j)

      call ykfctb(f,f,yk,r,dr,drkr,rkk,h,0,no,idr)

      do 20 j=1,no
         f(j)=p(nos+nop+1,j)
 20      ft(j)=ft(j)+2.d0*6.d0*yk(j)

      call ykfctb(f,f,yk,r,dr,drkr,rkk,h,0,no,idr)

      do 21 j=1,no
 21      ft(j)=ft(j)+2.d0*10.d0*yk(j)

      if(nos.eq.4)return

      end
c------------------------------------------------------------------------
      subroutine xchb(pr,p,yk,l,no,dr,nos,nop,ntc,h,r,idr)
      implicit real*8(a-h,o-z)
      include "parameter.inc"
C      parameter(np=2000, nl = 15)
      dimension pr(1),r(np), dr(np)
      dimension p(nl,np),f(np),yk(1),y1(np)
      dimension drkr(np),rkk(np)
C.........................................

      fs=2.0d0/dfloat(2*l+1)
      fpm=3.d0*fs*dfloat(l)/dfloat(2*l-1)
      fpp=3.d0*fs*dfloat(l+1)/dfloat(2*l+3)
      fdm=2.5d0*fpm*dfloat(l-1)/dfloat(2*l-3)
      fd=5.d0*fpm*dfloat(l+1)/dfloat(2*l+3)/3.d0
      fdp=2.5d0*fpp*dfloat(l+2)/dfloat(2*l+5)
      lpm=l-1
      lpp=l+1
      ldm=l-2
      ldp=l+2

      if(nos.eq.0) return

      do 10 j=1,no
      f(j)=p(1,j)
 10   yk(j)=0.d0

      call rkprep(r,drkr,rkk,l,no)

      call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,l,no,idr)

      if(nos.eq.0) return

      do 11 j=1,no
 11   yk(j)=yk(j)+fs*y1(j)*f(j)

      if(nos.eq.1) return

      do 30 j=1,no
 30      f(j)=p(2,j)

      call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,l,no,idr)

      do 12 j=1,no
      f(j)=p(nos+1,j)
 12   yk(j)=yk(j)+fs*y1(j)*p(2,j)

      call rkprep(r,drkr,rkk,lpp,no)

      call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,lpp,no,idr)

      do 13 j=1,no
 13   yk(j)=yk(j)+fpp*y1(j)*f(j)

       if(lpm.ge.0)then

          call rkprep(r,drkr,rkk,lpm,no)

          call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,lpm,no,idr)

          do 14 j=1,no
 14       yk(j)=yk(j)+fpm*y1(j)*f(j)
       endif

      if(nos.eq.2) return

      do 31 j=1,no
 31      f(j)=p(3,j)

      call rkprep(r,drkr,rkk,l,no)

      call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,l,no,idr)

      do 15 j=1,no
      f(j)=p(nos+2,j)
 15   yk(j)=yk(j)+fs*y1(j)*p(3,j)

      call rkprep(r,drkr,rkk,lpp,no)

      call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,lpp,no,idr)

      do 16 j=1,no
 16   yk(j)=yk(j)+fpp*y1(j)*f(j)

       if(lpm.ge.0)then

          call rkprep(r,drkr,rkk,lpm,no)
          call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,lpm,no,idr)

          do 17 j=1,no
 17       yk(j)=yk(j)+fpm*y1(j)*f(j)
       endif

      if(nos.eq.3) return

      do 32 j=1,no
 32      f(j)=p(4,j)

      call rkprep(r,drkr,rkk,l,no)
      call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,l,no,idr)

      do 20 j=1,no
      f(j)=p(nos+3,j)
 20   yk(j)=yk(j)+fs*y1(j)*p(4,j)

      call rkprep(r,drkr,rkk,lpp,no)
      call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,lpp,no,idr)

      do 21 j=1,no
 21   yk(j)=yk(j)+fpp*y1(j)*f(j)

       if(lpm.ge.0)then

          call rkprep(r,drkr,rkk,lpm,no)
          call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,lpm,no,idr)

          do 22 j=1,no
 22       yk(j)=yk(j)+fpm*y1(j)*f(j)
       endif

      do 33 j=1,no
 33      f(j)=p(nos+nop+1,j)

      call rkprep(r,drkr,rkk,ldp,no)
      call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,ldp,no,idr)

      do 23 j=1,no
 23   yk(j)=yk(j)+fdp*y1(j)*f(j)

      call rkprep(r,drkr,rkk,l,no)
      call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,l,no,idr)

      do 24 j=1,no
 24   yk(j)=yk(j)+fd*y1(j)*f(j)

       if(ldm.ge.0)then
          call rkprep(r,drkr,rkk,ldm,no)
          call ykfctb(pr,f,y1,r,dr,drkr,rkk,h,ldm,no,idr)

          do 25 j=1,no
 25       yk(j)=yk(j)+fdm*y1(j)*f(j)
       endif

      return
      end
c -----------------------------------------------------------------
      subroutine ykfctb(p,pp,yk,r,dr,drkr,rkk,h,k,no,idr)
      implicit real*8(a-h,o-z)
      include "parameter.inc"
C      parameter(np=2000)
      dimension p(np),pp(np),r(np),dr(np),drkr(np),rkk(np)
      dimension yk(np)

      if(idr.eq.0) call yksinb(p,pp,yk,r,dr,drkr,rkk,h,k,no)

      if(idr.eq.1) call ykexpb(p,pp,yk,r,dr,drkr,rkk,h,k,no)

      if(idr.gt.1) stop

      return
      end
c -----------------------------------------------------------------
c --- yk function in  r = ri * exp( (x+xi)**(1/nrp) )  scale
c --- x = (n-1) * h
c --- drkr = 1.0/r**kr
c --- rkk  = r**kk
c -
      subroutine ykexpb(p,pp,yk,r,dr,drkr,rkk,h,k,no)
      implicit real*8(a-h,o-z)
      include "parameter.inc"
C      parameter(nomx=2000)
      dimension p(np),pp(np),r(np),dr(np),drkr(np),rkk(np)
      dimension yk(np),zk(np),wy(np),wz(np)
      n=1
      kr=k+n-1
      kk=k+1

      if(kr .eq. 0) go to 30

      do 10 j=2,no
      wz(j)=p(j)*pp(j)*drkr(j)/dr(j)
   10 wy(j)=p(j)*pp(j)*rkk(j)/dr(j)

      call yint(wy,wz,yk,zk,no,h)

      do 20 j=2,no
   20 yk(j)=yk(j)*drkr(j)+zk(j)*rkk(j)

      return

   30 do 40 j=2,no
      wz(j)=p(j)*pp(j)/dr(j)
   40 wy(j)=p(j)*pp(j)*rkk(j)/dr(j)

      call yint(wy,wz,yk,zk,no,h)

      do 50 j=2,no
   50 yk(j)=yk(j)+zk(j)*rkk(j)

      return
      end

c ------------------------------------------------------------------
 
c --- drkr = 1.0/r**kr
c --- rkk  = r**kk
c -
      subroutine yksinb(p,pp,yk,r,dr,drkr,rkk,h,k,no)
      implicit real*8(a-h,o-z)
      include "parameter.inc"
C      parameter(np=2000)
      dimension p(np),pp(np),r(np),dr(np),drkr(np),rkk(np)
      dimension yk(np),zk(np),wy(np),wz(np)

      n=1
      kr=k+n-1
      kk=k+1
      wz(1)=0.0d0
      wy(1)=0.0d0

      if(kr .eq. 0) go to 30
      do 10 j=2,no
      wz(j)=p(j)*pp(j)*drkr(j)/dr(j)
   10 wy(j)=p(j)*pp(j)*rkk(j)/dr(j)

      call yint(wy,wz,yk,zk,no,h)
      yk(1)=0.0d0
      do 20 j=2,no
   20 yk(j)=yk(j)*drkr(j)+zk(j)*rkk(j)
      return

   30 do 40 j=2,no
      wz(j)=p(j)*pp(j)/dr(j)
   40 wy(j)=p(j)*pp(j)*rkk(j)/dr(j)
      call yint(wy,wz,yk,zk,no,h)
      yk(1)=0.0d0
      do 50 j=2,no
   50 yk(j)=yk(j)+zk(j)*rkk(j)
      return
      end
c ------------------------------------------------------------------
      subroutine rkprep(r,drkr,rkk,k,no)
      implicit real*8(a-h,o-z)
** this routine prepares r**k and 1/r**(k+1) for use in routine ykfctb.
c -----------------------------
      include "parameter.inc"
C      parameter(nomx=2000)
c -----------------------------
      dimension drkr(np),rkk(np),r(np)
      dimension rpw(10,np)
      do 28 kk=1,10
         do 28 j=1,no
            rpw(kk,j)=r(j)**kk
 28   continue

      n=1
      kr=k+n-1
      kk=k+1
      if(kr.eq.0) go to 30
      drkr(1)=0.0d0
      if(r(1) .ne. 0.0d0) drkr(1)=1.0d0/rpw(kr,1)
      do 20 j=2,no
   20 drkr(j)=1.0d0/rpw(kr,j)
   30 continue
      do 21 j=1,no
   21 rkk(j)=rpw(kk,j)
      return
      end
c--------------------------------------------------------
      subroutine corewf(c,p,r,dr,t,n,k,ni,nf,nx,no,ntc,h)
      implicit real*8(a-h,o-z)
      include "parameter.inc"
C      parameter(ns=62, np = 2000,nl= 15, kx = 15)
      dimension c(nx,nx),p(nl,np),f(np)
      dimension r(np),dr(np),coe(nx),t(np+kx)
C.....................................
      in = n - 1
      coe(1) = 0.d0
      coe(n) = 0.d0

      do 20 npp = ni, nf
        in = in - 1

        do 10 i = 1, n - 2
 10        coe(i+1)=c(i,in)

        p(npp,1) = 0.0D+00
        do 12 j = 2, no
 12        p(npp, j) = bvalue(t,coe,n,k,r(j),0)

c        if(p(npp,6).gt.0.d0) goto 14
c        do 13 j=1,no
c 13        p(np,j)=-p(np,j)
c 14     continue

        f(1) = 0.0D+00
        do 15 j = 2, no
 15        f(j) = p(npp,j) * p(npp,j) * r(j) / dr(j)

        oth = rint(f, 1, no, 14, h)

        write(*,*) ' Core WF #',npp,' norm =',oth

        sq2 = dsqrt( oth )

        do 16 j = 1, no
 16        p(npp,j) = p(npp,j) / sq2

 20   continue

      return
      end
c$$$c********************************************************
c$$$      function bvalue(t,bcoef,n,k,x,jderiv)
c$$$      implicit real*8(a-h,o-z)
c$$$      parameter(KMAX=15)
c$$$      dimension bcoef(1),t(1),aj(KMAX),dl(KMAX),dr(KMAX)
c$$$c      dimension bcoef(n),t(1),aj(KMAX),dl(KMAX),dr(KMAX)
c$$$
c$$$      bvalue = 0.0D0
c$$$      if(jderiv .ge. k)                go to 99
c$$$      call interv (t,n+k,x,i,mflag)
c$$$      if(mflag .ne. 0)                 go to 99
c$$$
c$$$      km1 = k - 1
c$$$      if(km1 . gt. 0)                  go to 1
c$$$      bvalue = bcoef(i)
c$$$                                       go to 99
c$$$
c$$$ 1    jcmin = 1
c$$$      imk = i - k
c$$$      if(imk .ge. 0)                    go to 8
c$$$      jcmin = 1 - imk
c$$$
c$$$      do 5 j = 1,i
c$$$         dl(j) = x - t(i+1-j)
c$$$ 5    continue
c$$$      do 6 j = i,km1
c$$$         aj(k-j)  = 0.0D0
c$$$         dl(j) = dl(i)
c$$$ 6    continue
c$$$                                      go to 10
c$$$
c$$$ 8    do 9 j = 1,km1
c$$$         dl(j) = x - t(i+1-j)
c$$$ 9    continue
c$$$
c$$$ 10   jcmax = k
c$$$      nmi = n - i
c$$$      if(nmi .ge. 0)                  go to 18
c$$$      jcmax = k + nmi
c$$$      do 15 j = 1,jcmax
c$$$          dr(j) = t(i+j) - x
c$$$ 15   continue
c$$$      do 16 j = jcmax,km1
c$$$          aj(j+1) = 0.0D0
c$$$          dr(j) = dr(jcmax)
c$$$ 16   continue
c$$$                                     go to 20
c$$$
c$$$ 18   do 19 j = 1,km1
c$$$        dr(j) = t(i+j) - x
c$$$ 19   continue
c$$$
c$$$ 20   do 21 jc = jcmin,jcmax
c$$$          aj(jc) = bcoef(imk + jc)
c$$$ 21   continue
c$$$
c$$$      if(jderiv .eq. 0)             go to 30
c$$$      do 23 j = 1,jderiv
c$$$         kmj = k - j
c$$$         fkmj = DFLOAT(kmj)
c$$$         ilo = kmj
c$$$         do 22 jj = 1,kmj
c$$$            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
c$$$            ilo = ilo - 1
c$$$ 22      continue
c$$$ 23   continue
c$$$
c$$$ 30   if(jderiv .eq. km1)           go to 39
c$$$      do 33 j = jderiv+1,km1
c$$$         kmj = k - j
c$$$         ilo = kmj
c$$$         do 32 jj = 1,kmj
c$$$           aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
c$$$           ilo = ilo - 1
c$$$ 32      continue
c$$$ 33   continue
c$$$ 39   bvalue = aj(1)
c$$$ 
c$$$ 99   RETURN
c$$$      end
c$$$
c$$$
c$$$      subroutine interv (xt,lxt,x,left,mflag)
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension xt(lxt)
c$$$      data ilo /1/
c$$$
c$$$      ihi = ilo + 1
c$$$      if(ihi .lt. lxt)                go to 20
c$$$      if(x .ge. xt(lxt))              go to 110
c$$$      if(lxt . le. 1)                 go to 90
c$$$      ilo = lxt - 1
c$$$      ihi = lxt
c$$$
c$$$ 20   if(x .ge. xt(ihi))              go to 40
c$$$      if(x .ge. xt(ilo))              go to 100
c$$$
c$$$      istep = 1
c$$$ 31   continue
c$$$         ihi = ilo
c$$$         ilo = ihi - istep
c$$$         if (ilo .le. 1)              go to 35
c$$$         if (x .ge. xt(ilo))          go to 50
c$$$         istep = istep*2
c$$$      go to 31
c$$$ 35   ilo = 1
c$$$      if(x .lt. xt(1))                go to 90
c$$$                                      go to 50
c$$$
c$$$ 40   istep = 1
c$$$ 41   continue
c$$$         ilo = ihi
c$$$         ihi = ilo + istep
c$$$         if(ihi .ge. lxt)             go to 45
c$$$         if(x .lt. xt(ihi))           go to 50
c$$$         istep = istep*2
c$$$      go to 41
c$$$
c$$$ 45   if (x .ge. xt(lxt))             go to 110
c$$$      ihi = lxt
c$$$ 50   continue
c$$$            middle = (ilo + ihi)/2
c$$$            if(middle .eq. ilo)       go to 100
c$$$            if(x .lt. xt(middle))     go to 53
c$$$            ilo = middle
c$$$          go to 50
c$$$ 53       ihi = middle
c$$$      go to 50
c$$$
c$$$ 90   mflag = -1
c$$$      left = 1
c$$$      RETURN
c$$$
c$$$ 100  mflag = 0
c$$$      left = ilo
c$$$      RETURN
c$$$
c$$$ 110  mflag = 1
c$$$      left = lxt
c$$$      RETURN
c$$$      end
c$$$
c$$$
c$$$      subroutine bsplvd(t,k,x,left,dbiatx,nderiv)
c$$$      implicit real*8(a-h,o-z)
c$$$*******  changes from de Boor
c$$$      parameter(KX=15)
c$$$      dimension a(KX,KX),dbiatx(KX,nderiv)
c$$$*  also a missing in arg list
c$$$*****************
c$$$      dimension t(1)
c$$$
c$$$      mhigh = max0(min0(nderiv,k),1)
c$$$      kp1 = k + 1
c$$$      call bsplvb(t,kp1-mhigh,1,x,left,dbiatx)
c$$$      if(mhigh .eq. 1) go to 99
c$$$
c$$$      ideriv = mhigh
c$$$      do 15 m = 2,mhigh
c$$$          jp1mid = 1
c$$$          do 11 j = ideriv,k
c$$$             dbiatx(j,ideriv) = dbiatx(jp1mid,1)
c$$$             jp1mid = jp1mid + 1
c$$$ 11       continue
c$$$          ideriv = ideriv - 1
c$$$          call bsplvb(t,kp1-ideriv,2,x,left,dbiatx)
c$$$ 15   continue
c$$$
c$$$      jlow = 1
c$$$      do 20 i = 1,k
c$$$         do 19 j = jlow,k
c$$$            a(j,i) = 0.0D0
c$$$ 19      continue
c$$$         jlow = i
c$$$         a(i,i) = 1.0D0
c$$$ 20   continue
c$$$
c$$$      do 40 m = 2,mhigh
c$$$         kp1mm = kp1 - m
c$$$         fkp1mm = DFLOAT(kp1mm)
c$$$         il = left
c$$$         i = k
c$$$
c$$$         do 25 ldummy = 1,kp1mm
c$$$            factor = fkp1mm/(t(il+kp1mm) - t(il))
c$$$            do 24 j = 1,i
c$$$               a(i,j) = (a(i,j) -a(i-1,j))*factor
c$$$ 24         continue
c$$$            il = il - 1
c$$$            i = i - 1
c$$$ 25      continue
c$$$
c$$$         do 36 i = 1,k
c$$$              sum = 0.0D0
c$$$              jlow = max0(i,m)
c$$$              do 35 j = jlow,k
c$$$                 sum = a(j,i)*dbiatx(j,m) + sum
c$$$ 35           continue
c$$$              dbiatx(i,m) = sum
c$$$ 36      continue
c$$$ 40   continue
c$$$
c$$$ 99   RETURN
c$$$      end
c$$$********************************************************
c$$$      subroutine bsplvb(t,jhigh,iindex,x,left,biatx)
c$$$      implicit real*8(a-h,o-z)
c$$$      parameter(JMAX=20)
c$$$      dimension biatx(jhigh),t(1),deltal(jmax),deltar(jmax)
c$$$      data j/1/
c$$$
c$$$      go to (10,20),iindex
c$$$
c$$$ 10   j = 1
c$$$      biatx(1) = 1.0D0
c$$$      if(j .ge. jhigh)                 go to 99
c$$$
c$$$ 20   continue
c$$$         jp1 = j + 1
c$$$         deltar(j) = t(left+j) - x
c$$$         deltal(j) = x - t(left+1-j)
c$$$         saved = 0.0D0
c$$$         do 26 i = 1,j
c$$$             term = biatx(i)/(deltar(i) + deltal(jp1-i))
c$$$             biatx(i) = saved + deltar(i)*term
c$$$             saved = deltal(jp1-i)*term
c$$$ 26      continue
c$$$         biatx(jp1) = saved
c$$$         j = jp1
c$$$         if(j .lt. jhigh)
c$$$     *go to 20
c$$$
c$$$ 99   RETURN
c$$$      end
c$$$
c$$$
c$$$c ------------------------------------------------------------------
c$$$      subroutine yint (v,w,y,z,m,h)
c$$$c
c$$$c  this program calculates the indefinite integrals y and z using the
c$$$c  lagrange integration formula
c$$$c
c$$$c  y(r) = integral of v from 0 to r
c$$$c  z(r) = integral of w from r to infinity
c$$$c  m is the maximum tabulation point of v and w (virtual infinity)
c$$$c  h is the step size of the radial grid
c$$$c
c$$$      IMPLICIT REAL*8(A-H,O-Z)
c$$$      dimension v(1),w(1),y(1),z(1)
c$$$c ------------------------------------------------------------------
c$$$c                lagrange 10 point integration formula
c$$$c               --------------------------------------
c$$$      dimension aa(5,10),a(10,5),b(5)
c$$$c     equivalence (b(1),a(6,5))
c$$$      data ia/5/, ja/10/, da/7257600.0/
c$$$      data aa/2082753.0, -57281.0,   10625.0,   -3969.0,   2497.0,
c$$$     &       9449717.0, 2655563.0, -163531.0,   50315.0, -28939.0,
c$$$     &     -11271304.0, 6872072.0, 3133688.0, -342136.0, 162680.0,
c$$$     &      16002320.0,-4397584.0, 5597072.0, 3609968.0,-641776.0,
c$$$     &     -17283646.0, 3973310.0,-2166334.0, 4763582.0,4134338.0,
c$$$     &      13510082.0,-2848834.0, 1295810.0,-1166146.0,4134338.0,
c$$$     &      -7394032.0, 1481072.0, -617584.0,  462320.0,-641776.0,
c$$$     &       2687864.0, -520312.0,  206072.0, -141304.0, 162680.0,
c$$$     &       -583435.0,  110219.0,  -42187.0,   27467.0, -28939.0,
c$$$     &         57281.0,  -10625.0,    3969.0,   -2497.0,   2497.0/
c$$$      data b/ 4134338.0,-641776.0,  162680.0,  -28939.0,   2497.0/
c$$$      data h0/0.0/
c$$$c
c$$$c  note that a different even order method can be used by replacing the
c$$$c  dimension and data statements in this block (also in do loop 20)
c$$$c ------------------------------------------------------------------
c$$$      if(h.eq.h0) go to 5
c$$$      hd=h/da
c$$$      do 2 i=1,ia
c$$$         do 1 j=1,ja
c$$$            a(j,i)=aa(i,j)*hd
c$$$ 1       continue
c$$$         b(i) = b(i)*hd
c$$$ 2    continue
c$$$      h0=h
c$$$ 5    y(1)=0.0D0
c$$$      z(m)=0.0D0
c$$$      do 10 i=2,ia
c$$$      k=m-i+1
c$$$      y(i)=y(i-1)
c$$$      z(k)=z(k+1)
c$$$      ii=i-1
c$$$      do 10 j=1,ja
c$$$      y(i)=y(i)+a(j,ii)*v(j)
c$$$ 10   z(k)=z(k)+a(j,ii)*w(m-j+1)
c$$$      im=ia+1
c$$$      in=m-ia+1
c$$$      do 20 i=im,in
c$$$      k=m-i+1
c$$$      y(i)=y(i-1)+b(1)*(v(i  )+v(i-1))+b(2)*(v(i+1)+v(i-2))
c$$$     &           +b(3)*(v(i+2)+v(i-3))+b(4)*(v(i+3)+v(i-4))
c$$$     &           +b(5)*(v(i+4)+v(i-5))
c$$$ 20   z(k)=z(k+1)+b(1)*(w(k  )+w(k+1))+b(2)*(w(k-1)+w(k+2))
c$$$     &           +b(3)*(w(k-2)+w(k+3))+b(4)*(w(k-3)+w(k+4))
c$$$     &           +b(5)*(w(k-4)+w(k+5))
c$$$      in=in+1
c$$$      do 30 i=in,m
c$$$      k=m-i+1
c$$$      y(i)=y(i-1)
c$$$      z(k)=z(k+1)
c$$$      do 30 j=1,ja
c$$$      y(i)=y(i)+a(j,k)*v(m-j+1)
c$$$ 30   z(k)=z(k)+a(j,k)*w(j)
c$$$      return
c$$$      end
c$$$C ------------------------------------------------------------------
c$$$      subroutine gauss(k,x,w)
c$$$************************************************
c$$$*
c$$$*   x(i) = gaussian coordinates for [0,1]
c$$$*   w(i) = gaussian coordinates for [0,1]
c$$$*   1<= k <= 15   for k point case
c$$$*
c$$$*************************************************
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension x(k),w(k)
c$$$
c$$$      if(k.lt.1.or.k.gt.15) go to 901
c$$$      goto(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150),k
c$$$ 10   x(1) = 0.5D 00
c$$$      w(1) = 1.0D0
c$$$      go to 900
c$$$ 20   x(1) = .211324865405187
c$$$      x(2) = .788675134594813
c$$$      w(1) = 0.5D 00 
c$$$      w(2) = 0.5D 00
c$$$      go to 900
c$$$ 30   x(1) = .112701665379258 
c$$$      x(2) = 0.5D 00               
c$$$      x(3) = .887298334620742
c$$$      w(1) = .277777777777778 
c$$$      w(2) = .444444444444444 
c$$$      w(3) = .277777777777778
c$$$      go to 900
c$$$ 40   x(1) = .0694318442029737 
c$$$      x(2) = .330009478207572  
c$$$      x(3) = .669990521792428  
c$$$      x(4) = .930568155797026
c$$$      w(1) = .173927422568727 
c$$$      w(2) = .326072577431273 
c$$$      w(3) = .326072577431273 
c$$$      w(4) = .173927422568727
c$$$      go to 900
c$$$ 50   x(1) = .046910077030668 
c$$$      x(2) = .230765344947158 
c$$$      x(3) = 0.5D 00               
c$$$      x(4) = .769234655052842
c$$$      x(5) = .953089922969332
c$$$      w(1) = .118463442528095 
c$$$      w(2) = .239314335249683 
c$$$      w(3) = .284444444444444 
c$$$      w(4) = .239314335249683
c$$$      w(5) = .118463442528095
c$$$      go to 900
c$$$ 60   x(1) = .033765242898424 
c$$$      x(2) = .169395306766868 
c$$$      x(3) = .380690406958402 
c$$$      x(4) = .619309593041598
c$$$      x(5) = .830604693233132 
c$$$      x(6) = .966234757101576
c$$$      w(1) = .0856622461895852 
c$$$      w(2) = .180380786524069  
c$$$      w(3) = .233956967286346  
c$$$      w(4) = .233956967286346
c$$$      w(5) = .180380786524069
c$$$      w(6) = .0856622461895852
c$$$      go to 900
c$$$ 70   x(1) = .0254460438286207 
c$$$      x(2) = .129234407200303  
c$$$      x(3) = .297077424311301  
c$$$      x(4) = 0.5D 00
c$$$      x(5) = .702922575688699
c$$$      x(6) = .870765592799697
c$$$      x(7) = .974553956171379
c$$$      w(1) = .0647424830844348 
c$$$      w(2) = .139852695744638  
c$$$      w(3) = .19091502525256   
c$$$      w(4) = .208979591836735
c$$$      w(5) = .19091502525256
c$$$      w(6) = .139852695744638
c$$$      w(7) = .0647424830844348
c$$$      go to 900
c$$$ 80   x(1) = .0198550717512319 
c$$$      x(2) = .101666761293187  
c$$$      x(3) = .237233795041835  
c$$$      x(4) = .408282678752175
c$$$      x(5) = .591717321247825  
c$$$      x(6) = .762766204958164  
c$$$      x(7) = .898333238706813  
c$$$      x(8) = .980144928248768
c$$$      w(1) = .0506142681451881 
c$$$      w(2) = .111190517226687  
c$$$      w(3) = .156853322938944  
c$$$      w(4) = .181341891689181
c$$$      w(5) = .181341891689181  
c$$$      w(6) = .156853322938944  
c$$$      w(7) = .111190517226687  
c$$$      w(8) = .0506142681451881
c$$$      go to 900
c$$$ 90   x(1) = .015919880246187  
c$$$      x(2) = .0819844463366821 
c$$$      x(3) = .193314283649705  
c$$$      x(4) = .337873288298095
c$$$      x(5) = 0.5D 00                
c$$$      x(6) = .662126711701904  
c$$$      x(7) = .806685716350295  
c$$$      x(8) = .918015553663318
c$$$      x(9) = .984080119753813
c$$$      w(1) = .0406371941807872 
c$$$      w(2) = .0903240803474287 
c$$$      w(3) = .130305348201468  
c$$$      w(4) = .156173538520001
c$$$      w(5) = .16511967750063   
c$$$      w(6) = .156173538520001  
c$$$      w(7) = .130305348201468  
c$$$      w(8) = .0903240803474287
c$$$      w(9) = .0406371941807872
c$$$      go to 900
c$$$ 100  x(1) = .0130467357414141 
c$$$      x(2) = .0674683166555077 
c$$$      x(3) = .160295215850488  
c$$$      x(4) = .283302302935376
c$$$      x(5) = .425562830509184  
c$$$      x(6) = .574437169490816  
c$$$      x(7) = .716697697064624  
c$$$      x(8) = .839704784149512
c$$$      x(9) = .932531683344492
c$$$      x(10)= .986953264258586
c$$$      w(1) = .0333356721543441 
c$$$      w(2) = .0747256745752903 
c$$$      w(3) = .109543181257991  
c$$$      w(4) = .134633359654998
c$$$      w(5) = .147762112357376  
c$$$      w(6) = .147762112357376  
c$$$      w(7) = .134633359654998
c$$$      w(8) = .109543181257991
c$$$      w(9) = .0747256745752903
c$$$      w(10)= .0333356721543441
c$$$      go to 900
c$$$ 110  x(1) = .0108856709269715 
c$$$      x(2) = .0564687001159523 
c$$$      x(3) = .134923997212975  
c$$$      x(4) = .240451935396594
c$$$      x(5) = .365228422023827  
c$$$      x(6) = 0.5D 00                
c$$$      x(7) = .634771577976172  
c$$$      x(8) = .759548064603406
c$$$      x(9) = .865076002787025  
c$$$      x(10)= .943531299884048  
c$$$      x(11)= .989114329073028
c$$$      w(1) = .0278342835580868 
c$$$      w(2) = .0627901847324523 
c$$$      w(3) = .0931451054638672 
c$$$      w(4) = .116596882295995
c$$$      w(5) = .131402272255123  
c$$$      w(6) = .13646254338895   
c$$$      w(7) = .131402272255123  
c$$$      w(8) = .116596882295995
c$$$      w(9) = .0931451054638672 
c$$$      w(10)= .0627901847324523 
c$$$      w(11)= .0278342835580868
c$$$      go to 900
c$$$ 120  x(1) = .00921968287664038 
c$$$      x(2) = .0479413718147626  
c$$$      x(3) = .115048662902848
c$$$      x(4) = .206341022856691   
c$$$      x(5) = .31608425050091    
c$$$      x(6) = .437383295744266
c$$$      x(7) = .562616704255734   
c$$$      x(8) = .68391574949909 
c$$$      x(9) = .793658977143309
c$$$      x(10)= .884951337097152
c$$$      x(11)= .952058628185237   
c$$$      x(12)= .99078031712336
c$$$      w(1) = .0235876681932559 
c$$$      w(2) = .0534696629976592 
c$$$      w(3) = .0800391642716731 
c$$$      w(4) = .101583713361533
c$$$      w(5) = .116746268269177  
c$$$      w(6) = .124573522906701  
c$$$      w(7) = .124573522906701  
c$$$      w(8) = .116746268269177
c$$$      w(9) = .101583713361533  
c$$$      w(10)= .0800391642716731 
c$$$      w(11)= .0534696629976592 
c$$$      w(12)= .0235876681932559
c$$$      go to 900
c$$$ 130  x(1) = .00790847264070593 
c$$$      x(2) = .041200800388511   
c$$$      x(3) = .099210954633345
c$$$      x(4) = .17882533027983    
c$$$      x(5) = .275753624481777   
c$$$      x(6) = .384770842022433
c$$$      x(7) = 0.5D 00                 
c$$$      x(8) = .615229157977567   
c$$$      x(9) = .724246375518223
c$$$      x(10)= .82117466972017    
c$$$      x(11)= .900789045366655   
c$$$      x(12)= .958799199611489
c$$$      x(13)= .992091527359294
c$$$      w(1) = .0202420023826579 
c$$$      w(2) = .0460607499188642 
c$$$      w(3) = .0694367551098937 
c$$$      w(4) = .0890729903809729
c$$$      w(5) = .103908023768444  
c$$$      w(6) = .113141590131449  
c$$$      w(7) = .116275776615437  
c$$$      w(8) = .113141590131449
c$$$      w(9) = .103908023768444  
c$$$      w(10)= .0890729903809729 
c$$$      w(11)= .0694367551098937 
c$$$      w(12)= .0460607499188642
c$$$      w(13)= .0202420023826579
c$$$      go to 900
c$$$ 140  x(1) = .00685809565159383 
c$$$      x(2) = .0357825581682132  
c$$$      x(3) = .0863993424651175
c$$$      x(4) = .156353547594157   
c$$$      x(5) = .242375681820923   
c$$$      x(6) = .340443815536055
c$$$      x(7) = .445972525646328   
c$$$      x(8) = .554027474353672   
c$$$      x(9) = .659556184463945
c$$$      x(10)= .757624318179077   
c$$$      x(11)= .843646452405843   
c$$$      x(12)= .913600657534882
c$$$      x(13)= .964217441831787   
c$$$      x(14)= .993141904348406
c$$$      w(1) = .0175597301658759 
c$$$      w(2) = .0400790435798801 
c$$$      w(3) = .0607592853439516 
c$$$      w(4) = .0786015835790968
c$$$      w(5) = .092769198738969  
c$$$      w(6) = .102599231860648  
c$$$      w(7) = .107631926731579  
c$$$      w(8) = .107631926731579
c$$$      w(9) = .102599231860648  
c$$$      w(10)= .092769198738969  
c$$$      w(11)= .0786015835790968 
c$$$      w(12)= .0607592853439516
c$$$      w(13)= .0400790435798801 
c$$$      w(14)= .0175597301658759
c$$$      go to 900
c$$$ 150  x(1) = .00600374098975728 
c$$$      x(2) = .031363303799647   
c$$$      x(3) = .0758967082947864
c$$$      x(4) = .137791134319915   
c$$$      x(5) = .214513913695731   
c$$$      x(6) = .302924326461218
c$$$      x(7) = .399402953001283   
c$$$      x(8) = 0.5D 00                 
c$$$      x(9) = .600597046998717
c$$$      x(10)= .697075673538782   
c$$$      x(11)= .785486086304269   
c$$$      x(12)= .862208865680085
c$$$      x(13)= .924103291705214   
c$$$      x(14)= .968636696200353
c$$$      x(15)= .993996259010243
c$$$      w(1) = .0153766209980586 
c$$$      w(2) = .0351830237440541 
c$$$      w(3) = .053579610233586  
c$$$      w(4) = .0697853389630772
c$$$      w(5) = .083134602908497  
c$$$      w(6) = .0930805000077812 
c$$$      w(7) = .0992157426635559 
c$$$      w(8) = .101289120962781
c$$$      w(9) = .0992157426635559 
c$$$      w(10)= .0930805000077812 
c$$$      w(11)= .083134602908497  
c$$$      w(12)= .0697853389630772
c$$$      w(13)= .053579610233586  
c$$$      w(14)= .0351830237440541 
c$$$      w(15)= .0153766209980586
c$$$
c$$$ 900  return
c$$$
c$$$ 901  stop 'error in gauss'
c$$$      end
