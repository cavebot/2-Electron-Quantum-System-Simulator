c------------------------------------------------------------------ 
      subroutine vhfbsp(t,r,dr,p,ad,ae,h,n,k,no,lang,nos,nop,ntc,idr)
      implicit real*8(a-h,o-z)
C      parameter(ns=62, np=2000, nl = 15, kx = 15)
      include "parameter.1e.inc"
      dimension t(ns+nl),r(np),dr(np),p(nl,np),c(ns),vd(np), yk(np)
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
      include "parameter.1e.inc"
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
      include "parameter.1e.inc"
C      parameter(np=2000, nl = 15)
      dimension pr(np),r(np), dr(np)
      dimension p(nl,np),f(np),yk(np),y1(np)
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

C      call ykfctb(f,f,yk,r,dr,drkr,rkk,h,0,no,idr)

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
      subroutine ykfctb(p, pp, yk,r,dr,drkr,rkk,h,k,no,idr)
      implicit real*8(a-h,o-z)
      include "parameter.1e.inc"
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
      include "parameter.1e.inc"
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
      include "parameter.1e.inc"
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
      include "parameter.1e.inc"
c -----------------------------
      dimension drkr(np),rkk(np),r(np)
      dimension rpw(10,np)
C
      do  kk = 1, 10
        do  j = 1,no
          rpw(kk,j) = r(j)**kk          
        enddo
      enddo

      n  = 1
      kr = k + n - 1
      kk = k + 1

      if(kr.eq.0) go to 30
      
      drkr(1)=0.0d0
      if(r(1) .ne. 0.0d0) drkr(1)=1.0d0/rpw(kr,1)

      do  j = 2, no
        drkr(j) = 1.0d0 / rpw(kr,j)
      enddo

 30   continue
      do  j = 1, no
        rkk(j)=rpw(kk,j)
      enddo

      return
      end
c--------------------------------------------------------
      subroutine corewf(c,p,r,dr,t,n,k,ni,nf,nx,no,ntc,h)
      implicit real*8(a-h,o-z)
C      parameter(ns=62, np = 2000,nl= 15, kx = 15)
      include "parameter.1e.inc"
      dimension c(nx,nx),p(nl,np),f(np)
      dimension r(np),dr(np),coe(nx),t(np+nl)
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
