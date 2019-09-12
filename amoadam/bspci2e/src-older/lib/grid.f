c**  mkgrid   === 1 ===
c**  step     === 2 ===

c*********************************************************
c*
c*
c*  input :
c*         n = number of collacation points
c*         k = order of spline > 2
c*        a,b = 1st and last points on interval
c*  output :
c*         t(1) = t(2) = ... t(k) = a
c*         t(n+1) = t(n+2) = ... t(n+k)
c*         t(i) = a + (i-k)*h   i = k+1,n
c*           with  h = (b-a)/(n-k+1)
c*
c***********************************************************
C#  sub1
      subroutine mkgrid(n,k,rmax,dr0,t)
      implicit doubleprecision(a-h,o-z)
      dimension t(1)

      do 50 i = 1,k
         t(i) = 0.0D+00
 50   continue

      Hh=RMax/DFLOAT(N-k+1)
      PI=3.14159265358979324D 00
      Y=-DLOG(2.0D0*DASIN(dr0/RMax)/PI)/DLOG(DFLOAT(N-k+1))
c      write(6,1000) ' t0 =',dr0,'  h =',hh
      do 101 i = k+1,n

         FH=DFLOAT(i-k)*Hh
         AUG=(PI/2.0D0)*(FH/Rmax)**Y
         t(i)=RMax*DSIN(AUG)
 101  continue
      

 1000 format(/a,1pd14.6,a,1pd14.6/)

      do 150 i = n+1,n+k
         t(i) = rmax
 150  continue

      return
      end
c*****************************************************************
C#  sub2
      subroutine step(n,k,rmax,dr,t0,h)

      implicit doubleprecision(a-h,o-z)
      data small /1.0d-10/,nmax/54/

         f(x) = 1.0D0 + con*x - (1.0D0+x)**ipow

       ipow = n-k+1
       con = rmax/dr

       fa = f(small)
       a = small
       if(fa.le.0.0D0) then
          write(6,*) 'con=',con,'  ipow=',ipow
          write(6,*) ' a =',a,  '    fa=',fa
          write(6,*) ' problem with fa !!!'
          stop
       endif
       fb = f(1.0D0)
       b = 1.0D0

       if(fb.ge.0.0D0) then
            write(6,*) ' b =',b,'   fb=',fb
            write(6,*) ' problem with fb !!!'
            stop
       endif
       dx = 1.0d0
       i = 0
 200   i = i+1
       dx = 0.5D0*dx
       c = 0.5D0*(a+b)
       fc = f(c)
       if(fb*fc.gt.0.0D0) then
          b = c
          fb = fc
       else
          a = c
          fa = fc
       endif
c      write(6,*) ' i=',i,' c=',c,'  f=',fc
       if(fc.ne.0.0D0.and.i.lt.nmax) go to 200

       h = dlog(1.0D0+c)
       t0 = dr/c

       return
       end
********************************************************
      subroutine rin(r,dr,h,no,id)
      implicit real*8(a-h,o-z)
      dimension r(10),dr(1)
    1 format(2x,'r(2) & rmx=?'/)
    2 format(2x,'no=?'/)
    3 format(2x,'sine scale:  r(2),h,rmx & no=',1p3e14.6,i5/)
    7 format(2x,'xi,r(1) & rmx=?'/)
    8 format(2x,'no & nrp=?'/)
    9 format(2x,'r in sin scale: type 0 & r in exp scale: type 1'/)
   12 format(2x,'xi,r(1),rmx & h=',1p4e13.5,', no & nrp=',2i5/)
c------------

      open(10,file='inp/grid.inp')

      read(10,*) id
      if(id .eq. 1) go to 30
      if(id .ne. 0) stop

      read(10,*) xi, r(2), rmx
      read(10,*) no

      call rsin(r,dr,rmx,h,no)

c      write(*,3)  r(2), h, rmx, no
      write(16,3) r(2), h, rmx, no

      close(10)
      return
   30 continue

      read(10,*) xi, r(1), rmx
      read(10,*) no, nrp

      call rexp(r,dr,xi,rmx,h,no,nrp)

      write(16,12) xi,r(1),rmx,h,no,nrp

      close(10)

      return
      end
c ------------------------------------------------------------------
      subroutine rexp(r,dr,xi,rmx,h,no,nrp)
  
c --- r = ri*dexp( (x+xi)**(1/nrp))
c --- dr = dx  r/[nrp * (ln(r/ri))**(nrp-1)]
   
      implicit real*8(a-h,o-z)
      dimension r(10),dr(1)

      rho=1.0d0/dfloat(nrp)
      ri=r(1)
      if(xi.ne.0.d0) ri=r(1)/dexp((xi)**rho)
      xmx=(dlog(rmx/ri))**nrp-xi
      h=xmx/dfloat(no-1)
      do 25 j=1,no
      x=(j-1)*h
   25 r(j)=ri*dexp((x+xi)**rho)
      do 26 j=1,no
   26 dr(j)=1.0d0
      if(nrp.eq.1) return
      do 27 j=1,no
   27 dr(j)=dfloat(nrp)*(dlog(r(j)/ri))**(nrp-1)
      return
      end
c --------------------------------------------------
      subroutine rsin(r,dr,rmx,h,no)
 
c --- r = rmax sin[(pi/2) (x/rmax)**y]
c --- dr = rdx/dR
c --- dR = [2 rmax/(pi*y)] [rmax/x]**(y-1) tan[(pi/2) (x/rmax)**y]
 
      implicit real*8(a-h,o-z)
      dimension r(10),dr(1)
      if(no .le. 1) stop
      h=rmx/dfloat(no-1)
      pi=3.14159265358979324d 00
      y=-dlog(2.0d0*dasin(r(2)/rmx)/pi)/dlog(dfloat(no-1))
      drmx=2.0d0*rmx**y/(pi*y)
      yy=y-1.0d0
      r(1)=0.0d0
      dr(1)=0.0d0
      do 10 j=2,no
      fh=dfloat(j-1)*h
      aug=(pi/2.0d0)*(fh/rmx)**y

      r(j)=rmx*dsin(aug)
   10 dr(j) = drmx * dtan(aug) / fh**yy

      return
      end
c ------------------------------------------------------------------
      function rint (f,na,nb,nq,h)
      implicit real*8(a-h,o-z)
c
c  this program calculates the integral of the function f from point na
c  to point nb using a nq points quadrature ( nq is any integer between
c  1 and 14 ).  h is the grid size.
c                                      written by c. c. j. roothaan
c
      dimension c(105),c1(25),c2(80),d(14),f(1)
      equivalence (c1(1),c(1)),(c2(1),c(26))
      data c1/1.d0,2.d0,1.d0,23.d0,28.d0,9.d0,25.d0,20.d0,31.d0,8.d0,
     &1413.d0,1586.d0,1104.d0,1902.d0,475.d0,1456.d0,1333.d0,1746.d0,
     &944.d0,1982.d0,459.d0,119585.d0,130936.d0,89437.d0,177984.d0/
      data c2/54851.d0,176648.d0,36799.d0,122175.d0,11108.d1,156451.d0,
     & 46912.d0,220509.d0,29336.d0,185153.d0,35584.d0,7200319.d0,
     & 7783754.d0,5095890.d0,12489922.d0,-1020160.d0,16263486.d0,
     &261166.d0,11532470.d0,2082753.d0,7305728.d0,6767167.d0,9516362.d0,
     & 1053138.d0,18554050.d0,-7084288.d0,20306238.d0,-1471442.d0,
     & 11965622.d0,2034625.d0,952327935.d0,1021256716.d0,636547389.d0,
     & 1942518504.d0,-1065220914.d0,3897945600.d0,-2145575886.d0,
     & 3373884696.d0,-454944189.d0,1637546484.d0,262747265.d0,
     & 963053825.d0,896771060.d0,1299041091.d0,-196805736.d0,
     & 3609224754.d0,-3398609664.d0,6231334350.d0,-3812282136.d0,
     & 4207237821.d0,-732728564.d0,1693103359.d0,257696640.d0,
     & 5206230892907.d0,5551687979302.d0,3283609164916.d0,
     & 12465244770050.d0,-13155015007785.d0,39022895874876.d0,
     & -41078125154304.d0,53315213499588.d0,-32865015189975.d0,
     & 28323664941310.d0,-5605325192308.d0,9535909891802.d0,
     & 1382741929621.d0,5252701747968.d0,4920175305323.d0,
     & 7268021504806.d0,-3009613761932.d0,28198302087170.d0,
     & -41474518178601.d0,76782233435964.d0,-78837462715392.d0,
     & 81634716670404.d0,-48598072507095.d0,34616887868158.d0,
     & -7321658717812.d0,9821965479386.d0,1360737653653.d0/
      data d/2.0d0,2.0d0,24.0d0,24.0d0,1440.0d0,1440.0d0,120960.0d0,
     &120960.0d0,7257600.d0,7257600.d0,958003200.d0,958003200.d0,
     &5230697472000.0d0,5230697472000.0d0/
      a=0.0d0
      l=na
      m=nb
      i=nq*(nq+1)/2
      do 1 j=1,nq
      a=a+c(i)*(f(l)+f(m))
      l=l+1
      m=m-1
    1 i=i-1
      a=a/d(nq)
      do 2 n=l,m
    2 a=a+f(n)
      rint=a*h
      return
      end
C############################################################



