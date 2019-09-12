C**  mkgrid   
C**  step     
C#   rint
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
C#  
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
C##
C#eof


