C#######################################################################
C#      subroutines : rin, rsin, rexp, rint
C#          
C#     rin    : reads "rin.din" input and grid parameters
C#              (nitialize grid)
C#     rsin   : sine - like  knot sequence
C#     rexp   : exp  - like knot sequence           
C#     rint   : performs the integral of a function 
C#                f(x) inside the  interval [a,b]
C#              
C#
C# 
C#     Authors    : X. Tang (1985 -1990) Jian Zhang (1993 -1996)
C#     
C#  
C#    Modified by : L.AA. Nikolopoulos                (Sept 1999)
C#
C#
C#
C#    Purpose :: reads "rin.din" input , 
C#               selects knot sequence  (sine or exp)  and then
C#               reads grid parameters
C#
C#
C#    Used by :: modules.f, norm.f ,..(?)
C#
C#
C#                                               iesl - days
C#######################################################################


C#######################################################################
C#                "rin.din   input File " 
C#
C#######################################################################
C#      idwf,idBsp      :: select the knot sequence
C#        
C#      idwf == 0          sine-like knot points
C#      idwf == 1          exp-like  knot points
C#      otherwise        stop
C#       
C#      xi, rmin, rmx, no ::  
C#                   
C#          For the interval [a,b] :
C#                  xi   =  a, 
C#                  rmx  =  b,     
C#                  rmin =  First non-zero knot point,
C#                  no   =  number of points where the 
C#                          function is going to be evaluated
C#
C#
C#      note :: 
C#             nrp input is used only in case of exp-like grid
C#                 see   "rexp" function
C#
C#######################################################################

      subroutine rin(r, dr, h ,no, idr)
      parameter(np = 2000)
      implicit real*8(a-h,o-z)
      dimension r(np),dr(np)
C....................................

      open(10,file='INP/GRID.INP')

      read(10,  '(I3)' )           idr  
      read(10,  '(E10.2,1X,I5)')   xi,   nrp
      read(10,  '(2E10.2,1X,I5)' ) rmin, rmx, no

      close(10)

C.....................................
C#    idBsp == 0     use sine-like knot sequence     (rsin) 
C#
C#    idwf  == 1      use exponent-like knot sequence (rexp)
C#                 
C#   otherwise 
C#                   stop


      if(idr.eq.0) then

         R(2) = rmin

         call rsin( r, dr, rmx, h, no)

C         write(*,*) 'Sine - like knot for Wavefunction'
c         write(*,3) r(2),h,rmx,no

         write(16, 3) r(2), h, rmx, no

      else if (idr.eq.1) then

         r(1) = rmin
         call rexp(r,dr,xi,rmx,h,no,nrp)
         
C         write(*,*) ' Exp - like knot for Wavefunction'
c         write(*,12) xi,r(1),rmx,h,no,nrp
         write(16,12) xi,r(1),rmx,h,no,nrp
         
      else if(idr.gt.1) then
      
         write(*,*) 'No available knot sequence for idwf = ', idwf
         write(*,*) 'Supported idwf = 0 (Sine) , idwf = 1 (exp) '
         stop
         
      endif
      
C........................................

    1 format(2x,'r(2) & rmx=?'/)
    2 format(2x,'no=?'/)
    3 format(2x,'SINE-LIKE GRID:  R(2), H, R & NP =', 1p3e14.6, i5/)
    7 format(2x,'xi,r(1) & rmx=?'/)
    8 format(2x,'no & nrp=?'/)
    9 format(2x,'r in sin scale: type 0 & r in exp scale: type 1'/)
   12 format(2x,'xi,r(1),rmx & h=',1p4e13.5,', no & nrp=',2i5/)
C.......................................


      return

      end

C#######################################################################
C#
C#    r = ri*dexp( (x+xi)**(1/nrp))
C#    dr = dx  r/[nrp * (ln(r/ri))**(nrp-1)]
C#
C#
C#
C#######################################################################

      subroutine rexp(r,dr,xi,rmx,h,no,nrp)
      parameter(np = 2000)
      implicit real*8(a-h,o-z)
      dimension r(np),dr(np)
C.................................
   
      rho = 1.0d0 / dfloat(nrp)
      ri  = r(1)

      if(xi.ne.0.d0)   ri = r(1) / dexp((xi)**rho)

      xmx = (dlog(rmx/ri))**nrp - xi
      h   = xmx / dfloat(no-1)

      do j = 1,no
         
         x    =  (j - 1) * h
         r(j) = ri * dexp((x + xi)**rho)

      enddo

      do  j = 1, no

         dr(j)= 1.0D+00

      enddo

      if(nrp.eq.1) return

      do  j = 1,no

         dr(j) = dfloat(nrp) * ( dlog( r(j)/ri ) )**(nrp-1)

      enddo

      return
      end
C#######################################################################
C#
C#
C#    r = rmax sin[(pi/2) (x/rmax)**y]
C#    dr = rdx/dR
C#    dR = [2 rmax/(pi*y)] [rmax/x]**(y-1) tan[(pi/2) (x/rmax)**y]
C#
C#
C#######################################################################

      subroutine rsin(r, dr, rmx, h, no)
      parameter(np = 2000) 
      implicit real*8(a-h,o-z)
      dimension r(np),dr(np)
C........................................................
 
      if(no .le. 1) stop

      h     =  rmx/dfloat(no-1)
      pi    =  3.14159265358979324d00
      y= -dlog(2.0d0*dasin(r(2)/rmx)/pi)/dlog(dfloat(no-1))
      drmx  =  2.0d0 * rmx**y / (pi*y)
      yy    =  y - 1.0d0
      r(1)  =  0.0D+00
      dr(1) =  0.0D+00

      do  j  =  2, no

         fh    =  dfloat(j-1) * h
         aug   =  (pi/2.0d0)  * (fh/rmx)**y
         r(j)  =  rmx * dsin(aug)
        dr(j)  =  drmx * dtan(aug) / fh**yy

      enddo
 
      return
      end
C#######################################################################
C#
C#
C#  this program calculates the integral of the function f 
C#  from point na to point nb using a nq points quadrature 
C#  ( nq is any integer between  1 and 14 ).  
C#   h is the grid size.
C#                                 Written by c. c. j. Roothaan
C#
C#
C#######################################################################
      function rint (f,na,nb,nq,h)
      implicit real*8(a-h,o-z)
      dimension c(105),c1(25),c2(80),d(14),f(1)
      equivalence (c1(1),c(1)),(c2(1),c(26))
C------------------------------------------------------------------


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

      a  = 0.0d0
      l  = na
      m  = nb
      i  = nq*(nq+1)/2

      do 1 j=1,nq
      a = a + c(i)*( f(l) + f(m) )
      l = l + 1
      m = m - 1
    1 i = i - 1

      a = a/d(nq)

      do 2 n = l, m
    2 a = a + f(n)

      rint = a * h

      return
      end

C#######################################################################
C#  EOF


