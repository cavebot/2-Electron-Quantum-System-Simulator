C#######################################################################
C#
C#    input :
C#              n =  number of collacation points
C#              k =  order of spline > 2
C#            a,b =  1st and last points on interval
C#
C#
C#    output :
C#         t(1)   =  t(2) = ... t(k) = a
C#         t(n+1) =  t(n+2) = ... t(n+k)   
C#         t(i)   =  a + (i-k)*h   i = k+1,n (linear sampling)
C#            h   =  (b-a)/(n-k+1)
C#
C#######################################################################
C#
C#       The difference between  subroutine RIN(...) and MKGRID(,...)
C#       is the following :
C#
C#       mkgrid : Generates the knot -sequence t(i) , i = 1, 2,..., n+k
C#                on which the the diagonalization of the matrix equations 
C#                takes place.
C#
C#
C#                     ___ 
C#      p_{m,l}(r)  =  \     C^{n,l}_{i,k}  * B_{i,k}(r)          
C#                     /__ 
C#                      i
C#                          i = 1,2,....no of B-Splines set
C#                          k   Order of B-splines
C#                          m   Eigenstate |m,l >  of 
C#                              energy e_n anf momentum l 
C$                              out of n - 2 eigensolutions
C#      
C#     R_{m.l} = p_{ml}(r) / r    
C#
C#     Therefore the produced Coeficients C_{i,k} are calculated based 
C#               on the paramaters passed to mkgrid subroutine.
C#
C#     However: 
C#
C#     Having obtained the Coefficients C^m_{i,k} for the 
C#     wavefunction (wf)      p_{m.l}(r) 
C#     we are now able to calculate the wf not only on the knot 
C#     points t(i) generated from the mkgrid subroutine based on 
C#     paramaters (n,k,Rmax,knot-sequence) but everywhere in between 
C#     [0,Rmax].
C#      
C#  rin:
C#     This is achieved through the RIN(...) subroutine which can 
C#     ask for evaluating the wf in 'no' points with a selected 
C#     knot sequence.
C#     This is desired since usually is needed to calculated dipole 
C#     matrix elements (dme) between states which have been  calculated 
C#     with different grid parameters.
C#
C#     i.e. :
C#
C#     |1s H > :
C#               Exponential knot points n = 100 , k = 8, Rmax = 100
C#
C#     |2p H > : 
C#                Linear  Knot points,   n = 80,   k = 9, Rmax = 100
C#
C#     And then you can ask for the dme :
C#
C#                     _Rmax
C#                    /
C#     < 1s|r|2p > = /  dr  p_1s(r) * r * p_2p(r)
C#                 _/
C#                 0
C#    
C#     (Note that the Rmax should be equal in both case otherwise
C#      orthogonality  < 1s|2p> = 0 is not satisfied.)
C#
C#     You can ask, therefore for evaluation of the integral in 
C#     'no' of points with whatever knot you wish.
C#     The 'bvalue' subroutine (de' Boor) is used to evaluate the 
C#     expanded on B-splines set wavefuntion in a point  r  different
C#     than the  t(i).   ( t(i) != r for all i)
C#
C#
C#           _________
C#          /         \         f(t(i))         
C#         /           \                           ____________
C#        /             \_______|___              /            \      
C#       /                          \            /              \     |
C#     |/                            \          /                \____|
C#                                    \________/                     
C#                  
C#  r: 0                        t(i)              t(i+k)            t(n+k)  
C#     .... . . . . .  .   .    .   ,    .         .        .         .
C#                                  |
C#                                  |
C#                                  r
C#
C#
C#
C#######################################################################

      subroutine mkgrid(n, k, rmax, dr0, t)
      implicit real*8(a-h,o-z)
      dimension t(1)
C------------------------------------------------------------------

      open(10,file = 'INP/GRID.INP')

      read(10, '(I3)')  idbsp 
      read(10, '(E10.2,1X,I5)' )  xi,    nrp
      read(10, '(2E10.2,1X,I5)')  rmin,  rmx, no

      close(10)
      

      if(rmx.ne.rmax)  then 
         write(*,*) ' Incosistent Input for Radius box in bspeg.din
     1        and rin.din file'
         
         write(*,*) 'rmx(rin.din) = ',rmx, ' rmax(bspeg.din) =', rmax
         stop

      endif

C####      First   k   knot-points are zero 

      do i = 1, k

         t(i) = 0.0D0

      enddo

C###   Define some Sine - like   quantities

      if(idbsp.eq.0) then

C         write(*,*) 'Sine-like knots used for Diagonalization'

         Hh = RMax / DFLOAT(n - k + 1)
         PI = 3.14159265358979324D00
         Y  = -DLOG(2.0D0*DASIN(dr0/RMax)/PI)/DLOG(DFLOAT(n-k+1))

C         write(6,1000) ' t0 =',dr0,'  h =',hh


C###  Generate Sine - like Knot points
         do  i = k + 1, n
         
            FH   = DFLOAT(i-k) * Hh             
            AUG  = (PI/2.0D0) * (FH/Rmax)**Y   
            t(i) = RMax * DSIN(AUG)             

         enddo
C###  Generate Exp - like Knot points      
      else if (idBsp.eq.1) then

C         write(*,*) 'Exp-like knots used for Diagonalization'

         rho = 1./nrp 
         ri  = rmin
         xmx = (dlog(Rmax/ri))**nrp - xi
         h   = xmx / dfloat( n - k + 1)


C         write(*,*)  ri,xi,xmx,rmax,h,nrp,rho
         
C         write(6,1000) ' t0 =',dr0,'  h =',hh
      

         do  i = k + 1, n

            x    =  (i - k) * h
            t(i) = ri * dexp( (x + xi)**rho )
            
         enddo
c
c      else if (idBsp.eq.2) then
c
c
c         do  i = k + 1, n
c
c            x    =  (i - k) * h
c            t(i) = ri * dexp( ( x + xi )**rho )
c            
c         enddo
c
c
c

      else if (idBsp.gt.1) then
         
         write(*,*) 'No available knot sequence for idBsp = ', idBsp
         write(*,*) 'Supported idBsp = 0 (Sine) , idBsp = 1 (exp) '
         stop
         
      endif

 1000 format(/a,1pd14.6,a,1pd14.6/)


C####        Last   k   knot-points are zero 
      do i = n + 1, n + k

         t(i) = rmax

      enddo
             
      return
      end
C#######################################################################
C#   EOF
