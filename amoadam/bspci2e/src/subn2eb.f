C####################################################
      subroutine nom(n2e_bin, kl, e, l, no, r, p, z1, ethr, ad, vp)
      implicit real*8(a-h,o-z)
      include "parameter.2e.inc"      
      INCLUDE "units.inc"
      integer  n2e_asc 
      dimension del(np)
      dimension p(np),r(np)
C............................

C     in Ryd

      n2e_asc = n2e_bin + 1


      dk2 = E - Ethr


         w1 = 0.0D+00
         an = 1.0D+00        


         IF(dk2.GT.0.0D+00) THEN 

         ek = SQRT( dk2 )

         print*, "#               k_e = " ,  ek/2.0 


C#   Fitting to the last 10 points


C....................  start matching

C#     method No 1              (phi_1, w_1, a_1)
C#   Fitting to the last NP0 points with NP_STEP
C#   for getting the short-range scattering phase shift
C#   
C#   w1 = del
C#  


        NPOINTS = no
        NP0     = 10
        NP_STEP = 1
        NP1     = 5
        npoint_1 = npoints - np1

      IF( NP0 > NPOINTS.OR.NP1 > NPOINTS) THEN 
         WRITE(*,*) '# ERROR IN SELECTING MATCHING POINTS  '
         WRITE(*,*) '#  NP_STEP < NP1 < NP0 < NPOINTS      '
         WRITE(*,*) '#                            NPOINTS =', NPOINTS
         WRITE(*,*) '#                                NP0 =', NP0
         WRITE(*,*) '#                                NP1 =', NP1
         STOP    
      ENDIF

        id1 = 0
C        DO  j = no - 10, no-2,  1

       DO  j = npoints - np0, npoints - 2,  np_step 

           id1 = id1 + 1

           del(id1) = delta(r(j-1),r(j),p(j-1),p(j),l,ek,z1,ad,vp)
           aa2 = zeta( r(j), ek, l, z1, vp)

C     print 222,r(j), del(id1), aa2

        ENDDO

c# average value on the last id1 points
            
        w1 = 0.0D+00
        do  i = 1, id1
           w1 = w1 + del(i)
        enddo
       
C#       w1 : short range (scattering) phase shift, method 1

        w1 = w1/id1

C#       w2 : short range (scattering) phase shift,  method 2

C        write(*,*) l,z1,vp

        wr  = delta_r( kl, ek, l, z1, r(no))

        w2  = delta(R(NPOINT_1-1),R(NPOINT_1),P(NPOINT_1-1),P(NPOINT_1)
     &             ,L, EK,Z1,AD, VP )
 


!        function delta_R( kl, ek, l, zeff, r, r_b)

        ph1 = phi(  r(npoint_1), ek, l, z1, ad, vp)
        aa1 = zeta( r(npoint_1), ek, l, z1, vp)

C#   An is determined such a way as  :
C#                 ______
C#   A_k * P_k  =  /2/pi*k * sqrt( k / z(r) ) * sin[ phi + delta ] ===>
C#
C#                   ________  
C#                  /    2
C#              =  / -------- *   sin[ phi + delta ] ===>
C#                /  pi * z(r)
C#                      ___________                                  
C#                     /2/(pi*z(R)) * sin[ phi(R) + delta(R) ]      
C#        1./A_k  =  --------------------------------------- 
C#                              P_k(R)                           
C#
C#   note that    z(R) --> k
C#

        g1 = sqrt( 2.0D+00/ ( pi * aa1 ) ) * sin ( ph1 + w1 ) 
        an  = g1 / p( npoint_1 )

        g2 = sqrt( 2.0D+00/ ( pi * aa1 ) ) * sin ( ph1 + w2 ) 
        an2 = g2 / p( npoint_1 )

!
        if(w2.lt.0.0d+00) w2 = w2 + pi
!        w2 = kl*pi + w2
       
        write(n2e_asc, '(6(1PE20.8))')  0.5*ek, cos(w1), cos(w2), an,an2

!        write(nnormascii, '(6(1PE20.8))')  0.5*dk2, w2, mod(wr,pi)
!     &                                              ,w1,an2,an

C         write(*,*) '#           E, k^2, = ', E, dk2 
C         write(*,*) '#          w1,  w2  = ', w1, w2
C         write(*,*) '#          an, an2  = ', an,an2
  
      ENDIF
            
      write(n2e_bin)  e, w1, an

C      write(nnorm)  e, w2, an2


      return
      end
C#
C####################################################################
      function delta_r( kl, ek, l, zeff, r_b)
      implicit none
      include "units.inc"
      double precision ek,zeff,r_b,qz
      integer l,kl
      double precision phi_c_R, delta_R
      double precision pscoul
C....................

      phi_c_R = zeff * dlog( 2 * ek * r_b)/ek - dble(l) * pi/2.0d+00       
     &         + pscoul(l,zeff,ek)

      delta_R = kl * pi - ek * r_b  - phi_c_R


      end
    
C###########################################################
      function delta(r1, r2, p1, p2, l, ek, znuc, ad, vp)
      implicit real*8(a-h,o-z)
 
C..................
            ph1 = phi( r1, ek, l, znuc, ad, vp)
            aa1 = zeta(r1, ek, l, znuc, vp)

C            a1  = 1.d0/dsqrt(aa1)
C....
            ph2  = phi( r2, ek,  l, znuc, ad, vp)
            aa2  = zeta(r2, ek, l, znuc, vp)

C            a2   = 1.d0/dsqrt(aa2)

C..................

            A =  p2 * sin(ph1)/sqrt(aa1) 
            B =  p1 * sin(ph2)/sqrt(aa2)  

            C =  p1 * cos(ph2)/sqrt(aa2) 
            D =  p2 * cos(ph1)/sqrt(aa1)
            
            delta =  atan( ( A-B) / ( C - D ) )

c            delta=datan( ( a1 * dsin(ph1) * p2 - a2 * dsin(ph2) * p1 ) /
c     1     ( a2 * dcos(ph2) * p1 - a1 * dcos(ph1) * p2 ) )


            return
            end
c  this function caculate arg of gamm function
      function pscoul(l, qz, qk)
      implicit complex*16(a-h,r-z),real*8(o-q)

      q1=l+1
      ak=dcmplx(0.d0,qz/qk)
      ga=q1-ak
      pi=dasin(1.d0)*2.d0
      cc=cdlog(ga)+cdlog(1.d0+ga)+cdlog(2.d0+ga)+cdlog(3.d0+ga)
      cc=cc+cdlog(4.d0+ga)+cdlog(5.d0+ga)+cdlog(6.d0+ga)+cdlog(7.d0+ga)
      z=ga+8.d0
      gam=-cc+(z-0.5d0)*cdlog(z)-z+0.5d0*dlog(2.d0*pi)+1.d0/(12.d0*z)
      gam=gam-1.d0/(360.d0*z**3)+1.d0/(1260.d0*z**5)-1.d0/(1680.d0*z**7)
      pscoul=dimag(gam)
      return
      end
C*=======================================================================
C#     x -- r,  q -- momentum(k),  z -- charge at large distance,
C#     alpha -- polarization potential,  l -- angular momentum,
C#     zeta approaches q when x becomes very large

      function zeta(x, q, l, z, vp)
      double precision zeta, x, q, z, c, zd
      double precision w0,w1,w2
      double precision vp

      c      = dble(l*(l+1))/x**2
      zd     = z/x

      w0  = q**2 - c + 2.0D+00*zd - vp


      w1 = (4.0D+00/x**2) * ( c + 2.d0 * vp - zd)**2

C#
C# make again the analytical calculation to see what is the 
C# correct formula.
C#
C#      w1 = (2.0D+00/x) * ( c + 2.d0 * vp - zd)**2

      w2 = -(2.0D+00/x**2) * (3.0D+00*c + 10.0D+00*vp - 2.d0 * zd)

      zeta = sqrt(w0)+0.15625d0*w1*w0**(-2.5)-0.125d0*w2*w0**(-1.5)

      return
      end
C##################################################################

C   x  --  r  
C   q  -- momentum(k)  
C   z  -- charge at large distance
C   vp -- polarization potential  
C   l  -- angular momentum

      double precision function phi( x, q, l, z, ad, vp)
      INTEGER L

      DOUBLE PRECISION ad, dln, invcos 
      DOUBLE PRECISION x, q, z, vp, halfpi, c, qr, xi 
      DOUBLE PRECISION cc, fk, rho, qqc, theta, vptheta
      DOUBLE PRECISION  phi_r1, phi_r2, phi_r4
      DOUBLE PRECISION  pscoul
C....................

      qr = q * x

      halfpi = 1.5707963268D+00

      c  = dble( l * ( l + 1 )) + ad 

C$$$      WRITE(*,*) ' l    = ', l
C$$$      WRITE(*,*) ' ad   = ', ad
c$$$      WRITE(*,*) ' c    = ', c
c$$$      WRITE(*,*) ' z    = ', z
c$$$      WRITE(*,*) ' q-c  = ', qr**2 - c 

      IF(z.EQ.0.D+00) THEN 

C#
C#  negative ions :   Z_eff = 0.0D+00 
C#
C#           V(r) = - [ l(l+1) + a ] / r^2  - vp / r^4  + O( r^{-x} )
C#

C#****  NEGATIVE ION CASE  : H- , ...

         xi =  sqrt( qr**2 - c )

         IF(xi.LT.0.0D+00) THEN

            WRITE(*,*) ' (k*r)^2 -c  NEGATIVE '
            WRITE(*,*) '     k = ', q**2
            WRITE(*,*) '     r = ', x**2
            WRITE(*,*) ' l,a,c = ', l,a,c

            STOP
            
         END IF

Cxxxx   0.208333 ==  5/24
Cxxxx   0.125    ==  1/8

         IF(c.EQ.0.0D+00) THEN

            theta   = 0.125D+00 / qr

            vptheta = vp * 0.5D+00 / (q*q*x )
 

         ELSE IF(c.GT.0.0D+00) THEN

            cc = sqrt(c)

            theta   = ( ( c + 0.125D+00) / cc) * asin( cc / qr ) 

            vptheta = vp / ( x * ( xi + q*x) )  

         ELSE

            write(*,*) " c is negative. Is that ok? "
            stop

            cc = sqrt(-c)
            
            dln = log( ( cc  + xi ) / qr )

            theta   = ( ( c + 0.125D+00) / cc) * dln

            vptheta = vp / ( x * ( xi + q*x) )  


         ENDIF

C#  2*Z/r             contribution 

         phi_r1 = 0.D+00


C#   - [ l(l+1) + a ] / r^2          contribution 

         phi_r2 = xi - halfpi*dble(l) - 0.125D+00 / xi 
     1               - 0.208333333D+00 * c/ xi**3
     1               + theta

C#  -vp/r^4         contribution 


         phi_r4 =  0.5D+00 * vptheta            
           
         phi = phi_r1 + phi_r2 + phi_r4 

      ELSE 

C#
C#  potentials of type : 
C#
C#           V(r) = 2*Z/r - [ l(l+1) + a ] / r^2  - vp / r^4  + O( r^{-x} )
C#

         fk  = q / z
         rho = x * z
         xi  = sqrt( qr**2 + 2.0D+00 *rho - c )
         qqc = fk * fk * c


         IF(xi.LT.0.0D+00) THEN

            WRITE(*,*) ' (k*r)^2 -c  NEGATIVE '
            WRITE(*,*) '     fk  = ', qr**2
            WRITE(*,*) '     rho = ', rho**2
            WRITE(*,*) ' l,a,c   = ', l,a,c

            STOP
            
         END IF

         IF(c.EQ.0.0D+00) THEN

            theta   = 0.25D+00/ ( xi + qr )
            
            vptheta = vp *  2 * z * ( xi + 2 * fk * rho )
     1                   / ( 3* rho *(xi + fk * rho )**2 )

         ELSE IF(c.GT.0.0D+00) THEN

            cc = sqrt( c )

            invcos = acos( ( rho - c + fk * c *xi )
     1                  / ( (1.0D+00 + qqc) * rho ) )


            theta = ( ( c + 0.125D+00) / cc ) * invcos 
            
            vptheta = vp * (  - z * ( 2 * rho - c ) 
     1                   / ( c * rho *(xi + fk * rho ))
     1                   + ( z/ c * cc) * invcos  )

         ELSE

            cc = dsqrt( - c )

            dln = log( ( rho - c + xi * cc )
     1                / ( (1.0D+00 + fk *cc ) * rho ) )

            theta = ( ( c + 0.125D+00) / cc ) * dln 

            vptheta = vp * (  - z * ( 2 * rho - c )
     1                   / ( c * rho *(xi + fk * rho ))
     1               + (z/ c * cc) * dln  )
            
         ENDIF

!!!  2*Z/r             contribution 


         phi_r1 =  - ( 1.0D+00 + log(fk)) / fk + pscoul( l, z, q)


!!!   - [ l(l+1) + a ] / r^2          contribution 

         phi_r2 = xi 
     &            + log( 1.0D+00 + qr * fk + fk * xi )/fk 
     &            -  halfpi * dble(l)                                 
     &            - ( xi * ( 3.0D+00 * qqc + 4.0D+00 ) 
     &            + qr *( 3.0D+00 * fk*fk + 2.0D+00) + fk * c ) 
     &            / ( 24.0D+00 * (1.d0 + qqc) * xi * (xi + qr ) )   
     &            + 0.208333333D+00*( rho - c ) / xi**3 
     &            + theta

!!!  -vp/r^4         contribution 

              phi_r4 = 0.5D+00 * vptheta        

              phi = phi_r1 + phi_r2 + phi_r4 

           ENDIF

      
c$$$           WRITE(*,*)  ' Z/r            COULOMBIC PHASE = ', phi_r1
c$$$           WRITE(*,*)  ' [l(l+1) + a ]/r^2 CI     PHASE = ', phi_r2
c$$$           WRITE(*,*)  '  b/r^4      POLARIZATION PHASE = ', phi_r4
c$$$           WRITE(*,*)  '                    TOTAL PHASE = ', phi

         END 
    
C###########################################################
C#
      subroutine wfnlin(p, l, no, nmin, nmx)
      implicit real*8(a-h,o-z)
      include "parameter.2e.inc"

      dimension p(np)
      common/pca/pc(ns, np)
C...........................................


      NWF1E = 2

C..........................................

      call wf1efile(NWF1E, l) 

      kmx  = nmx - nmin + 1
      kmin = nmin - 1

      print *, no, kmin, kmx

      if(kmin.eq.0) go to 20

      DO  k = 1, kmin

      read(NWF1E) ( p(j), j = 1, no)

      ENDDO

 20   DO  k = 1, kmx

         READ(NWF1E) ( p(j), j = 1, no )

         DO  j = 1, no

            pc(k,j) = p(j)

         ENDDO

      ENDDO

      CLOSE(NWF1E)

      RETURN
      END
C#######################################################

c**** x -- r,  q -- momentum(k),  z -- charge at large distance,
c**** vp -- polarization potential,  l -- angular momentum

c$$$      function phi(x,q,z,vp,l)
c$$$      double precision phi,x,q,z,vp,pi,c,qr,xi,cc,fk,rho,qqc,theta,thi
c$$$      double precision pscoul
c$$$      qr=q*x
c$$$      pi=1.5707963268d0
c$$$      c=float(l*(l+1))
c$$$
c$$$      if (z .ne. 0.d0) go to 100
c$$$c**** v-n potential
c$$$      xi=dsqrt(qr**2-c)
c$$$      phi=xi-pi*float(l)-0.125d0/xi-(5.d0/24.d0)*c/xi**3
c$$$      if (l .gt. 0) go to 1
c$$$      phi=phi+0.125d0/qr+vp*x/(6.d0*q)
c$$$      return
c$$$    1 cc=dsqrt(c)
c$$$      phi=phi+(c+0.125d0)/cc*dasin(cc/qr)+0.5d0*((qr*x)**2/(2.d0*c)*
c$$$     1    vp/cc*dasin(cc/qr)-xi*x**2*vp/(2.d0*c))
c$$$      return
c$$$  100 continue
c$$$c**** v(n-z) potential
c$$$      fk=q/z
c$$$      rho=x*z
c$$$      xi=dsqrt(qr**2+2.d0*rho-c)
c$$$      qqc=fk*fk*c
c$$$      if(l .gt. 0) go to 2
c$$$      theta=0.25d0/(xi+qr)
c$$$      go to 3
c$$$    2 cc=dsqrt(c)
c$$$      theta=(c+0.125d0)/cc*dacos((rho-c+fk*c*xi)/(1.d0+qqc)/rho)
c$$$    3 phi=xi+dlog(1.d0+qr*fk+fk*xi)/fk-pi*float(l)-(xi*(3.d0*qqc+4.d0)
c$$$     1    +qr*(3.d0*qqc+2.d0)+fk*c)/(24.d0*(1.d0+qqc)*xi*(xi+qr))
c$$$     2    +(5.d0/24.d0)*(rho-c)/xi**3+theta
c$$$      phi=phi-(1.d0+dlog(fk))/fk+pscoul(l,z,q)
c$$$      if (vp .eq. 0.d0) return
c$$$      if (l .gt. 0) go to 4
c$$$      thi=(2.d0/15.d0)*vp*x**2*(9.d0*qr*x+11.d0*qr**2+6.d0*rho)/(xi+
c$$$     1    qr)**3
c$$$      go to 5
c$$$    4 thi=z**2*vp*x**4/(2.d0*c**2)*(3.d0+qqc)*dacos((rho-c+fk*c*xi)/
c$$$     1    (1.d0+qqc)/rho)/cc-z*vp*x**3/(2.d0*c**2)*(3.d0*(2.d0*rho-c)/
c$$$     2    (xi+qr)+c*xi/rho)
c$$$    5 phi=phi+0.5d0*thi
c$$$      return
c$$$      end
c$$$C=====================================================================
c$$$C*EOF











