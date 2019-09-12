
!
!    MODULE WF_1E : 
!
!    SUBROUTINES :
!                  READ_TARGET_WF
!                  P_ASYMPTOTIC
!                  PSCOUL
!                  FGWKB
!                  WKB
!
!      FUNCTIONS :
!                  DELTA
!                  PHI
!                  ZETA
!
MODULE WF_1E

  USE PRECISION, ONLY: DPK
!  IMPLICIT NONE
  PUBLIC

  CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE P_ASYMPTOTIC(l, k_e, zeff, r, p_a, dp_a)
  
  USE precision, ONLY: DPK
  USE units, ONLY: M_PI, M_PI_2

  IMPLICIT NONE

  INTEGER   L
  REAL(DPK) k_e, zeff, r
  REAL(DPK) AMP_COULOMB, DAMP_COULOMB
  REAL(DPK) PHASE_COULOMB
  REAL(DPK), INTENT(OUT) :: P_A, DP_A


!..............

  AMP_COULOMB  = SQRT( 2.0D+00/(M_PI*K_E) ) 
  DAMP_COULOMB = AMP_COULOMB *( K_E + ZEFF /(k_E*R) )

  PHASE_COULOMB =  ZEFF * LOG(2 * K_E * R) / K_E  + PSCOUL(L, ZEFF, K_E) 


  P_A  =  AMP_COULOMB * SIN( K_E * R - M_PI_2 * L  + PHASE_COULOMB )
  DP_A = DAMP_COULOMB * COS( K_E * R - M_PI_2 * L  + PHASE_COULOMB )


END SUBROUTINE p_asymptotic
!  CONTAINS



!
!  this function caculate arg of gamm function
!
FUNCTION pscoul(l, qz, qk)

  USE units, ONLY: M_PI
  
  IMPLICIT COMPLEX*16(a-h,r-z), REAL*8(o-q)
  INTEGER L
  
  q1 = l + 1
  ak = CMPLX(0.0D+00, qz / qk )
      
  ga = q1 - ak
      
  cc=LOG(ga) + LOG(1.d0+ga)+ LOG(2.d0+ga)+LOG(3.d0+ga)
  cc=cc+LOG(4.d0+ga)+ LOG(5.d0+ga)+LOG(6.d0+ga)+LOG(7.d0+ga)
  
  z = ga + 8.0D+00
  
  gam = -cc + (z - 0.5d0)* LOG(z)-z+0.5d0* LOG(2.d0*M_PI)+1.d0/(12.d0*z)
  gam = gam-1.d0/(360.d0*z**3)+1.d0/(1260.d0*z**5)-1.d0/(1680.d0*z**7)
  pscoul= AIMAG(gam)
  
  RETURN
END FUNCTION pscoul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!   
!  Gamma function for integer       G(z+1) = z * G(z)
!
!
REAL(DPK) FUNCTION GAMMA_I(N)
  
  IMPLICIT NONE 
  INTEGER N, NP, I
  REAL(DPK) GAMMA
  
  IF(N.LE.0) THEN 

     WRITE(*,*) '# FUNCTION GAMMA_I : ERROR '
     WRITE(*,*) '# NEGATIVE OR ZERO ARGUMENT PASSED'
     WRITE(*,*) '#                             N =', n
     STOP
  ENDIF

!     WRITE(*,*) ' NP, GAMMA_I', I, NP

  NP = 0
  GAMMA = 1
  DO I = 1, N

     NP = NP + 1


     GAMMA = GAMMA * NP

  ENDDO

  GAMMA_I = GAMMA
    RETURN 
END FUNCTION GAMMA_I
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


REAL(dpk) FUNCTION delta(r1, r2, p1, p2, l, ek, znuc, a, vp)

 !     USE units, ONLY: DPK

      IMPLICIT NONE
      
      INTEGER L
      REAL(DPK) ZNUC, A, VP, EK
      REAL(DPK) AA1, AA2,A1,A2
      REAL(DPK) R1, R2, P1, P2
      REAL(DPK) AA, BB, CC, DD
      REAL(DPK) PH1, PH2

!................      
!                  1
!    P_ks(r) = ----------- *  sin( phi(r) )    
!              sqrt( z(r) ) 
!                  1
!    P_kc(r) = ----------- *  cos( phi(r) ) 
!              sqrt( z(r) ) 
!
!   P_k(r_1 ) =  a_k * P_ks(r_1) + b_k P_kc(r_1) 
!   P_k(r_2 ) =  a_k * P_ks(r_2) + b_k P_kc(r_2)
!
!   
!   P_k(r)  calculated function at point r
!
!                   b_k 
!   tan( delta ) = ----- =
!                   a_k
!
!                   P_k(r_1) * P_ks(r_2) - P_k(r_2) * P_ks(r_1) 
!                = -------------------------------------------
!                   P_k(r_2) * P_kc(r_2) - P_k(r_1) * P_kc(r_2)
!................

      ph1 = phi( r1, ek, l, znuc, a, vp)
      aa1 = zeta(r1, ek, l, znuc, vp)
      a1  = 1.0D+00/SQRT(aa1)
      
      ph2  = phi( r2, ek,  l, znuc, a, vp)   ! coulombing
      aa2  = zeta(r2, ek, l, znuc, vp)       !   
      a2   = 1.0D+00/SQRT(aa2)
           
      AA = a1 * SIN(ph1) * p2 
      BB = a2 * SIN(ph2) * p1 
      
      CC = a2 * COS(ph2) * p1
      DD = a1 * COS(ph1) * p2 
      
      delta = atan( ( AA - BB ) / ( CC - DD ) )


      !            delta=datan( ( a1 * dsin(ph1) * p2 - a2 * dsin(ph2) * p1 ) /
      !     1     ( a2 * dcos(ph2) * p1 - a1 * dcos(ph1) * p2 ) )
      

      RETURN
    END FUNCTION delta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!#     x     --> r,  
!#     q     --> momentum(k),
!#     z     --> charge at large distance,
!#     alpha --> polarization potential,
!#     l     --> angular momentum,
!#     zeta  --> q when x ---> very large


    REAL(DPK) FUNCTION zeta(x, q, l, z, vp)

!      USE precision, ONLY: DPK
      IMPLICIT NONE
           
      INTEGER L
      REAL(DPK) x, q, z, c, zd
      REAL(DPK) omega, omega1, omega2
      REAL(DPK) VP
!.....................

      
         c      = DBLE( l * ( l + 1 ) ) / x**2
         zd     = z / x
      

         omega  = q**2 - c + 2.0D+00 * zd - vp
           
         omega1 =   (4.0D+00/x**2) * ( c + 2.0D+00 * vp - zd)**2
         omega2 = - (2.0D+00/x**2) * (3.0D+00*c + 10.0D+00*vp - 2.0D+00 * zd)
           
         zeta  =    SQRT(omega) &
              &    + 0.15625D+00 * omega1 * omega**(-2.5) &
              &    - 0.125D+00   * omega2 * omega**(-1.5)
      

      RETURN
    END FUNCTION zeta
!############################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!#     x     --> r,  
!#     q     --> momentum(k),
!#     z     --> charge at large distance,
!#     alpha --> polarization potential,
!#     l     --> angular momentum,
!#     zeta  --> q when x ---> very large


!    REAL(DPK) FUNCTION zeta_bound(n, l, z)
      !      USE precision, ONLY: DPK
!      IMPLICIT NONE
           
!      INTEGER   n, l
!      REAL(DPK) z
!.....................

!      WRITE(*,*)" G1,G2 = " , Gamma_I(n+l+1), Gamma_I(n-l) 
!
!      zeta_bound = n**2 * Gamma_I(n+l+1) * Gamma_I(n-l) 
!      
!      RETURN
!    END FUNCTION zeta_bound
!############################################################
!   x  --  r  
!   q  -- momentum(k)  
!   z  -- charge at large distance
!   vp -- polarization potential  
!   l  -- angular momentum

    
    REAL(dpk) FUNCTION phi( x, q, l, z, a, vp)

!      USE PRECISION, ONLY: DPK 
      USE UNITS, ONLY :M_PI_2
      IMPLICIT NONE

      INTEGER L

      REAL(DPK) a, dln, invcos 
      REAL(DPK) x, q, z, vp, c, qr, xi 
      REAL(DPK) cc, fk, rho, qqc, theta, vptheta
      REAL(DPK) phi_r1, phi_r2, phi_r4

           !....................
           
      qr = q * x
      c  = dble( l * ( l + 1 )) + a 

      IF(z.EQ.0.D+00) THEN 

!#
!#  potentials of type :   (Z_eff = 0.0D+00) 
!#
!#           V(r) = - [ l(l+1) + a ] / r^2  - vp / r^4  + O( r^{-x} )
!#

!#****  NEGATIVE ION CASE  : H- , ...

         xi =  sqrt( qr**2 - c )

         IF(xi.LT.0.0D+00) THEN

            WRITE(*,*) ' (k*r)^2 -c  NEGATIVE '
            WRITE(*,*) '     k = ', q**2
            WRITE(*,*) '     r = ', x**2
            WRITE(*,*) ' l,a,c = ', l,a,c

            STOP
            
         END IF

!xxxx   0.208333 ==  5/24
!xxxx   0.125    ==  1/8

         IF(c.EQ.0.0D+00) THEN

            theta   = 0.125D+00 / qr
            vptheta = vp * 0.5D+00 / (q*q*x )
 
         ELSE IF(c.GT.0.0D+00) THEN

            cc      = sqrt(c)
            theta   = ( ( c + 0.125D+00) / cc) * asin( cc / qr ) 
            vptheta = vp / ( x * ( xi + q*x) )  

         ELSE

            cc = sqrt(-c)
            
            dln = LOG( ( cc  + xi ) / qr )

            theta   = ( ( c + 0.125D+00) / cc) * dln

            vptheta = vp / ( x * ( xi + q*x) )  


         ENDIF

!#  2*Z/r             contribution 

         phi_r1 = 0.0D+00


!#   - [ l(l+1) + a ] / r^2          contribution 

         phi_r2 = xi - M_PI_2*DBLE(l) - 0.125D+00 / xi  &
              &      - 0.208333333D+00 * c / xi**3      &
              &      + theta

!#  -vp/r^4         contribution 


         phi_r4 =  0.5D+00 * vptheta            
           
         phi = phi_r1 + phi_r2 + phi_r4 

      ELSE 

!#
!#  potentials of type : 
!#
!#           V(r) = 2*Z/r - [ l(l+1) + a ] / r^2  - vp / r^4  + O( r^{-x} )
!#

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
            
            vptheta = vp *  2 * z * ( xi + 2 * fk * rho )&
                 &             / ( 3* rho *(xi + fk * rho )**2 )

         ELSE IF(c.GT.0.0D+00) THEN

            cc = sqrt( c )

            invcos = acos( ( rho - c + fk * c *xi ) / ( (1.0D+00 + qqc) * rho ) )


            theta = ( ( c + 0.125D+00) / cc ) * invcos 
            
            vptheta = vp * (  - z * ( 2 * rho - c ) / ( c * rho *(xi + fk * rho )) & 
     &                   + ( z/ c * cc) * invcos  )

         ELSE

            cc = sqrt( - c )

            dln = LOG( ( rho - c + xi * cc )              &
     &                / ( (1.0D+00 + fk *cc ) * rho ) )

            theta = ( ( c + 0.125D+00) / cc ) * dln 

            vptheta = vp * (  - z * ( 2 * rho - c )      &
     &                   / ( c * rho *(xi + fk * rho ))  &
     &                  + ( z/ c * cc) * dln  )
         ENDIF

!         fk  = q / z
!         rho = x * z
!         qr = q * x
!         c  = dble( l * ( l + 1 )) + a 
!         xi  = sqrt( qr**2 + 2.0D+00 *rho - c )
!         qqc = fk * fk * c

!!!  2*Z/r             contribution 


         phi_r1 =  - ( 1.0D+00 + log(fk)) / fk + pscoul( l, z, q)

!!!   - [ l(l+1) + a ] / r^2          contribution 

         phi_r2 = xi                                               &
     &            + log( 1.0D+00 + qr * fk + fk * xi )/fk          & 
     &            -  M_PI_2 * dble(l)                              &    
     &            - ( xi * ( 3.0D+00 * qqc + 4.0D+00 )             & 
     &            + qr *( 3.0D+00 * fk*fk + 2.0D+00) + fk * c )    & 
     &            / ( 24.0D+00 * (1.d0 + qqc) * xi * (xi + qr ) )  & 
     &            + 0.208333333D+00*( rho - c ) / xi**3            &
     &            + theta    

!!!  -vp/r^4         contribution 

              phi_r4 = 0.5D+00 * vptheta

              phi = phi_r1 + phi_r2 + phi_r4 

           ENDIF

      
!           WRITE(*,*)  ' Z/r            COULOMBIC PHASE = ', phi_r1
!           WRITE(*,*)  ' [l(l+1) + a ]/r^2 CI     PHASE = ', phi_r2
!           WRITE(*,*)  '  b/r^4      POLARIZATION PHASE = ', phi_r4
!           WRITE(*,*)  '                    TOTAL PHASE = ', phi

         END FUNCTION phi
!###########################################################

!*     -----------------------------------------------------
!*           F G W K B
!*     -----------------------------------------------------
!*
!*	This routine determines the energy normalized 
!*	Coulomb functions F , and G,
!*	and their derivatives F' , G' , with respect to r
!*
!*	Normalization obtained using the  WKB method as 
!*	proposed by Liu, Xi, and Li, PRA,48, 228,(1993)
!*
!*       Written my Jinhua Xi, December, 1994
!*     ----------------------------------------------------
!*
!*	F  =  sqrt(2/pi*k)  * sin(phi)
!*

!C  wkb

 SUBROUTINE fgwkb(ek, z, l, r, yr, dyr, f, df, g, dg, ierr)

   IMPLICIT DOUBLE PRECISION (a-h,o-z)
   DATA pi,pii/3.141592653589793d0,6.283185307179586d0/
  
  ierr = 0
      
  CALL  wkb(ek,z,l,r,zeta_wkb,dz,deltaz)


  IF( dabs(deltaz) .GT. 1.d-4 ) ierr = 1
  
  IF( dabs(deltaz) .GT. 1.d-1) THEN

     ierr = 2
     RETURN
  ENDIF

!C*
!C*	determine the phase function phi(r) and the
!C*	normalization constant 
!C*

  dzz = 0.5d0*dz/zeta_wkb
  pn  = dsqrt(2.d0/pi/zeta_wkb)

!C*

  phi_wkb = ATAN( zeta_wkb/( dyr / yr + dzz) )

  IF( SIN(phi_wkb)* yr.LT.0.d0 ) phi_wkb = phi_wkb + pi 


  IF( phi_wkb.LT.0.d0 )  phi_wkb = phi_wkb + pii
  IF( phi_wkb.GT.pii  )  phi_wkb = phi_wkb - pii  
!c
!c	for Coulomb potential , get the F, G, F',G'
!c
  f  =   pn * dsin(phi_wkb)
  g  =   pn * dcos(phi_wkb)

  df =   zeta_wkb * g - dzz * f
  dg = - zeta_wkb * f - dzz * g

  RETURN 

END SUBROUTINE fgwkb

!!$C#######################################################################
!!$*
!!$*     ------------------------------------------------------------------
!!$*           W K B
!!$*     ------------------------------------------------------------------
!!$*
!!$*	This routine performs the WKB iteration proposed by 
!!$*	Liu, Xi, and Li, PRA48, 228(1993)
!!$*
!!$*       The exact formulas for the iterative procedure were
!!$*       derived by C. F. Fischer, using the MAPLE Symbol 
!!$*       manipulation package.
!!$*	
!!$*	This routine is called by asympn and fgwkb.
!!$*	asympn: determines phase and normalization,
!!$*		as proposed by LXL paper
!!$*	fgwkb:  computes only f,g, f',g' 
!!$*
!!$*       Written by C. F. Fischer, July, 1994
!!$*       Modified by Jinhua Xi, December, 1994
!!$*     ------------------------------------------------------------------

SUBROUTINE wkb(ek, z, l, r, zeta_wkb, dz, deltaz)
  
  IMPLICIT DOUBLE PRECISION (a-h,o-z)
  DOUBLE PRECISION w(0:8), u(0:8,0:4)

!!$cxi
!!$cxi	it is ABSOLUTELY necessary to initiate the arrays
!!$cxi	if you would like to have a correct result. 
!!$cxi	

          
  DO  i = 0, 8

     w(i)= 0.0d+00
     
     DO  j = 0,4 

        u(i,j)= 0.0d+00

     ENDDO
  ENDDO


!cxi
!cxi	...................................................
!cxi

  a    = 2*z/r
  b    = -l*(l+1)/r/r
  ekk  = ek*ek
  w(0) = ekk + a + b
  
  IF( w(0) .LT. 0.3d0*ekk ) THEN

     WRITE(*,*) '  in WKB: r-value too small, r=', r
     deltaz=99.d0
     RETURN
  ENDIF

  u(0,0) = w(0)

  DO  i = 1,8

     a      = -i*a/r
     b      = -(i+1)*b/r
     w(i)   = a+b
     u(i,0) = w(i)/w(0)
  ENDDO


  DO  j = 0,3

     u(0,j+1) = w(0) + (5*u(1,j)**2 -4*u(2,j))/16.d0

     DO  i = 1, 6-2*j

        IF (i .EQ. 1) THEN

           u(1,j+1) = 7*u(1,j)*u(2,j)-5*u(1,j)**3 -2*u(3,j)

        ELSE IF (i .EQ. 2) THEN

           u(2,j+1) = 7*u(2,j)**2 -29*u(1,j)**2*u(2,j)            &
          &         + 9*u(1,j)*u(3,j) +15*u(1,j)**4 -2*u(4,j)

        ELSE IF (i .EQ. 3) THEN

           u(3,j+1) = 23*u(2,j)*u(3,j) -72*u(2,j)**2*u(1,j)       &
                &   +  147*u(1,j)**3*u(2,j) - 47*u(1,j)**3*u(2,j) &
                &  +  11*u(1,j)*u(4,j) - 15*u(1,j)**5 -2*u(5,j)

        ELSE IF (i .EQ. 4) THEN

           u(4,j+1) = 23*u(3,j)**2 -284*u(2,j)*u(3,j)*u(1,j)     &  
                &                 + 34*u(2,j)*u(4,j) + 657*(u(2,j)*u(1,j))**2   &
                &                 - 72*u(2,j)**3 - 888*u(1,j)**4*u(2,j)         &
                &                 + 288*u(1,j)**3*u(3,j) -69*u(1,j)**2*u(4,j)   &
                &                 + 13*u(1,j)*u(5,j) + 300*u(1,j)**6 - 2*u(6,j)
           
        ELSE IF (i .EQ. 5) THEN
           
           u(5,j+1) = -2*u(7,j) + 3030*u(2,j)*u(3,j)*u(1,j)**2       &
                &                   - 490*u(2,j)*u(4,j)*u(1,j)                      &  
                &                   - 500*u(2,j)**2*u(3,j) +47*u(2,j)*u(5,j)        & 
                &                   - 2040*u(1,j)**4*u(3,j) +495*u(1,3)**3*u(4,j)   & 
                &                   - 95*u(1,j)**2*u(5,j) +15*u(1,j)*u(6,j)         &
                &                   - 1800*u(1,j)**7 +80*u(3,j)*u(4,j)              & 
                &                   - 330*u(3,j)**2*u(1,j)                          &   
                &                   - 6180*u(2,j)**2*u(1,j)**3                      &   
                &                   + 1530*u(2,j)**3*u(1,j) + 6240*u(1,j)**5*u(2,j) 
           
        ELSE IF (i .EQ. 6) THEN
           
           u(6,j+1) = -2*u(8,j) + 62100*u(2,j)**2*u(1,j)**4           &
                &                   - 24660*u(2,j)**3*u(1,j)**2                      & 
                &                   + 4020*u(3,j)**2*u(1,j)**2                       &  
                &                   - 32640*u(2,j)*u(3,j)*u(1,j)**3                  &
                &                   + 5985*u(2,j)*u(4,j)*u(1,j)**2                   &
                &                   + 12150*u(2,j)**2*u(3,j)*u(1,j)                  &
                &                   - 1310*u(3,j)*u(4,j)*u(1,j)                      &
                &                   - 774*u(2,j)*u(5,j)*u(1,j)                       &
                &                   - 990*u(2,j)**2*u(4,j) -1330*u(2,j)*u(3,j)**2    &
                &                   + 16440*u(1,j)**5*u(3,j) +17*u(7,j)*u(1,j)       &
                &                   + 127*u(3,j)*u(5,j) +62*u(2,j)*u(6,j)            &
                &                   - 4020*u(1,j)**4*u(1,4) + 780*u(1,j)**3*u(5,j)   &  
                &                   - 125*u(1,j)**2*u(6,j) -50040*u(1,j)**6*u(2,j)   & 
                &                   + 12600*u(1,j)**8 +80*u(4,j)**2 + 1530*u(2,j)**4 
        END IF
        
        u(i,j+1) = (w(i) + u(i,j+1)/8.d0)/u(0,j+1)
        
     ENDDO
  ENDDO
  
  !.......
  
  IF( u(0,3) .LE. 0.d0 .OR. u(0,4) .LE. 0.d0 ) THEN
     
     WRITE(*,*) '  in WKB: r-value too small, r=', r
     
     deltaz =99.0D+00
     RETURN
     
  ENDIF
  
  !.......
  
  zeta_wkb = SQRT( u(0,4) )
  
  dz = (w(1)+(7*u(1,3)*u(2,3)-5*u(1,3)**3-2*u(3,2))/8)/(2*zeta_wkb)
  
  deltaz = zeta_wkb - SQRT( u(0,3) )
  
  RETURN 
  
END SUBROUTINE WKB
!####################################################

END MODULE WF_1E
!!!##########################################################
!!!EOF
             
