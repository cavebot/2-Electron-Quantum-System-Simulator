!
! modules: set_grid, one_e_matrix, wf
!
!       set_grid: 
!             subroutines: rin, rsin, rexp, rexplin, cal_rpw
!
!   one_e_matrix:
!             subroutines: mkgrid, mkrgidmixed, setdipole, setmat, 
!    
!
!     potentials:
!            subroutines: init_hydrogenic_orbitals, corewf
!            functions  : hydrogenic, morse, woods_saxon, hwf
!
!       data: 
!       subroutines: write_v, write_mx, write_v_mx, read_v_mx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE functions
  !
  USE PRECISION, ONLY: dpk
  USE UNITS,     ONLY: M_PI, ZI
  !
  IMPLICIT NONE
  !
CONTAINS

!
! THE GAUSS DISTRIBUTION
!
  COMPLEX(dpk) FUNCTION fx_gauss(x,t)          ! cosine
    IMPLICIT NONE
    REAL(dpk),    INTENT(in) :: x
    REAL(dpk),      OPTIONAL :: t
    !locals
    REAL(dpk) k0
    REAL(dpk) x0,xt
    REAL(dpk) s0,st
    REAL(dpk) a0,at
    !
    k0 = 10.0D+00*M_PI
!    k0 = 0.0d+00
    
    x0 = 1.25D+00
    s0 = 0.125D+00
    a0 = 1.0D+00

    IF(.NOT.PRESENT(t))  t = 0.0D+00

    st =  SQRT( s0**2 + (t/s0)**2)
    at =  a0 * SQRT(s0/st)
    xt = x0 + k0*t

    fx_gauss = at * EXP( zi*k0*x - ( (x - xt) /(2.*st) ) )**2  !analytical 

  END FUNCTION fx_gauss
!
! real version
!
  REAL(dpk) FUNCTION fx_gauss_real(x)          ! cosine
    IMPLICIT NONE
    REAL(dpk),  INTENT(in) :: x
    REAL(dpk) :: t
    !locals
    REAL(dpk) k0
    REAL(dpk) x0,xt
    REAL(dpk) s0,st
    REAL(dpk) a0,at
    !
    k0 = 10*M_PI
!    k0 = 0.0d+00
    x0 = 1.25D+00
    s0 = 0.125D+00
    a0 = 1.0D+00

    !
    t =  0.0D+00
    !    IF(.NOT.PRESENT(t))  t = 0.0D+00

    st =  SQRT( s0**2 + (t/s0)**2)
    at =  a0 * SQRT(s0/st)
    xt = x0 + k0*t

!    fx_gauss_real = at * COS(k0*x)* EXP(- ( (x - xt) /(2.*st) ) )**2  !analytical 

    fx_gauss_real = COS(k0*x)* EXP( - ( x - x0 )**2.0/ (2.*s0**2) )
  END FUNCTION fx_gauss_real
  !
  !
  !
  REAL(dpk) FUNCTION fx_gauss_imag(x)           ! sine
    IMPLICIT NONE
    REAL(dpk),   INTENT(in) :: x
    REAL(dpk)               :: t
    !locals
    REAL(dpk) k0
    REAL(dpk) x0,xt
    REAL(dpk) s0,st
    REAL(dpk) a0,at
    !
    k0 = 10.0*M_PI
!    k0 = 0.0d+00
    x0 = 1.25D+00
    s0 = 0.125D+00
    a0 = 1.0D+00
    !
    t = 0.0D+00
    !    IF(.NOT.PRESENT(t))  

    st =  SQRT( s0**2 + (t/s0)**2)
    at =  a0 * SQRT(s0/st)
    xt = x0 + k0*t

!    fx_gauss_imag = at * sin(k0*x)* EXP(- ( (x - xt) /(2.*st) ) )**2  !analytical 

    fx_gauss_imag = SIN( k0*x )*EXP( - ( x - x0 )**2.0/ (2.*s0**2) )
  END FUNCTION fx_gauss_imag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! THE COULOMB POTENTIAL
!
  REAL(dpk) FUNCTION coulomb(x,z,l)
    IMPLICIT NONE
    REAL(dpk), INTENT(in) :: x,z,l
    coulomb = -z/x +  l*(l+1)/(2.0_dpk*x*x)
  END FUNCTION coulomb
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! THE HARMONIC OSCILLATOR
!
      REAL(dpk) FUNCTION harmonic(x,harm,l)
        IMPLICIT NONE
        REAL(dpk), INTENT(in) :: x,harm,l
        harmonic = harm*x*x + l*(l+1)/(2.0D+00*x*x)
      END FUNCTION harmonic
      !
      !
      ! A MORSE POTENTIAL
      !
      !
      REAL(dpk) FUNCTION morse(x,a,b)
        !
        REAL(dpk), INTENT(in) :: x,a,b
        REAL(dpk)             :: test
        !
        IF ((4.0D+00*(x-3.0D+00)).GT.690D+00) THEN
           morse = 0.0D+00
        ELSE
           morse = 25.0D+00*(EXP(-4.0D+00*(x - 3.0D+00) ) -&
               2.0D+00*EXP(-2.0D+00*(x-3.0D+00)) )
        ENDIF
      END FUNCTION morse
      
    END MODULE functions
    !M
    !M
    !M
    !M
    !M    
    MODULE set_grid
      !
      USE PRECISION, ONLY:dpk
      USE param  
      USE UNITS,    ONLY:M_PI 
      !
      IMPLICIT NONE
      !
      PUBLIC
      REAL(dpk)                :: hh
      REAL(dpk), DIMENSION(np) :: r, dr
      !
      INTERFACE r_grid
         MODULE PROCEDURE rin
      END INTERFACE r_grid

      PRIVATE rin, rsin, rlin, rexp, rlog

    CONTAINS
      !
      !##################################################
      !
      SUBROUTINE rin
        !
        IMPLICIT NONE
        !
        INTEGER :: i, nfile
        !
        
        SELECT CASE(idbsp)
        CASE(0)
           CALL rsin
        CASE(1)
           CALL rlin
        CASE(2)
           CALL rexp
        CASE(3)
           CALL rlog
        END SELECT
        

        ! store grid
        nfile = 12
        OPEN(nfile, file='out/grid.out')
        DO i = 1, SIZE(r)
           WRITE(nfile,*) i, r(i), dr(i), SIZE(r)
        ENDDO
        CLOSE(nfile)


      END SUBROUTINE rin

      !
      ! r = rmax sin[(pi/2) (x/rmax)**y]
      ! dr = rdx/dr
      ! dr = [2 rmax/(pi*y)]*[rmax/x]**(y-1) tan[(pi/2) (x/rmax)**y]
      !
      !      

        !  idbsp = 0
        !

      SUBROUTINE rsin
        !
        IMPLICIT NONE
        !LOC!
        INTEGER j
        REAL(dpk)  yy, y, drmax, fh, aug
        !EXE!
        
        hh     = rmax/DBLE(no-1)
        y      = -LOG( 2.0_dpk * ASIN(rmin/rmax) / m_pi) / LOG(DBLE(no-1) )
        drmax  = 2.0_dpk * rmax**y / ( m_pi * y )
        yy     = y - 1.0_dpk
        print*, 'y = ', y

        r(1)   = 0.0_dpk
        dr(1)  = 0.0_dpk
        DO  j  = 2, no
           fh    = DBLE(j-1) * hh
           aug   = ( m_pi / 2.0_dpk ) * ( fh / rmax)**y
           r(j)  = rmax  *  SIN(aug)
           dr(j) = drmax * TAN(aug) / fh**yy
        END DO
        

      END SUBROUTINE rsin
      !S
      !S       linear knots (idbsp = 1)
      !S
      SUBROUTINE rlin
        !
        IMPLICIT NONE
        !LOC!
        INTEGER j
        !EXE!
        
        hh      = rmax / DBLE(no-1)
        DO  j   = 1, no
           r(j) = hh *  (j-1) 
        END DO
        dr     = hh
        
      END SUBROUTINE rlin
      !
      !
      !   idbsp = 3
      !
      !   r = ri*dexp( (x+xi)**(1/nrp))
      !  dr = dx  r/[nrp * (ln(r/ri))**(nrp-1)]
      !
      SUBROUTINE rexp
        !
        IMPLICIT NONE
        INTEGER j
        REAL(dpk)  xmax,rho, x, ri
        !
        

        rho = 1.0_dpk / DBLE(nrp)
        ri  = rmin

        IF(xi.GT.0.0_dpk)  THEN
           ri = r(1) / dexp((xi)**rho)
        ENDIF
        
    
        xmax = (LOG(rmax/ri))**nrp - xi
        
        hh    = xmax/DBLE(no-1)
        
        DO j = 1, no
           x    = (j-1) * hh
           r(j) = ri * dexp( (x+xi)**rho )           
        END DO

        dr = DBLE(nrp) * ( LOG(r/ri) )**( nrp - 1 )

        
      END SUBROUTINE rexp
      !
      !
      ! rho_i = alpha * r_i + beta * ln(Z*r_i)
      !
      SUBROUTINE rlog
        !
        IMPLICIT NONE
        !
        INTEGER    :: j
        REAL(dpk)  :: alpha, beta
        REAL(dpk)  :: ri, dri
        !

        dri = rmax/DBLE(no-1)
        alpha = (rmax-LOG(znuc*rmax))/(rmax-rmin*LOG(znuc*rmax)/LOG(znuc*rmin))
        beta = 1.0_dpk - alpha * rmin/LOG(znuc*rmin)
       
        r(1) = 0.0_dpk
        DO j = 2, no
           ri    = (j-1) * dri 
           r(j)  = alpha * ri + beta * LOG(znuc*ri)
           dr(j) = dri * ( alpha + beta/r(j) ) 
        END DO
        

      END SUBROUTINE rlog
      !
      !
      !
    END MODULE set_grid
    !M
    !M
    !M       model potentials
    !M
    !M
MODULE potentials
  !
  USE PRECISION, ONLY:dpk
  USE param
  !
  PUBLIC v1e
  !  RRIVATE v_h, v_morse, v_he, v_hn, v_model
  
CONTAINS
!
!
!
  SUBROUTINE COREWF(c, ns, p, f, r, dr, t, npo, n, k, ni, nf, h)
    !
    USE param
    !
    IMPLICIT NONE

    
    !    IMPLICIT DOUBLEPRECISION(A-H,O-Z)
    REAL(dpk), INTENT(in), DIMENSION(:,:) :: c
    INTEGER                               :: ns
    REAL(dpk), INTENT(out), DIMENSION(:,:):: p
    REAL(dpk), INTENT(out), DIMENSION(:)  :: f
    REAL(dpk), INTENT(in),  DIMENSION(:)  :: r
    REAL(dpk), INTENT(in),  DIMENSION(:)  :: dr
    REAL(dpk), INTENT(in),  DIMENSION(:)  :: t
    INTEGER                               :: npo
    INTEGER                               :: n
    INTEGER                               :: k
    INTEGER                               :: ni
    INTEGER                               :: nf
    real(dpk)                             :: h
    !local
    REAL(dpk)                             :: coe(nb),work(3*k)
    REAL(dpk)                             :: oth, sq2
    INTEGER                               :: in
    INTEGER                               :: i, j, npp, inv
    REAL(dpk)                             ::  dbvalu, rint





!      DIMENSION C(NS,*), P(NL,*)
!      DIMENSION F(*), R(*), DR(*), T(*)
!      DIMENSION COE(NS), WORK(3*K) 
!..............................................................
    

    coe    = 0.0_dpk
    coe(1) = 0.0_dpk
    coe(n) = 0.0_dpk
    

    in = n - 1
    DO  npp = ni, nf
       in = in - 1
       DO  i = 1, n - 2
          coe( i + 1 ) = c(i, in)
       ENDDO
       
       !#         D'BOOR ROUTINE
        
       p = 0.0_dpk       !       p(npp, 1 ) = 0.0_dpk

       DO  j = 2, npo
            inv = 1 
            p(npp, j) = DBVALU(t, coe, n, k, 0, r(j), inv, work)
         ENDDO

!..... ?
         IF(p(npp, 6).GT.0.0_dpk) GOTO 14

         DO  j = 1, npo
           p(npp,j) = - P(npp,j)
        ENDDO

14      CONTINUE

        f(1) = 0.0_dpk

        DO  j = 2, npo

           f(J) = p(npp,j) * p(npp,j) * r(j) / dr(j)
        ENDDO
        


        oth = rint( f, 1, npo, 14, h)

        sq2= DSQRT(OTH)

        IF(sq2.EQ.(0.0_dpk)) THEN 
           WRITE(*,*) '# corewf :  problem for sq2 = ', sq2
           stop           
        ENDIF

        p(npp,:) = p(npp,:)/sq2

        !        DO j = 1, npo           
        !           P(npp,j) = p(NPP,J) / SQ2
        !        ENDDO

     ENDDO



    END SUBROUTINE COREWF



!
!  generate the initial estimates for the orbitals (hydrogenic estimates)
!
!
  SUBROUTINE init_hydrogenic_orbitals(pr,r)
    !
    use param
    !
    IMPLICIT NONE
    !
    REAL(dpk), INTENT(out), DIMENSION(:,:) :: pr
    REAL(dpk), INTENT(in),  DIMENSION(:)   :: r
    !locals
    INTEGER                                :: nfile
    INTEGER                                :: i, j, j_init
    !EXE!

    j_init = 2
    IF(r(1).NE.(0.0_dpk)) j_init = 1
    
    pr = 0.0_dpk

    !s-orbitals

    hydrogenic_s_orbitals:DO  i = 1, no_s

       DO j = j_init, no
          pr(i, j)  = hwf( i, lk(i), zk(i), r(j), ek(i) )
       ENDDO
       
    ENDDO hydrogenic_s_orbitals
    
    !p-orbitals
    hydrogenic_p_orbitals:IF(no_s.NE.ntc) THEN
       
       DO  i = no_s + 1, nsp
          
          DO j = j_init, no
             pr(i, j)  = hwf( i, lk(i), zk(i), r(j), ek(i) )
          ENDDO
       ENDDO
       
    ENDIF hydrogenic_p_orbitals
    
    !d-orbital
    hydrogenic_d_orbitals:IF(nsp.NE.ntc) THEN
       
       DO  i = nsp + 1, ntc

          DO j = j_init, no
             pr(i, j)  = hwf( i, lk(i), zk(i), r(j), ek(i) )
          ENDDO
       ENDDO
       
    ENDIF hydrogenic_d_orbitals
    
    !
    ! store the hydrogenic orbitals
    !

    nfile = 31
    OPEN(nfile, file='out/hwf.out')
    DO j = 1, SIZE(pr,dim=2)
       WRITE(nfile,*) r(j), (pr(i,j), i=1,SIZE(pr,dim=1))
    ENDDO
    CLOSE(nfile)
    

  END SUBROUTINE init_hydrogenic_orbitals
  !
  !
  !
  !
  !
  REAL(dpk) FUNCTION HWF(N, L, Z, R, E )
      
    INTEGER,   INTENT(in)   :: N 
    INTEGER,   INTENT(in)   :: L
    REAL(dpk), INTENT(in)   :: Z
    REAL(dpk), INTENT(in)   :: R
    REAL(dpk), INTENT(out)  :: E
    !LOC
    REAL(dpk)   :: FN
    REAL(dpk)   :: RK, RM, RLL, FK 
    REAL(dpk)   :: FM, FLL, A, P ,X
    !EXE
    
    M   = N + L
    K   = N - L - 1
    LL  = 2*L + 1
    FN  = N
    RK  = K
    RM  = M
    RLL = LL
    FK  = 1.0_dpk
    FM  = 1.0_dpk
    FLL = 1.0_dpk
    P   = 1.0_dpk
    A   = 1.0_dpk

    !......................


    X = - 2.0_dpk * Z * R / FN

      

    IF (ABS(X).GT.160_dpk) THEN
       
       HWF = 0.0_dpk 
         
    ELSE
       

       
       DO  I = 1, M
          
          FM = FM * RM
            RM = RM - 1.0_dpk
         ENDDO
         
         RM = M

         DO  I = 1, LL            
            FLL = FLL * RLL
            RLL = RLL - 1.0_dpk
         ENDDO
         
!#           2L+1
!#          L_(N+L)(X)  LAGUERRE POLYNOMIAL
!#
!#                    K = N + L - 1
!#
!#    F( A,B,X) = 1 + ( A/B * 1!) * X + (A(A+1) /B(B+1)* 2!) * X^2 + ...
!#
!#                    A * (A + 1)*  ... * (A + N)     1
!#            C_N =  -----------------------------   ---
!#                    B * (B + 1)* .... * (B+N)     N!
!#
!#
!#     A = - (N  - L - 1), B = 2 * L + 2
!#
!#
!#

!#
!#         N-L-1 = 1,2,3 
!#

!      IF (K) 1,2,3

         DO  I = 1, K
           
            P   = 1.0_dpk + A / RK * P / RM * X

            FK  = FK * RK
            A   =  A + 1.0_dpk
            RK  = RK - 1.0_dpk
            RM  = RM - 1.0_dpk
         ENDDO

         E   = - Z**2 / FN**2
         
         hwf = SQRT( Z * FM/FK) * P * EXP(X/2.0_dpk)*(-X)**(L+1)/FLL/FN
         
      ENDIF

    END FUNCTION HWF

!! potential

  REAL(dpk) FUNCTION v1e(x, l)
    !
    IMPLICIT NONE
    !
    INTEGER,   INTENT(in):: l
    REAL(dpk), INTENT(in):: x
    !
    REAL(dpk)            :: v
                 
    v1e = 0.0_dpk
    SELECT CASE(potential) 
       case('h') 
          v1e = v_h(x,l)                     
       case('he') 
          v1e = v_he(x,l)                     
       case('hn') 
          v1e = v_hn(x,l)                     
       case('morse') 
          v1e = v_morse(x)                     
       case('model')                 !various atoms 
          v1e = v_model(x,l)                     
       case('qdot')                 !various atoms 
          v1e = v_qdot(x)                     
       END SELECT



!    IF(potential.EQ.'h') THEN 
!       v1e = v_h(x,l)
!    ELSE IF(potential.EQ.'morse') THEN                      
!       v1e = v_morse(x) 
!    ELSE IF(potential.EQ.'he') THEN 
!       v1e = v_he(x,l)       
!    ENDIF


  END FUNCTION v1e
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! coulomb (hydrogenic ~ -Z/r) 
  REAL(dpk) FUNCTION v_h(x, l)
    !
    USE param                  !, ONLY:znuc
    USE units, ONLY:m_pi
    !
    IMPLICIT NONE
    
    INTEGER,   INTENT(in):: l
    REAL(DPK), INTENT(in):: x
    REAL(DPK)            :: w_x, v_model
    !


    v_h =  -znuc/x - (1.0_dpk  -  EXP( - ( x/rc )**6 )) * (ap/x**4)


  END FUNCTION v_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
  REAL(dpk) FUNCTION v_morse(x)
    
    IMPLICIT NONE
    
    REAL(dpk), INTENT(in):: x
    REAL(dpk)            :: D, a, xe
    
    

!! morse potential for diatomic molecules.

    D  = 0.1026_dpk
    xe = 1.9972_dpk
    a  = 0.732_dpk
    
    v_morse =  D * ( 1.0_dpk  -  EXP( - a*( x - xe ) ) )**2 - D
                  

  END FUNCTION v_morse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#
!# Model potential. Form and data taken by erc thesis (pg 195)
!#                  Also see: PRA,41,3534,1990 by Bachau H, Galan P, Martin F.
!#        
!# V(r) = -z/r + (zc/r)*(1 - (1+ar)*exp(-2ar)
!#
!# ap(1) == zc
!# rc(1) == a
!
!#       #    znuc    #     zc      #   a 
!..............................................
!# Li    #     3      #      2      #   1.65
!# Be+   #     4      #      2      #   2.35 
!# B++   #     5      #      2      #   3.04 
!# Mg+   #    12      #     10      #   2.13 
!# He    #     2      #      1      #   1.688 
!..............................................

  REAL(dpk) FUNCTION v_he(x, l)
    
    IMPLICIT NONE
    
    INTEGER l
    REAL(dpk) x
    REAL(dpk) w_x, v_model

! Coulombing (-Z/r) + model potential (imodel = 1)  

!
    IF(l.EQ.0) THEN 
       v_model = (0.4D+00/x) * EXP(-0.26*x*x) + (0.8/x**4)*(1 - EXP(-7.5*x**6)) -8.3*EXP(-4.5*x**2)       
    ELSE IF (l.EQ.1) THEN
       v_model =  (2.74267/x)*EXP(-0.943D+00*x**2) - 21.08D+00 * EXP(-2.5D+00*x**2) 
    ELSE IF (l.EQ.2) THEN
       v_model =  (2.74267D+00/x)*EXP(-0.943D+00*x**2) - 7.0D+00 * EXP(-2.0D+00*x**2)                         
    ELSE 
       v_model =  (2.74267D+00/x)*EXP(-0.943D+00*x**2) 
    ENDIF


!    v1e_he = 2.0D+00/x - l*(l+1)/x**2 + 2 * v_model

    v_he = -2.0_dpk/x - 2 * v_model

  END FUNCTION v_he
  !
  !
  !
  REAL(dpk) FUNCTION v_model(x, l)
    !    
    IMPLICIT NONE
    !
    REAL(dpk) :: x
    INTEGER   :: l
    REAL(dpk) :: w_x
    !



! cormier
    w_x = - 1.0_dpk + ( 1.0_dpk + rc * x ) * EXP( - 2.0_dpk * rc * x )
    v_model = w_x * ap / x 

  END FUNCTION v_model

!#
!# Negative hydrogen model potential. 
!#        
!# V(r) = 1.1*exp(-r)
!#
!..............................................

  REAL(dpk) FUNCTION v_hn(x, l)
    
    IMPLICIT NONE
    
    INTEGER l
    REAL(dpk) x
    REAL(dpk) w_x, v_model
    
    v_model = 0

!                  WRITE(*,*)  "# v1e_he::   zc, a = ", ap(1),rc(1) 
!                  stop
!    v1e_hn = 2.0D+00 * 1.1*EXP(-x) / x - l*(l+1)/x**2 + v_model


    v_hn = -2.0_dpk * 1.1_dpk*EXP(-x) / x - v_model

    
  END FUNCTION v_hn

!#
!# Quantum dot model potential.
!#        
!#  pra, 78, 075316, 2008
!#
!# V(r) = v0 * exp(-beta*r^2)
!#
!#   here v0 = znu
!#   rmax = sqrt(ln2.0/beta)
!..............................................

  REAL(dpk) FUNCTION v_qdot(x)
    !
    !
    USE param, ONLY: ap
!    USE param, ONLY: rc
    !
    IMPLICIT NONE
    !
    REAL(dpk)  :: x
    !exe
 

    !    beta = 0.1         !znuc == beta = LOG(2.0_dpk)/rmax**2

    v_qdot = - ap * EXP( - znuc * x**2 ) 

  END FUNCTION v_qdot



END MODULE potentials
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !*******************************************************
  !
  !  input :
  !         n = number of collocation points (no. B-splines)
  !         k = order of spline > 2
  !         a, b = 1st and last points on interval
  !  output :
  !         t(1) = t(2) = ... t(k) = a
  !         t(n+1) = t(n+2) = ... t(n+k)
  !         t(i) = a + (i-k)!h   i = k+1,n
  !         with  h = (b-a)/(n-k+1)
  !
  !**********************************************************

MODULE one_e_matrix
  
  USE param,      ONLY: dpk, nb, kb, rmax, rs, idbsp, znuc, nrp, xi
  USE potentials, ONLY: v1e 
  USE units,      only:m_pi  
  IMPLICIT NONE
  PUBLIC
CONTAINS 

  !##########################################################
!fff
  REAL(dpk) FUNCTION f(x,r,p)

    USE PRECISION, ONLY:dpk
    IMPLICIT NONE
    REAL(dpk), INTENT(in):: x,r
    INTEGER p
    !
    !executable statements!
    !
    
    f = 0.0_dpk 
    
    IF(r == 0.0_dpk) THEN     ! atomic case

!     f = x**p

     f = xp(x,p)
  ELSE                      ! multipole expansions

     IF (x == r) THEN 
!        f = 1/r                
       f = xp(r, -1)
     ELSE IF (x < r) THEN 

        !f = x**p/ r**(p+1)     !
        f = xp(x, p) / xp(r, p+1)
     ELSE IF(x > r) THEN     

        !f = r**p / x**(p+1)    !
        f = xp(r,p) / xp(x,p+1)
     ENDIF
  
  ENDIF
  
CONTAINS

!!!xp = x**p       
  REAL(dpk) FUNCTION xp(x,p)
!
    USE PRECISION, ONLY:dpk
    IMPLICIT NONE
    REAL(dpk), INTENT(in):: x
    INTEGER p
    !
    IF(p == 0) THEN
       xp = 1.0_dpk
    ELSE
       xp = x**p
    ENDIF
    !
  END FUNCTION xp

END FUNCTION f

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE mkgrid(t)
    !
    IMPLICIT NONE
    !
    REAL(dpk), INTENT(out), DIMENSION(:) :: t
    !locals
    INTEGER                              :: nknotfile
    REAL(dpk)                            :: gamma
    REAL(dpk)                            :: dr
    REAL(dpk)                            :: r_i, ri
    REAL(dpk)                            :: xmax
!    REAL(dpk)                            :: xi
    INTEGER                              :: i
    !EXE!


!    xi = 0.0_dpk
    !

    nknotfile  = 31

    !
    dr = rmax/DBLE( nb - kb + 1 )
    !    dr = rmax/DBLE( nb - kb + 1 )

    t = 0.0_dpk

    select_grid: IF(idbsp.EQ.0) THEN     ! sine-like 


       gamma  = -dlog( 2.0_dpk * dasin( rs/rmax ) /m_pi ) /  dlog( DBLE( nb - kb + 1))
       
       DO i = kb + 1, nb
          
          r_i =  DBLE( i - kb ) * dr
          ri  =  0.5_dpk * m_pi * (r_i/rmax)**gamma
          
          t(i) =  rmax * SIN(ri)

       END DO
       WRITE(*,*) '# set_grid::    sine grid,   rs(k+1), gamma = ', t(kb+1), gamma

    ELSE IF(idbsp.EQ.1) THEN    !linear

       
       DO i = kb, nb-1
          t(i+1) = t(i) + dr
       END DO

       !       WRITE(*,*) '# set_grid::    linear grid, rs(k+1) = ', t(kb+1)


    ELSE IF(idbsp.EQ.2) THEN    !exp


        gamma = 1.0_dpk / DBLE(nrp)

        ri = rs / EXP(xi**gamma)
        
        xmax = ( LOG(rmax/ri) )**nrp  - xi

        
        dr = xmax / DBLE( nb - kb + 1 )
                

        DO i = kb + 1, nb
           r_i  =  DBLE( i - kb + 1) * dr
           t(i) = ri * EXP( ( r_i + xi )**gamma)           
        END DO

  

     ENDIF select_grid



     t(nb + 1: nb + kb) = rmax

     OPEN(nknotfile, file='dat/knots-bs.dat')
     DO i = 1, SIZE(t)
        WRITE(nknotfile,*) i, t(i)
     ENDDO
     CLOSE(nknotfile)
     
   END SUBROUTINE mkgrid
   
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$!=====================================================================
!!$!
!!$!     sets up the knots for splines in semi-exponential grid:
!!$!
!!$!-----------------------------------------------------------------------
!!$      SUBROUTINE mkgrid1(z)
!!$        !
!!$      IMPLICIT NONE
!!$      !
!!$      REAL(dpk), INTENT(out), DIMENSION(:) :: t
!!$      !locals
!!$      INTEGER                              :: nknotfile
!!$      INTEGER(4)                           :: i,nt
!!$      Real(dpk)                            :: z
!!$      !
!!$
!!$      rmax = rmax * znuc
!!$      hmax = hmax * znuc
!!$
!!$      h = rmax/(nb-1)
!!$
!!$!      nt = rmax / h + 1 + kb      ! t must have dimensions of nt
!!$      
!!$ 
!!$      ! ... determine 'ml', the number of intervals from 0 to 1
!!$      ! ... and make sure that h = 1/n as recomended
!!$ 
!!$
!!$      
!!$      ml = 1.0_dpk/h + 0.5_dpk;  
!!$      h  = 1.0_dpk/ml
!!$
!!$
!!$      t = 0.0_dpk    !init
!!$
!!$      
!!$      ns = kb
!!$      linear_origin:DO i = 1, ml
!!$         ns    = ns + 1; 
!!$         t(ns) = t(ns-1) + h
!!$      END DO linear_origin
!!$      
!!$      ! determine 'me', the number of intervals in "exponential" grid
!!$
!!$      IF(hmax.LT.h) hmax = h
!!$      
!!$      me = 0
!!$      exponential:DO
!!$
!!$         t( ns + 1 ) = t(ns) * (1.0_dpk + h )
!!$         
!!$         IF( t(ns+1) - t(ns).GT.hmax) EXIT
!!$
!!$         ns = ns + 1
!!$
!!$         IF(t(ns).GT.rmax) EXIT
!!$
!!$         me = me + 1
!!$
!!$      END DO exponential
!!$
!!$! ... rest of interval with step = hmax
!!$
!!$      IF( t(ns).LT.rmax ) THEN
!!$         DO
!!$            t(ns+1) = t(ns) + hmax
!!$            ns = ns + 1
!!$            IF(t(ns).GE.rmax) EXIT
!!$         END DO
!!$      END IF
!!$
!!$      IF(t(ns).GT.rmax) THEN
!!$         ns = ns - 1; 
!!$         me = me - 1
!!$         t(ns+1:ns+ks) = rmax
!!$         t(ns) = (t(ns+1)+t(ns-1))/2.0_dpk
!!$      END IF
!!$      
!!$      nv = ns - ks + 1
!!$      
!!$      ! ... scale to the R variable
!!$
!!$      t = t/z
!!$      hmax = hmax / znuc
!!$      rmax = rmax / zncu
!!$
!!$    END SUBROUTINE mkgrid1
!!$    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mkgrid_mixed( i_0, n, k, R, t)

    IMPLICIT NONE
    INTEGER i_0, n, k, n_s, i
    INTEGER nknotfile
    REAL(dpk)  R, r_0, alpha, beta
    REAL(dpk), INTENT(out), DIMENSION(:) :: t

!......................

    nknotfile  = 31
!.....................

    open(nknotfile, file='dat/knot.dat')

    
    n_s = n - k + 2    ! number of points in [0,R] (included the boundaries) 
                       ! n - k + 1 : number of segments 

    r_0 = R * ( i_0 - 1 ) / dble( 2 * n_s - i_0 - 1)

    alpha = r_0 / dble(i_0 - 1)**2
    beta  = ( R - r_0 ) / ( n_s - i_0 ) 


    t = 0.0_dpk
    
    do i = k + 1, n

       if(I.lt.(I_0 + K )) then
          t(i) =  alpha * dble( i - k )**2 
       else        
          t(i) =  r_0 + beta * ( i - k + 1 - i_0 )
       endif

       write(nknotfile, * ) t(i), i  , abs(t(i) - t(i-1))

    end do
    close(nknotfile)
    do i = n + 1, n + k
       t(i) = R
    end do

  end subroutine mkgrid_mixed
  !
  !

  !######################################################################
  ! 
!!!!!!               note the sign changed from the original routine
  !
  !
  !#######################################################################
  !  input:
  !         n  = number of collocation points (number of B-splines)
  !         k  = order of B-splines >2
  !         t  = collocation points (?)
  !         AD = direct Hartree-Fock potential
  !         AE = exchange Hartree-Fock potential
  !  output: 
  !         a  = single-electron Hamiltonian matrix in B-splines basis
  !         b  = overlapping matrix of the B-splines basis
  !#######################################################################
  !
  !      b (i,j) = Int[ Bi(x) Bj(x) dx ] 
  !
  !      b(j,i) = b(i,j)
  !      
  !      for i<=j  :
  !                       t(i+k)  
  !               b(i,j) =  Int [Bi(x) Bj(x) dx]  j <i+k
  !                         t(j)
  !               b(i,j) =           0            j>=i+k
  !
  !      therefore
  !                        l=i+k-1  t(l+1)
  !               b(i,j)  =  Sum   { Int [ Bi(x) Bj(x) dx ] }
  !                          l=j     t(l)
  !
  !       where t(l), l = 1,n+k = knot sequence
  !
  !     t(l+1)                             k
  !      Int [ f(x) dx ]  = [t(l+1)-t(l)] Sum [ w(m) f(x(m)) ]
  !     t(l)                              m=1
  !
  !      x(m) = [t(l+1)-t(l)] z(m) + t(l) 
  !      z(m) = gaussian k-point coordinates for [0,1]
  !      w(m) = gaussian k-point weights     for [0,1]
  !      
  !#######################################################################


!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
! H = h_0 - l_bloch
!

  SUBROUTINE setmat(n, lang, a, b, t)
    !
    USE param
    USE UTILS,      ONLY : print_mx
    !
    implicit none
    !arg!
    INTEGER                             :: n
    INTEGER                             :: lang
    REAL(dpk), DIMENSION(:,:)           :: a 
    REAL(dpk), DIMENSION(:,:)           :: b 
    REAL(dpk), DIMENSION(:), intent(in) :: t
    !loc!
    INTEGER                             :: i, j, ii,jj
    INTEGER                             :: jhi, low, lhi, l, m, nderiv
    REAL(dpk)                           ::  seg,seh, dl, x0, r, fm, hm, v_l
    REAL(dpk), DIMENSION(kx,kx)         :: db
    REAL(dpk), DIMENSION(kx)            :: xg,wg
    REAL(DPK)                           :: a_s             ! boundary kinetic terms
    INTEGER                             :: n_first, n_last
    INTEGER                             :: nout 
    !
    data nderiv/2/
    external gauss, bsplvd

    !EX/

    nout = 11
    open(nout, FILE ='out/setmat.log') 

    !
    !
    !a(ndim,dim), b(ndim,ndim)     ndim = nb - 2 + first_spline + last_spline

    n_first = 2 - first_spline
    n_last  = n - 1 + last_spline


    WRITE(*,*) ' n_first = ', n_first, first_spline
    WRITE(*,*) '  n_last = ', n_last, last_spline
    WRITE(*,*) '  n,ndim = ', n, ndim

    CALL gauss(kb, xg, wg)


    fm = 0.0_dpk
    hm = 0.0_dpk

    a = 0.0_dpk
    b = 0.0_dpk
    loop_over_bsplines: do  i = n_first, n_last     !  set up non zero loops  i = 2 ... n  ,j = i,i+k-1

       jhi = min0(i + kb - 1, n_last)            

       loop_over_non_vanishing_matrix_elements:DO  j = i, jhi          ! < B_i | O |B_{i+k} >

          !   B-spline non-zero at [t(i),  t(i+k)], t(1->k) = t0

          low = max0(kb, j)
          lhi = min0( i + kb - 1, n)

          loop_over_segments:DO  l = low, lhi           !sum over intermediate segments

             seg = 0.0_dpk
             seh = 0.0_dpk

             x0 = t(l)
             dl = t(l+1) - x0


             DO  m = 1, kb                      !sum over gaussian weights

                r = dl * xg(m) + x0
                CALL bsplvd(t, kb, r, l, db, nderiv)
                
                fm = db( i - l + kb, 1 ) * db( j - l + kb, 1 )           !B_i  * B_j   
                hm = db( i - l + kb, 2 ) * db( j - l + kb, 2 )           !B'_i * B'_j

                !
                !        B_ij   = < B_i | B_j >
                !      h_l(i,j) = (1/2) * { B'_iB'_j  + B_i [ l(l+1)/r^2 + 2 V(r) ] B_j }
                !
                !
                ! V(r) = -Z/r (hydrogenic potential)

                v_l = lang * (lang + 1)/r**2

                seg = seg + wg(m) * fm
                seh = seh + wg(m) * ( 0.5_dpk *( hm + v_l * fm ) / mass  + v1e(r, lang) * fm )

             ENDDO

             
             ii = i - 1 + first_spline 
             jj = j - 1 + first_spline

             A(ii, jj) = A(ii, jj) +  dl * seh 
             B(ii, jj) = B(ii, jj) +  dl * seg 

          ENDDO loop_over_segments
       ENDDO loop_over_non_vanishing_matrix_elements
    ENDDO loop_over_bsplines
    !
    !    CALL bsplvd(t, kb, rmax, n_last, db_rmax, nderiv)
    !
    symmetrize:DO i = 1, n_last - n_first + 1
       DO j = 1, i                   
          a(i,j) =  a(j,i)
          b(i,j) =  b(j,i)
       END DO
    END DO symmetrize
    
    
    !    DO J = kb, 1, -1            
    !    WRITE(nout,'(a2,1x,a50,i4,1X,2E17.5)') "#","b(R), b'(R) :",&
    !     & n_last+j-kb, db_rmax(j, 1), db_rmax(j, 2)
    !ENDDO


    write(nout,*) '#'
    WRITE(nout,*) 'B'
    write(nout,*) '#'
    write(nout,*) 'b(n,n ) = ', B(n_last-1, n_last-1)
    write(nout,*) 'b(n-1,n) = ', B(n_last-2, n_last-1)
    write(nout,*) 'b(n,n-1) = ', B(n_last-1, n_last-2)
    WRITE(nout,*) 'h0'
    write(nout,*) ''
    write(nout,*) 'h(n,n) = ', A(n_last-1, n_last-1)
    write(nout,*) 'h(n-1,n) = ', A(n_last-2, n_last-1)
    write(nout,*) 'h(n,n-1) = ', A(n_last-1, n_last-2)
    WRITE(nout,*) 'h(n,n)- b(n,n) = ', A(n_last-1, n_last-1) - B(n_last-1, n_last-1)
    WRITE(nout,*) 'h(n,n) - a_s = ', A(n_last-1, n_last-1) - A_S


! S(i,j) = - (1/2) B_i(R) B'_j(R) 
!
! S(i,j)    =   0          , for all i,j except:
! S(n, n)   = - A_S        ,  A_S = (kb-1)/[2*(R-t_n)]
! S(n, n-1) =   A_S  
!

    A_S = -0.5D+00 * ( KB - 1 ) / ( t(n + 1 ) - t(n) )  

!    A( N - 1, N - 1 ) = A( N - 1, N - 1) + A_S         !     diagonal element of boundary term
!    A( N - 2, N - 1 ) = A( N - 1, N - 2) - A_S         ! non-diagonal element of boundary term
!    A( N - 1, N - 2 ) = A( N - 1, N - 2) - A_S         ! non-diagonal element of boundary term

    WRITE(nout,*) '#'
    WRITE(nout,*) '#'
    WRITE(nout,*) '#  bloch operator:     h_S = ', A_S
    WRITE(nout,*) '#'
    WRITE(nout,*) '# h0 + L(R)'
    WRITE(nout,*) '#'
    write(nout,*) 'h(n,n) = ', A(n_last-1, n_last-1)
    write(nout,*) 'h(n-1,n) = ', A(n_last-2, n_last-1)
    write(nout,*) 'h(n,n-1) = ', A(n_last-1, n_last-2)
    WRITE(nout,*) 'h(n,n)- b(n,n) = ', A(n_last-1, n_last-1) - B(n_last-1, n_last-1)
    WRITE(nout,*) 'h(n,n) - a_s = ', A(n_last-1, n_last-1) - A_S


    WRITE(nout,*) '#'
    WRITE(nout,*) '#'
    CLOSE(NOUT)
! print out matrices h0,B

!    CALL print_mx(n-1, a, 'h0', 'f')
!    CALL print_mx(n-1, b, 'b', 'f')
 
  END SUBROUTINE SETMAT
  !
  !
  ! inherited from h2p.f90
  !
  !
  !
  SUBROUTINE bsp_integral_f(n, t, mx_b, f, p, n_der, i_system)
    !
    USE param
    IMPLICIT NONE
    !
    INTEGER                             :: n
    REAL(dpk), INTENT(in), DIMENSION(:) :: t
    REAL(dpk), DIMENSION(:,:)           :: mx_b
    REAL(dpk)                           :: f
    INTEGER                             :: p
    INTEGER                             :: n_der
    INTEGER                             :: i_system ! 0 = atom, 1=molecule

!locals!
    INTEGER                             :: n_first, n_last
    INTEGER                             :: i, j, ii, jj
    INTEGER                             :: jhi, low, lhi 
    INTEGER                             :: l, m, nderiv
    real(dpk)                           :: dl, x0, r
    real(dpk),            dimension(kx) :: xg, wg
    REAL(dpk),         DIMENSION(kx,kx) :: bsp        !, bsp_rmax
    REAL(dpk)                           :: dm 
    REAL(DPK)                           :: h_s        ! boundary kinetic terms
    REAL(dpk)                           :: r_m_half
    integer                             :: nout
    !    data nderiv/2/
    EXTERNAL gauss, bsplvd
!................................

    nout = 16

    open(nout, FILE ='out/mat_bsp.log') 

    !................................
    
    !
    ! i_system = 0 (atomic case)      ! f(x,0,p) == x**p
    !          = 1 (molecular case)
    

    r_m_half = i_system * r_m/2.0_dpk

    !
    !

    n_first = 2 - first_spline
    n_last  = n - 1 + last_spline



    CALL gauss(kb, xg, wg)        ! prepare for gaussian integration over B-splines

    mx_b = 0.0_dpk


    DO  i = n_first, n_last     !  set up non zero loops  i = 2 ... n  , j = i,i+k-1

       jhi = min0(i + kb - 1, n_last )

       do  j = i, jhi

          low = max0(kb, j)                 ! B-spline none zero at [t(i),  t(i+k)], 
          lhi = min0( i + kb - 1, n_last)        !       t(1->k) = t0

          do  l = low, lhi           ! sum over intermediate segments

             dm  = 0.0_dpk

             dl = t(l+1) - t(l)
             x0 = t(l)

             DO  m = 1, kb              ! sum over gaussian weights

                r = dl * xg(m) + x0

                CALL bsplvd(t, kb, r, l, bsp, n_der+1)

                dm = dm + wg(m) * bsp(i-l+kb, n_der+1) * bsp(j-l+kb, n_der+1) * f(r,r_m_half,p) 


             enddo

             ii = i - 1 + first_spline 
             jj = j - 1 + first_spline

             mx_b(ii, jj) = mx_b(ii, jj) +  dl * dm 

          enddo
       enddo
    enddo
    !

    !    PRINT*, ndim, n_last
    !    stop
    DO i = 1, n_last-n_first+1              !     store remainder matrix part
       DO j = 1, i                   
          mx_b( i, j ) =  mx_b( j, i )
       END DO
    END DO
    

    h_s = 0.5_dpk * ( kb - 1 ) / ( t(n + 1 ) - t(n) )     ! S(i,j) = - (1/2) B_i(R) B'_j(R) 


!    mx_b( N - 1, N - 1 ) = mx_b( N - 1, N - 1) - n_der*h_s   ! last diagonal element
!    mx_b( N - 1, N - 2 ) = mx_b( N - 1, N - 2) + n_der*h_s   ! unsymmetric element of kinetic term

!       CALL bsplvd(t, kb, rmax, n, bsp_rmax, nderiv)         ! exactly on r = rmax
!    do J = kb, 1, -1           
!       WRITE(nout,'(a2,1x,a50,i4,1X,2E15.8)') "#","modules::mat_bsp:     B(R), B'(R) :",&
!                                                   &n+j-kb, db_rmax(j, 1), db_rmax(j, 2)
!    enddo

    WRITE(nout,*) '# modules::mat_bsp:    B-SPLINES MATRIX        '
    write(nout,*) '# modules::mat_bsp:           mx_b(   n,   n ) = ', mx_b(N-1, N-1)
    write(nout,*) '# modules::mat_bsp:           mx_b( n-1,   n ) = ', mx_b(N-2, N-1)
    write(nout,*) '# modules::mat_bsp:           mx_b(   n, n-1 ) = ', mx_b(N-1, N-2)
    write(nout,*) '# modules::mat_bsp:                       mx_s = ', h_s

    CLOSE(NOUT)
    
  END SUBROUTINE bsp_integral_f

  !
  !
  !
  !
  !  calculates b-splines dipole matrix elements <B_i|d|B_j>  
  !  
  !
  !
  SUBROUTINE bsp_dipole( n, d, t)
    !
    USE param
    USE UTILS,      ONLY : print_mx
    !
    IMPLICIT NONE
    INTEGER                             :: n
    REAL(dpk),           DIMENSION(:,:) :: d
    REAL(dpk), INTENT(in), DIMENSION(:) :: t
    !L
    REAL(dpk),         DIMENSION(kx,kx) :: db
    REAL(dpk),            DIMENSION(kx) :: xg,wg
    INTEGER i, j, jhi, low, lhi, l, m,  nderiv
    REAL(dpk)                seg, dl, x0, r, fm 
    INTEGER n_first, n_last
    !
    data nderiv/2/
    external gauss, bsplvd

    !EX/

    !....................

    WRITE(*,'(a60)') 'subroutine bsp_dipole in.'

    n_last  = n 
!    n_last  = n - 1 + last_spline
    n_first = 2 - first_spline

    d = 0.0D+00

    call gauss(kb, xg, wg)

    DO  i = n_first, n_last     !  set up non zero loops  i = 2 ... n  ,j = i,i+k-1

       jhi = min0(i + kb - 1, n_last )

       DO  j = i, jhi

          !   B-spline non-zero at [t(i),  t(i+k)], t(1->k) = t0

          low = max0(kb, j)
          lhi = min0( i + kb - 1, n_last)

          DO  l = low, lhi           !sum over intermediate segments

             seg = 0.0D+00

             dl = t(l+1) - t(l)
             x0 = t(l)

             DO  m = 1, kb                      !sum over gaussian weights

                r = dl * xg(m) + x0
                CALL bsplvd(t, kb, r, l, db, nderiv)
                
                fm = db( i - l + kb, 1 ) * r * db( j - l + kb, 1 )     !B_i *r* B_j

                seg = seg + wg(m) * fm


             ENDDO
             d(i - 1, j - 1 ) = d(i - 1, j - 1 ) +  dl * seg 
          ENDDO
       ENDDO
    ENDDO

    DO i = 1, n_last - 1
       DO J = 1, I
          d(i,j) =  d(j,i)
       END DO
    END DO
    

! print out matrices h0,B

!    CALL print_mx(n-1, a, 'd', 'f')

    WRITE(*,'(a60)') 'subroutine bsp_dipole out.'
 
  END SUBROUTINE BSP_DIPOLE
  !
  !  calculates 
  !      
  !   <B_i| d^q/dx^q | B_j>,    <B_i| r^q | B_j> 
  !
  !
  SUBROUTINE bsp_integral( n, d, t, deriv, q)
    !MOD!
    USE param
    USE UTILS,      ONLY : print_mx
    !PAR!
    IMPLICIT NONE
    !ARG!
    INTEGER                             :: n
    REAL(dpk),           DIMENSION(:,:) :: d
    REAL(dpk), INTENT(in), DIMENSION(:) :: t
    INTEGER                             :: deriv
    INTEGER                             :: q
    !LOC!
    REAL(dpk),         DIMENSION(kx,kx) :: db
    REAL(dpk),            DIMENSION(kx) :: xg,wg
    INTEGER                             :: i, j, ii, jj
    INTEGER                             :: jhi, low, lhi, l, m   
    REAL(dpk)                           :: seg, dl, x0, r, fm 
    INTEGER                             :: n_first, n_last
    INTEGER        nderiv
    DATA nderiv/2/
    external gauss, bsplvd
    !EXE!

    WRITE(*,'(a60)') 'subroutine bsp_integral in.'


    call gauss(kb, xg, wg)

    n_first = 2 - first_spline
    n_last  = n - 1 + last_spline      ! == n (last_spline = 1), 

    fm = 0.0_dpk
    d  = 0.0_dpk    
    sum_over_b_splines: DO  i = n_first , n_last 


       jhi = min0(i + kb - 1, n_last) 

       !       jhi = min0(i + kb - 1, n) ! set up non zero loops  i = 2 ... n  ,j = i, i + k - 1

       sum_over_the_overlaping_part_of_grid:DO  j = i, jhi 

          !   B-spline non-zero at [t(i),  t(i+k)], t(1->k) = t0

          low = max0(kb, j)
          lhi = min0( i + kb - 1, n)

          sum_over_segments:DO  l = low, lhi        !sum over intermediate segments

             seg = 0.0_dpk

             dl = t(l+1) - t(l)
             x0 = t(l)

             sum_over_gaussian_weights:DO  m = 1, kb                      !sum over gaussian weights

                r = dl * xg(m) + x0


                IF(deriv==0) THEN

                   CALL bsplvd(t, kb, r, l, db, 1)

                   fm = db( i - l + kb, 1 ) * db( j - l + kb, 1 ) 

                ELSE IF(deriv==1) then 

                   CALL bsplvd(t, kb, r, l, db, 1)
                   CALL bsplvd(t, kb, r, l, db, 2)

                   fm = db( i - l + kb, 1 ) * db( j - l + kb, 2 ) 

                ELSE IF(deriv==2) THEN

                   CALL bsplvd(t, kb, r, l, db, 2)

                   fm = db( i - l + kb, 2 ) * db( j - l + kb, 2 ) 
                ENDIF

                fm = r**q * fm
                  
                seg = seg + wg(m) * fm

             ENDDO sum_over_gaussian_weights

             ii = i - 1 + first_spline
             jj = j - 1 + first_spline

             d(ii, jj) = d(ii, jj)  +  dl * seg 

          ENDDO sum_over_segments
       ENDDO sum_over_the_overlaping_part_of_grid
    ENDDO sum_over_b_splines

!    assign_lower_matrix: DO i = 1, N - 1 
!       DO J = i, n - 1 
!          d(j,i) =   (1-2*deriv) * d(i,j) 
!       END DO
!    END DO assign_lower_matrix           
!    d(n-1,n-1) = -d(n-1,n-1) + deriv  ==> d(n-1,n-1) = deriv/2. !(velocity gauge) 

! after all  <B_i|d/dr|B_i> = 0(1/2) , i.ne.n (i==n)
    
       
! print out matrices h0,B

!    CALL print_mx(n-1, a, 'd', 'f')

  WRITE(*,'(a60)') 'subroutine bsp_integral out.'

  END SUBROUTINE BSP_INTEGRAL

  !
  ! calculates the B-splines expansion coefficients of function F(x) = S_j c(j) *B_j(x)
  ! note that B_splines is a non-orthogonal basis
  !
  ! (1) First calculate  by GL integration on the B-splines grid 
  !
  !              cb(i) = < B_i|F(x)> 
  !
  ! (2) then since B_i(x) non-orthogonal basis:
  !
  !             cb(i) = S_i c_i < B_i |B_j > = S_i c_i * B_ij ==>
  !
  !             CB = B * C  (1) 
  !
  ! therefore a solution of the above system (1) is performed by LU factorization
  !
  !             ==>   C = B^-1 * CB
  !
  !
  SUBROUTINE bsp_overlap(fx, b, cb)
    !
    USE param
    !
    IMPLICIT NONE
    INTEGER                                 :: n
    REAL(dpk), EXTERNAL                     :: fx
    REAL(dpk), DIMENSION(:,:)               :: b
    REAL(dpk), DIMENSION(:),  INTENT(inout) :: cb    
    ! locals                                !(b-splines integration)
    REAL(dpk), ALLOCATABLE, DIMENSION(:)    :: t
    REAL(dpk)                               :: segment 
    REAL(dpk), ALLOCATABLE,DIMENSION(:)     :: xg,wg
    REAL(dpk), ALLOCATABLE,DIMENSION(:,:)   :: db
    INTEGER    low, lhi
    INTEGER    deriv
    REAL(dpk)  dl, x0, r
    ! locals                         !( b*x=cb linear system solution)
    REAL(DPK), ALLOCATABLE, DIMENSION(:,:)  :: a 
    INTEGER,   ALLOCATABLE, DIMENSION(:)    :: iw
    INTEGER dim_1, dim_2, k
    INTEGER info
    real(dpk) norm
    INTEGER i,j,l,m
    !
    external gauss, bsplvd

    !EX/
!....................


    
    ALLOCATE( t(1:nb + kb)) ;    CALL mkgrid(t)           ! make b-splines grid
    ALLOCATE(xg(kb),wg(kb)) ;    CALL gauss(kb, xg, wg)   ! gauss weights 
    ALLOCATE(db(kb,kb))                                   ! b-splines matrix


    deriv = 1                                   ! only B_i(x) are needed
    loop_bsplines: DO  i = 2, nb                ! first-Bspline set to zero
       
       cb(i-1) = 0.0_dpk               
       low = max0(kb, i)                        ! B_i(x) non-zero only for x:
       lhi = min0( i + kb - 1, nb)              ! t(i) < x < t(i+k) 
       
       loop_segments:DO  l = low, lhi           !sum over intermediate segments
                                                ! t_i, t_(i+1), ..., t_(i+k)
          segment = 0.0D+00
          
          dl = t(l+1) - t(l)   ;  x0 = t(l)
          
          loop_over_k:DO  m = 1, kb             !sum over gaussian weights
             
             r = dl * xg(m) + x0
             CALL bsplvd(t, kb, r, l, db, deriv)
             segment = segment + wg(m) * db( i - l + kb, 1 ) * fx(r) 
          ENDDO loop_over_k
          
          cb(i-1) = cb(i - 1 ) +  dl * segment
       ENDDO loop_segments
    ENDDO loop_bsplines
       
    
    DEALLOCATE(t,xg,wg,db)

    k = kb - 1
    dim_1 = 3*kb - 2
    dim_2 = SIZE(b,1)
       
    ALLOCATE( a(dim_1, dim_2) )
           
    transform_full_symmetric_to_banded: DO j = 1, dim_2
       DO i = MAX(1, j - k), MIN( dim_2, j + k)

          IF(i>j) THEN             !because only the upper symmetric part is provided.
             a( 2*k + i - j + 1, j) = b(j,i)         
          ELSE
             a( 2*k + i - j + 1, j) = b(i,j)
          ENDIF
          
       ENDDO
    ENDDO transform_full_symmetric_to_banded
    
       !
    ALLOCATE(  iw(dim_2   ) )      

    iw   = 0   
    info = 0
    
    CALL F07BDF(dim_2, dim_2, k, k, a, dim_1, iw, info) ! LU factorization
    
    IF(info.NE.0) THEN
       WRITE(*,*) '# bsp_overlap: error in factorization subroutie  f07bdf.'
       STOP
    ENDIF
        
    CALL F07BEF('N', dim_2, k, k, 1, a, dim_1, iw, cb, dim_2, info)
    
    IF(info.NE.0) THEN
       WRITE(*,*) '# bsp_overlap: error in inversion subroutie  f07bef.'
       STOP
    ENDIF
    
    DEALLOCATE(iw)
    
  END SUBROUTINE BSP_OVERLAP
  
     
END MODULE one_e_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#######################################################################



