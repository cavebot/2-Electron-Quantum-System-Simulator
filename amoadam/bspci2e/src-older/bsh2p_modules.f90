!
! modules: set_grid, one_e_matrix, wf
!
!       set_grid: 
!             subroutines: rin, rsin, rexp, rexplin, cal_rpw
!
!   one_e_matrix:
!             subroutines: mkgrid, mkrgidmixed, setdipole,setmat
!    
!
!     potentials:
!            subroutines: 
!            functions  : hydrogenic, morse, woods_saxon, 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE set_grid

  USE param 

  IMPLICIT NONE
  
  PUBLIC

  REAL(dpk) hh
  REAL(dpk), DIMENSION(np) :: r, dr

  INTERFACE r_grid
     MODULE PROCEDURE rin
  END INTERFACE

  PRIVATE rin, rsin, rlin, rexp, rexplin
CONTAINS
!##################################################

  SUBROUTINE rin
    IMPLICIT NONE
    INTEGER NG1

    NG1 = 1

    SELECT CASE(idbsp)
    CASE(0)
       CALL rsin
    CASE(1)
       CALL rlin
    CASE(2)
       CALL rexp
    CASE(3)
       CALL rexplin
    END SELECT

    
    IF(code=='bs_fx') THEN
       OPEN(NG1, file='inp/grid.inp')
       WRITE(NG1,    '(I3)')   idbsp
       WRITE(NG1,'(3E10.2)')   rmin, rs, rmax
       WRITE(NG1,'(I6,1X,I3)') no, nrp
       CLOSE(NG1)
    ENDIF

    IF(problem=='2e') THEN
       OPEN(NG1, file='INP/GRID.INP')
       WRITE(NG1,'(I3)') idbsp
       WRITE(NG1,'(E10.2,1X,I3)' ) rmin, nrp
       WRITE(NG1,'(2E10.2,1X,I6)') rs, rmax, no
       CLOSE(NG1)
    ENDIF
    

  END SUBROUTINE rin

  !#################################################
!!! r = rmax sin[(pi/2) (x/rmax)**y]
!!! dr = rdx/dr
!!! dr = [2 rmax/(pi*y)]*[rmax/x]**(y-1) tan[(pi/2) (x/rmax)**y]
!!!

  SUBROUTINE rsin

    IMPLICIT NONE

    INTEGER j
    REAL(dpk)  pi, yy, y, drmax, fh, aug
!.......................


    hh      = rmax/DBLE(no-1)
    pi     = 3.14159265358979324D+00
    y      = -LOG( 2.0d+00 * ASIN(rmin/rmax) / pi) / LOG(DBLE(no-1) )
    drmax  = 2.0D+00 * rmax**y / ( pi * y )
    yy     = y - 1.0D+00
    r(1)   = 0.0D+00
    dr(1)  = 0.0D+00
    
    DO  j  = 2, no
       
!!!!             fh    = dfloat(j-1) * hh
       fh    = DBLE(j-1) * hh
       aug   = ( PI / 2.0D+00) * ( fh / rmax)**y
       r(j)  = rmax  *  SIN(aug)

!       r(j) =  hh * dble( j - 1 )
 
       dr(j) = drmax * TAN(aug) / fh**yy
       
    END DO
  END SUBROUTINE rsin
  !
  ! idbsp = 1, linear sequence
  !
  SUBROUTINE rlin

    IMPLICIT NONE
    INTEGER j
    REAL(dpk)  pi, yy, y, drmax, fh, aug
!.......................


    hh      = rmax/DBLE(no-1)
    pi     = 3.14159265358979324D+00
    y      = -LOG( 2.0d+00 * ASIN(rmin/rmax) / pi) / LOG(DBLE(no-1) )
    drmax  = 2.0D+00 * rmax**y / ( pi * y )
    yy     = y - 1.0D+00
    r(1)   = 0.0D+00
    dr(1)  = 0.0D+00
    
    DO  j  = 2, no
       
!!!!             fh    = dfloat(j-1) * hh
       fh    = DBLE(j-1) * hh
       aug   = ( PI / 2.0D+00) * ( fh / rmax)**y
       r(j)  = rmax  *  SIN(aug)

!       r(j) =  hh * dble( j - 1 )
 
       dr(j) = drmax * TAN(aug) / fh**yy
       
    END DO
  END SUBROUTINE rlin
  !#######################################################################

  SUBROUTINE rexp
    implicit none
    integer j
    real(dpk)  xmax,rho, x, ri
    !  r = ri*dexp( (x+xi)**(1/nrp))
    ! dr = dx  r/[nrp * (ln(r/ri))**(nrp-1)]

    rho = 1.0D+00 / dble(nrp)
    ri  = rmin

    if(xi.ne.0.d0)  ri = r(1) / dexp((xi)**rho)

    xmax = (log(rmax/ri))**nrp - xi

    hh    = xmax/dble(no-1)

    do j = 1, no

       x    = (j-1) * hh
       r(j) = ri * dexp( (x+xi)**rho )

    end do
    !        do j=1,no
    !          dr(j)=1.0d0
    !        end do
    !      if(nrp.eq.1) return
    do j = 1, no

       dr(j) = dble(nrp) * (log(r(j) / ri))**( nrp - 1 )

    end do

  end subroutine rexp

  !######################################################################

  subroutine rexplin

    implicit none
    integer j
    real(dpk)  xmax,rho, x, ri

    !  r = ri*dexp( (x+xi)**(1/nrp))
    ! dr = dx  r/[nrp * (ln(r/ri))**(nrp-1)]

    rho = 1.0d0 / dble(nrp)
    ri  = rmin

    if(xi.ne.0.d0)  ri = r(1) / dexp((xi)**rho)

    xmax = (log(rmax/ri))**nrp - xi

    hh   = xmax/dble(no-1)

    do j = 1, no

       x    = (j-1) * hh
       r(j) = ri * dexp( (x+xi)**rho )

    end do
    !        do j=1,no
    !          dr(j)=1.0d0
    !        end do
    !      if(nrp.eq.1) return
    do j = 1, no

       dr(j) = dble(nrp) * (log(r(j) / ri))**( nrp - 1 )

    end do

  end subroutine rexplin

end module set_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  module potentials
module potentials
  
  use precision, only:dpk
  use param

  PUBLIC v1e
!  RRIVATE v_h, v_morse, v_he, v_hn

contains

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  function
!  REAL(dpk) FUNCTION f(x,p)
!    USE precision, only:dpk
!    IMPLICIT NONE
!    REAL(dpk), INTENT(in):: x
!    integer p
!                  
!    f = x**p
!  END FUNCTION f
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  REAL(dpk) FUNCTION v1e(x)

!                  USE potentials
    IMPLICIT NONE
    REAL(dpk), INTENT(in):: x
    
           
    v1e = 0.0_dpk

    IF(potential.EQ.'h') THEN 
       
       v1e = v_h(x)
       
    ELSE IF(potential.EQ.'morse') THEN 
       
       v1e = v_morse(x) 
    ENDIF


  END FUNCTION v1e
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! coulomb (hydrogenic ~ -Z/r) 
  REAL(dpk) FUNCTION v_h(x)

    use param                  !, ONLY:za
    use units, only:m_pi
  
    implicit none
    real(DPK), intent(in):: x
    real(DPK) w_x, v_model
    
!...................

    w_x =  1.0D+00  -  EXP( - ( x/rc )**6 )
    v_model = (ap/x**4 ) * w_x
    

    v_h =  -za/x - v_model

  END FUNCTION v_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
  REAL(dpk) FUNCTION v_morse(x)

    IMPLICIT NONE
    
    REAL(dpk), INTENT(in):: x
    REAL(dpk) D, a, xe

    

!! morse potential for diatomic molecules.

    D  = 0.1026D+00
    xe = 1.9972D+00
    a  = 0.732D+00
    
    v_morse =  D * ( 1.0D+00  -  exp( - a*( x - xe ) ) )**2 - D
                  

  end function v_morse
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
!#       #    za    #     zc      #   a 
!..............................................
!# Li    #     3      #      2      #   1.65
!# Be+   #     4      #      2      #   2.35 
!# B++   #     5      #      2      #   3.04 
!# Mg+   #    12      #     10      #   2.13 
!# He    #     2      #      1      #   1.688 
!..............................................

  real(dpk) function v_he(x, l)

    implicit none
                  
    integer l
    real(dpk) x
    real(dpk) v_model
!    REAL(dpk) w_x

! Coulombing (-Z/r) + model potential (imodel = 1)  

! cormier
!                  w_x = - 1.0D+00 + ( 1.0D+00 + rc(1) * x) * EXP( - 2.0D+00 * rc(1) * x )
!                  v_model = ( ap(1) / x ) * w_x
!
    if(l.eq.0) then 
       
       v_model = (0.4D+00/x) * exp(-0.26*x*x) + (0.8/x**4)*(1 - exp(-7.5*x**6)) -8.3*exp(-4.5*x**2)
                     
    else if (l.eq.1) then

       v_model =  (2.74267/x)*exp(-0.943D+00*x**2) - 21.08D+00 * exp(-2.5D+00*x**2) 

    else if (l.eq.2) then

       v_model =  (2.74267D+00/x)*exp(-0.943D+00*x**2) - 7.0D+00 * exp(-2.0D+00*x**2) 
                        
    else 

       v_model =  (2.74267D+00/x)*exp(-0.943D+00*x**2) 

    endif


!    v1e_he = 2.0D+00/x - l*(l+1)/x**2 + 2 * v_model

    v_he = -2.0D+00/x - 2 * v_model

  end function v_he

!#
!# Negative hydrogen model potential. 
!#        
!# V(r) = 1.1*exp(-r)
!#
!..............................................

  real(dpk) function v_hn(x, l)
    
    implicit none
                  
    integer l
    real(dpk) x
    real(dpk) v_model
    
    v_model = 0.0_dpk

!                  WRITE(*,*)  "# v1e_he::   zc, a = ", ap(1),rc(1) 
!                  stop
    v_hn = 2.0D+00 * 1.1*EXP(-x) / x - l*(l+1)/x**2 + v_model


!    v_hn = -2.0D+00 * 1.1*exp(-x) / x - v_model


  end function v_hn


end module potentials
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

module one_e_matrix


  USE param,      ONLY: dpk, nb, kb, rmax, rs, idbsp, znuc, nrp, xi
  USE potentials, ONLY: v1e 
  USE units,      only:m_pi  

  implicit none
  public

contains 

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
!!!##############################################################

!  SUBROUTINE mkgrid_h2p( n, k, rmax, dr0, t)
  SUBROUTINE mkgrid_h1e(t)
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
    INTEGER                              :: i
    !EXE!


!    xi = 0.0_dpk
    !

    nknotfile  = 31

    !
    dr = rmax/DBLE( nb - kb + 1 )

    t = 0.0_dpk

    select_grid: IF(idbsp.EQ.0) THEN     ! sine-like 


       gamma  = -dlog( 2.0_dpk * dasin( rs/rmax ) /m_pi ) /  dlog( DBLE( nb - kb + 1))
       
       DO i = kb + 1, nb
          
          r_i =  DBLE( i - kb ) * dr
          ri  =  0.5_dpk * m_pi * (r_i/rmax)**gamma
          
          t(i) =  rmax * SIN(ri)

       END DO
       !       WRITE(*,*) '# set_grid::    sine grid,   rs(k+1) = ', t(kb+1)

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

 
     OPEN(nknotfile, file='dat/knot-bs.dat')
     DO i = 1, SIZE(t)
        WRITE(nknotfile,*) t(i), i
     ENDDO
     CLOSE(nknotfile)
     
   END SUBROUTINE mkgrid_h1e

  !##########################################################
  SUBROUTINE mkgrid_h2p( n, k, rmax, dr0, t)

    implicit none
    integer i, n, k
    integer nknotfile
    real(dpk) pi, y, dr0, hh, fh, aug, rmax
    real(dpk), intent(out), dimension(:) :: t

!......................

    nknotfile  = 31

!.....................

    open(nknotfile, file='dat/knot.dat')


    do i = 1, k

       t(i) = 0.0_dpk

    end do

    hh = RMAX/dble( n - k + 1 )

    PI = 3.14159265358979324_dpk
    y  = -dlog(2.0_dpk * dasin(dr0/rmax)/pi)/dlog(dble(n-k+1))

!!!  write(*,1000) ' t0 =',dr0,'  h =',hh


    !         write(*,*) 'Sine-knots is used for B-slpines'
    do i = k + 1, n

       fh   =  dble(i-k)*hh
       aug  =  (pi/2.0_dpk)*(fh/rmax)**y

       t(i) =  rmax*sin(aug)

!       t(i) =  hh * dble( i - k )

       write(nknotfile,*) t(i), i, abs(t(i) - t(i-1))

    end do

    close(nknotfile)

    do i = n + 1, n + k

       t(i) = RMAX

    end do

  END SUBROUTINE mkgrid_h2p


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mkgrid_mixed( i_0, n, k, R, t)

    implicit none
    integer i_0, n, k, n_s, i
    integer NKNOTFILE
    real(dpk)  R, r_0, alpha, beta
    real(dpk), intent(out), dimension(:) :: t

!......................


    nknotfile  = 31

!.....................

    open(nknotfile, file='dat/knot.dat')

    
    n_s = n - k + 2    ! number of points in [0,R] (included the boundaries) 
                       ! n - k + 1 : number of segments 

    r_0 = R * ( i_0 - 1 ) / dble( 2 * n_s - i_0 - 1)

    alpha = r_0 / dble(i_0 - 1)**2
    beta  = ( R - r_0 ) / ( n_s - i_0 ) 


    do i = 1, k

       t(i) = 0.0D+00

    end do

    
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

  SUBROUTINE mat_bsp(n, t, mx_b, f, p, n_der, i_system)
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
    INTEGER i, j
    INTEGER jhi, low, lhi 
    INTEGER l, m, nderiv
    real(dpk)  dl, x0, r
    real(dpk), dimension(kx) :: xg, wg
    REAL(dpk), DIMENSION(kx,kx) :: bsp   !, bsp_rmax
    REAL(dpk)  dm !, b_ij, bp_ij
    REAL(DPK) h_s                        ! boundary kinetic terms
    REAL(dpk) r_m_half
    integer nout
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

    CALL gauss(kb, xg, wg)        ! prepare for gaussian integration over B-splines

    mx_b = 0.0_dpk


    DO  i = 2, ndim + 1     !  set up non zero loops  i = 2 ... n  , j = i,i+k-1

       jhi = min0(i + kb - 1, ndim + 1 )

       do  j = i, jhi

          low = max0(kb, j)                 ! B-spline none zero at [t(i),  t(i+k)], 
          lhi = min0( i + kb - 1, n)        !       t(1->k) = t0

          do  l = low, lhi           ! sum over intermediate segments

             dm  = 0.0_dpk

             dl = t(l+1) - t(l)
             x0 = t(l)

             DO  m = 1, kb              ! sum over gaussian weights

                r = dl * xg(m) + x0

                CALL bsplvd(t, kb, r, l, bsp, n_der+1)

                   ! IF(n_der == 0 ) THEN                   
                   ! b_ij = bsp( i - l + kb, 1 ) * bsp( j - l + kb, 1 )     !B_i  * B_j                      
                   ! ELSE IF(n_der == 1) THEN
                   ! b_ij = bsp( i - l + kb, 2 ) * bsp( j - l + kb, 2 )     !B'_i * B'_j 
                   ! ENDIF
                   ! b_ij = bsp( i - l + kb, n_der+1 ) * bsp( j - l + kb, n_der+1)     !B_i  * B_j   
                   ! dm = dm + wg(m) * b_ij * f(r,r_m_half,p) 


                dm = dm + wg(m) * bsp(i-l+kb, n_der+1) * bsp(j-l+kb, n_der+1) * f(r,r_m_half,p) 


             enddo

             mx_b(i - 1, j - 1 ) = mx_b(i - 1, j - 1 ) +  dl * dm 

          enddo
       enddo
    enddo
    !

    do i = 1, ndim              !     store remainder matrix part
       do J = 1, I                   
            mx_b( i, j ) =  mx_b( j, i )
       end do
    end do


    h_s = 0.5_dpk * ( kb - 1 ) / ( t(N + 1 ) - t(N) )     ! S(i,j) = - (1/2) B_i(R) B'_j(R) 


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
    
  END SUBROUTINE mat_bsp


end module one_e_matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module data

!  USE PRECISION, ONLY:dpk


  PUBLIC write_en, write_v_mx, read_v_mx, write_wf1e 

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! store matrix elements
subroutine WRITE_MX(L, MX, FILENAME)

  use precision,  only : DPK
  use UTILS,      only : datafile, check
  use PARAM,      only : nb, nbs, ncs


  implicit none
!.... Arguments
  integer L
  real(DPK), dimension(:,:) :: MX
  character(LEN=*)  FILENAME
  !.... locals
  integer FILE
  integer I, J
  integer dim_1, dim_2 
  integer er
  !.........................


  FILE = 56

!........................

  dim_1 = size(mx,1)
  dim_2 = size(mx,2)

!  call check(dim_1, nb-1,er)
!  call check(dim_2, nb-1,er+1)
   
  CALL DATAFILE(FILE, L, FILENAME)
     
  WRITE(FILE) dim_1
  DO I = 1, dim_1
     WRITE(FILE) (MX(i,j), J = I, dim_2 )
     IF(i==10) THEN
        WRITE(*,'(10f10.3)') (MX(i,j),j=i,dim_2)
     ENDIF
  ENDDO
  CLOSE(FILE)

  write(*, '(a2,1X,a50,i3)') "#","write_mx: matrix elements for l = ",L
  write(*, '(a2,1X,a50,i3)') "#","write_mx: nstates = ", dim_1

  return

end subroutine WRITE_MX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! store vectors
SUBROUTINE WRITE_V(L, V, FILENAME)

  use precision,  only : DPK
  use UTILS,      only : ASCIIFILE
  use PARAM,      only : nb, nbs, ncs


  implicit none
!.... Arguments
  integer L
  real(DPK), dimension(:) :: V
  character(LEN=*)  FILENAME
  !.... locals
  integer FILE
  integer I, J
  integer dim 
  integer er
  !.........................


  FILE = 56

  !........................

  dim = size(V)

  !

  CALL ASCIIFILE(FILE, L, FILENAME)
     
  WRITE(FILE,'(i5)') dim

  DO I = 1, dim     
     WRITE(FILE,'(i5,2x,E15.8)') i, V(i) 
  ENDDO
  !
  CLOSE(FILE)
  !
  WRITE(*, '(a2,1X,a50,i3)') "#","Energy data for l = ",L
  WRITE(*, '(a2,1X,a50,i3)') "#","nstates = ", dim
  !
  RETURN
  
END SUBROUTINE WRITE_V
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write and read energies and coefficients for each l = 0,1,2
SUBROUTINE WRITE_V_MX(L, V, MX, FILENAME)


  use precision,  only: DPK
  USE UTILS,      ONLY: DATAFILE, CHECK
  USE param,      ONLY: nb, method, ndim
  USE PARAM,      ONLY: nbs, ncs, nl


  implicit none
!.... Arguments
  INTEGER                   :: L
  REAL(DPK), DIMENSION(:)   :: V
  REAL(DPK), DIMENSION(:,:) :: MX
  CHARACTER(LEN=*)            FILENAME
!... locals
  CHARACTER( LEN = 30 )       DFILE
  CHARACTER( LEN = 6  )       SL, SNL

  !.... locals
  INTEGER FILE
  INTEGER I, J
  INTEGER dim, dim_1, dim_2 
  INTEGER er
  !.........................


  FILE = 56

!........................
  dim = size(v)
  dim_1 = size(mx,1)
  dim_2 = size(mx,2)

!  WRITE(*,*) 'dim0,1,2 = ', dim,dim_1,dim_2

  CALL check(dim_1, dim, er)
  CALL check(dim_2, nl*ndim,er+1)

  IF(method =='l') THEN 
     CALL check(dim_1, nbs + ncs, er+2)
  ELSE IF(method == 'd') THEN 
     CALL check(dim_1, nl*(nb-2), er+2)
  ENDIF


  IF(nl==1) THEN
     CALL DATAFILE(FILE, L, FILENAME)
  ELSE

     WRITE(SL,'(I6)') L
     WRITE(SNL,'(I6)') NL
     DFILE = FILENAME//TRIM(ADJUSTL(SL))//"-"//TRIM(ADJUSTL(SNL))//".dat"

     CALL DATAFILE(FILE, NL, DFILE) !the DFILE is fixed-length character now
                                    !while subroutine  DATAFILE() expects 
                                    !unspecified length variable as it is 
                                    ! FILENAME.
                                    !thus no any further addition would work
                                    !in the body of sub DATAFILE.
  ENDIF

!  write(FILE) nbs, ncs

  WRITE(FILE) dim
  
  DO I = 1, dim_1
        !     WRITE(NFILE, '(I6,2X,E20.14)')  I, EN(I)
     WRITE(FILE)  V(I)  
     WRITE(FILE)( MX(I,J), j = 1, dim_2 )

  ENDDO
  
  CLOSE(FILE)

  WRITE(*, '(a2,1X,a50,i3)') "#","Energies/coefficients data for l = ",L
  WRITE(*, '(a2,1X,a50,i3)') "#","nstates = ", dim_1

  RETURN

END SUBROUTINE WRITE_V_MX
!.......................
SUBROUTINE READ_V_MX(L, V, MX, FILENAME)

  USE PRECISION,  ONLY : DPK
  use UTILS,      only : DATAFILE
  use PARAM,      only : nb


  implicit none
!.... Arguments
  integer L
  real(DPK), dimension(:), pointer   :: V
  real(DPK), dimension(:,:), pointer :: MX
  character(LEN=*)  FILENAME
  !.... local
  integer FILE
  integer I, J
  integer ns
  integer dim_1, dim_2 
  !.........................


  FILE = 56

!........................
 
  CALL DATAFILE(FILE, L, FILENAME)

  !  read(FILE) n_b, n_c

  READ(FILE) ns
     
!  call check(n_b, nbs, er)
!  call check(n_c, ncs, er+1)  
!  ns = n_b + n_c
     
  ALLOCATE( V( ns ) )
  ALLOCATE( MX( ns, nb-1 ) )
  
  dim_1 = SIZE(mx,1)
  dim_2 = SIZE(mx,2)
  
  DO I = 1, dim_1

     READ(FILE)  V(I)  
     READ(FILE)( MX(I,J), j = 1, dim_2 )
     
  ENDDO
  
  CLOSE(FILE)
  
  WRITE(*, '(a2,1X,a50,i3)') "#","Energies/coefficients data for l = ",L
  WRITE(*, '(a2,1X,a50,i3)') "#","nstates = ", dim_1
  
  RETURN

END SUBROUTINE READ_V_MX
!.......................
SUBROUTINE WRITE_WF1E(L, EN, CE)

  USE param
  USE set_grid
  USE one_e_matrix
  USE utils, ONLY: datafile
!  USE ioroutines
  
!!!......................
  
  implicit none
! arguments
  integer l                              ! angular number
  real(DPK), dimension(:)   :: EN        ! energy spectrum
  real(DPK), dimension(:,:) :: CE        ! coefficients
  
!  integer ns
! locals
  integer i, j, jj
  INTEGER FILE, OUTFILE
  integer NPLOT
  real(DPK), dimension(NB)      :: BCOEFF
  real(DPK), dimension(NB+KX)   :: T          !  r_i, knot sequence
  real(DPK), dimension(NP)      :: P          !  P(r_i)
  real(DPK), dimension(NP)      :: DP         !  dP(r_i)/dr
  real(DPK)  A_S
  real(DPK)  BVALUE
  CHARACTER*100 ARGV

!...............................

  CALL GETARG(2, ARGV)
  READ(ARGV,*)   NPLOT         

!...............................

  FILE = 56
  OUTFILE = 46

!...........................
  
  CALL R_GRID
  CALL MKGRID(NB, KB, RMAX, RS, T )

  
  CALL DATAFILE(FILE, L,"hwf1e-")


  WRITE(FILE) SIZE(EN)
  WRITE(FILE) NO, HH
  WRITE(FILE) ( R(J), J = 1, NO )
  WRITE(FILE) (DR(J), J = 1, NO )

!!!...................
         

  DO I = 1, SIZE(EN)   ! loop over the states of the symmetry 'l'
           
     bcoeff(1) = 0.0D+00
     do JJ = 1, NB - 1 

        bcoeff( jj + 1) = ce(i, jj)     
     end do



     P(1) = 0.0D+00
     do J = 2, NO

        P(J) = BVALUE( T, BCOEFF, NB, KB, R(J), 0 )
       DP(J) = BVALUE( T, BCOEFF, NB, KB, R(J), 1 )

     end do

     A_S = 0.5D+00 * ( KB - 1 ) / ( t(NB + 1 ) - t(NB) )

     DP(NO) = 2 * A_S * ( BCOEFF(NB) - BCOEFF(NB-1) )

!.....................

!  Plot radial wavefunction    P_nl( r)  [0, R]  

     WRITE(FILE)  EN(I) 
     WRITE(FILE) ( P(J), J = 1, NO)
     WRITE(FILE) (DP(J), J = 1, NO)

!....

     IF(I.EQ.NPLOT) THEN

        OPEN(OUTFILE, FILE="out/hwf1e.out")

        WRITE(*, '(a2,1X,a50,2i5,a2)') "#"," Writting p,dp in out/wf1e.inp for (",i,l,")"
        WRITE(*, '(a2,1X,a50,e15.8)') "#","e(n,l) = ", en(i)


        WRITE(*, *) '#'

        DO J = 2, NO

      WRITE(OUTFILE, *)  R(J), P(J), DP(J) 
      
   ENDDO

   CLOSE(OUTFILE) 

ENDIF

  ENDDO

  CLOSE(FILE) 
 
!!  !..............................

END SUBROUTINE WRITE_WF1E


end module data

!#######################################################################



