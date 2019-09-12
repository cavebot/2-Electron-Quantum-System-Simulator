!
! modules for k2e.f90 
!
! modules: 
!      get_space: set_space, set_config
!            rwf: cal_wf, bvalue
!           norm : nom, delta, phi, zeta, pscoul



MODULE GET_SPACE
  !
  USE PRECISION,ONLY:DPK
  !
  IMPLICIT NONE
  !
  PUBLIC
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: nhf,  lhf, ll,  nmin, nmax
  INTEGER, ALLOCATABLE, DIMENSION(:)   ::  is, noll, ndi, idcs
  REAL(DPK), ALLOCATABLE, DIMENSION(:) :: coe

  INTERFACE SET_SPACE
     MODULE PROCEDURE set_config
  END INTERFACE
  
  PRIVATE set_config
  
CONTAINS
  !
  SUBROUTINE set_config(nd, nbsp)
    !
    IMPLICIT NONE
    !
    INTEGER  nd, nbsp
    !..
    ALLOCATE( ndi(nd)  )
    ALLOCATE( nhf(nd)  )
    ALLOCATE( lhf(nd)  )
    ALLOCATE( ll(nd)   )
    ALLOCATE( nmin(nd) )
    ALLOCATE( nmax(nd) )
    ALLOCATE( noll(nd) )
    ALLOCATE( is(nd)   ) 
    ALLOCATE( idcs(nd) )      
    ALLOCATE( coe(nbsp) )
    !..
    WRITE(*,*) '# k2e:: set_config: seting dimensions for configurations done'
    WRITE(*,*) '# k2e:: set_config:                  nof cfg series    ncs = ', nd
    !
  END SUBROUTINE set_config
END MODULE GET_SPACE

!!####################################################################
MODULE rwf
  !
  USE PRECISION,ONLY: dpk
  !
  USE param
  USE set_grid
  USE one_e_matrix
  !
  IMPLICIT NONE
  !
  PUBLIC
  REAL(dpk), ALLOCATABLE, DIMENSION(:) :: wf
  REAL(dpk),          DIMENSION(ns+nk) :: t
  !
CONTAINS

!!!#############################

  SUBROUTINE cal_wf(lang, coe, np)
    !mod!
    USE PRECISION,ONLY: DPK
    USE bs_frb_2e, only: bvalue
    !
    IMPLICIT NONE
    !arg!
    INTEGER                  :: lang
    REAL(dpk), DIMENSION(ns) :: coe
    INTEGER                  :: np
    !loc!
    INTEGER                  :: j
 !   real(dpk)                :: bvalue
    !exe!
    

    CALL mkgrid( n2(lang+1), k, rmax, rs(lang+1), t)
    
    DO j = 1, np
       wf(j) = bvalue(t, coe, n2(lang+1), k, r(j), 0)
    END DO

  END SUBROUTINE cal_wf

END MODULE rwf

!!!#############################

MODULE norm
  !
  USE PRECISION,    only:dpk
  !
  USE PARAM,ONLY : ns
  !
  IMPLICIT NONE
  REAL(dpk) an , w1, pi, ph2, aa2
  !
  PUBLIC
  !
CONTAINS

  !----------------------------------------
  !   input:
  !        e  = energy 
  !        e1 = threshold energy
  !        l  = angular momentum q. number
  !        no = number of points for w.f.
  !        r() = r-grid points
  !        p() = wavefunction on grid
  !        z1  = nuclear charge
  !-----------------------------------------  

  SUBROUTINE nom(e, l, no, r, p, e1, z1, alpha, vp)
    !mod!
    USE PRECISION,ONLY : dpk
    !
    IMPLICIT NONE
    !
    REAL(dpk), PARAMETER :: EnAU = 27.211396181D+00
    !arg!
    real(dpk)                :: e
    integer                  :: l
    integer                  :: no
    REAL(dpk), DIMENSION(no) :: r
    REAL(dpk), DIMENSION(no) :: p
    real(dpk)                :: e1
    real(dpk)                :: z1
    real(dpk)                :: alpha
    real(dpk)                :: vp
    !
    REAL(dpk), DIMENSION(ns) :: del
    REAL(dpk)                :: ek, g1, aq, a1
    INTEGER                  :: id1
    INTEGER                  :: i,j
    !exe!

                
    ek = dsqrt( e - e1 )               !  k  =  sqrt( E - Eth )

    pi = dasin(1.0D+00) * 2.0D+00

    !write(*,*), ' E = ', e, ' Ethr = ', e1, ' Ee = ', ek, ' Zeff', z1, ' l = ', l
    !      write(*,*), ' k  = ', ek,  ' no = ', no
    !      do i=no-10, no
    !        write(*,*) r(i), p(i)
    !      end do
    
    !  
    !     Fitting  to the last 10 points 
    !

    id1 = 0
    DO  j = no - 10, no
       
       id1 = id1 + 1

       del(id1) = delta( r(j-1), r(j), p(j-1), p(j), l, ek, z1, alpha, vp)

       aa2 = zeta( r(j), ek, l, z1, vp)
                   
    ENDDO


    WRITE(*,301) r(no), del(id1), aa2

    


!!!............................

!
!      w1   = 0.d0
!      do i = 1, id1
!        w1 = w1 + del(id1)
!        w1 = w1 + del(i)
!      end do 
!      w1 = w1 / id1
!

    ph2 = phi( r(no), ek, z1, l, alpha, vp)
    aa2 = zeta( r(no),ek, l, z1, vp)
    w1  = delta( r(no-1), r(no), p(no-1), p(no), l, ek, z1, alpha, vp)

    g1  = dsin(ph2 + w1)




!C*
!C*
!C*         Evaluate the amplitude A1 of the discretised continuum
!C*         wavefunction by dividing the value of that discretised
!C*         wavefunction (at the {no-50}th grid point) by 
!C*                 sqrt{k/\zeta} * sin(\phi+\delta),
!C*         assuming that the continuum wavefunction \xi fulfiles 
!C*                 \xi = A1 * sqrt{k/\zeta} * sin(\phi+\delta) 
!C*         at large values of r:
!C*
!C*          A1 = p(no-50) / g1 / aq
!C*
!C*         Calculate the renormalisation factor from the amplitude A1,
!C*         i.e. calculate the prefactor required to change the 
!C*         normalisation of the continuum wavefunction from A1^2 
!C*         ("wrong" normalisation of the discretised wavefunction)
!C*         to the "correct" nomalisation 2/(\pi*k):
!C*         (In order to renormalise the wavefunction, it has to be 
!C*          multiplied by sqrt{correct factor/wrong factor}, and thus
!C*          by
!C* 
!C*            an = 1 / AN,   AN  = sqrt{ 2 /(\pi*k) / A1}) 
!C*

    an  = p(no) / dsqrt(2.0D+00 / pi/ aa2 )/g1 


    WRITE(*, 300) ek**2 * EnAU, w1, an**2

    WRITE(*,*) 'xxxxxxxxxxxxx'

!      output phase and zeta
!
!      ph2 = phi(r(no),ek,z1,0.d0,l ) + w1
!      ph2 = phi(r(no),ek,z1,0.d0,l)
!      aa2 = zeta(r(no),ek,z1,0.d0,l)


300 FORMAT(2x,' A(eV) = ',e15.7, ',   w = ',e15.7,', |An|^2  = ',e15.7)
301 FORMAT(2x,' R(no) = ',G15.3, ', del = ',e15.7,',   zeta  = ',e15.7)

  CONTAINS

!!!......................

    REAL(dpk) FUNCTION delta(r1, r2, p1, p2, l, ek, znuc, alpha, vp)
      !
      USE PRECISION, only: dpk
      !
      IMPLICIT NONE
      !arg!
      REAL(dpk)         :: r1
      REAL(dpk)         :: r2
      REAL(dpk)         :: p1
      REAL(dpk)         :: p2
      REAL(dpk)         :: ek
      REAL(dpk)         :: znuc
      REAL(dpk)         :: alpha
      REAL(dpk)         :: vp
      !loc!
      REAL(dpk)         :: ph1, ph2 
      REAL(dpk)         :: aa1, a1
      REAL(dpk)         :: aa2, a2
      !
      INTEGER l
      !
      ph1   = phi( r1, ek, znuc, l, alpha, vp)
      aa1   = zeta( r1, ek, l, znuc, vp)
      a1    = 1.0D+00 / dsqrt(aa1)
      
      ph2   = phi( r2, ek, znuc, l, alpha, vp)
      aa2   = zeta(r2, ek, l, znuc, vp)
      a2    = 1.0D+00 / dsqrt(aa2) 


      delta = datan((a1 * dsin(ph1) * p2 - a2 * dsin(ph2) * p1) / &
           &       ( a2 * dcos(ph2) * p1 - a1 * dcos(ph1) * p2))


    END FUNCTION delta

!!!   x    -- r,  
!!!   q    -- momentum(k),  
!!!   z    -- charge at large distance,
!!!  alpha -- polarization potential,  
!!!   l -- angular momentum,
!!!****   zeta approaches q when x becomes very large
!!!################################################################

    REAL(dpk) FUNCTION zeta(x, q , l, z, vp)
      !
      USE PRECISION, only: dpk
      !
      IMPLICIT NONE
      !
      INTEGER l
      REAL(dpk) x,q,z,vp,c,zd,omega,omega1,omega2
      !
!!................

      
      c = DBLE( l*(l+1) ) / x**2

      zd =  z / x

!!$      WRITE(*,*) ' l    = ', l
!!$      WRITE(*,*) ' q^2  = ', q**2
!!$      WRITE(*,*) ' c    = ', c * x**2
!!$      WRITE(*,*) ' z    = ', z
!!$      WRITE(*,*) ' vp   = ', vp

         omega   =   q**2 - c + 2.0D+00 * zd - vp

         omega1  =  (4.0D+00/x**2) * ( c + 2.0D+00 * vp  - zd)**2

         omega2  = -(2.0D+00/x**2) * (3.0D+00 * c + 10.0D+00 * vp - 2.0D+00 *zd)

      
         zeta    = dsqrt(omega) + 0.15625D+00 * omega1 * omega**(-2.5)&
              &    - 0.125D+00 * omega2 * omega**(-1.5)


    END FUNCTION zeta

!!!................................................


!!!   x  --  r  
!!!   q  -- momentum(k)  
!!!   z  -- charge at large distance
!!!   vp -- polarization potential  
!!!   l  -- angular momentum

    REAL(dpk) FUNCTION phi( x, q, z, l, a, vp)
      !
      USE PRECISION, only:dpk
      !
      IMPLICIT NONE
      !
      INTEGER l
      REAL(dpk) a, dln, invcos 
      REAL(dpk) x, q, z, vp, halfpi, c, qr, xi 
      REAL(dpk) cc, fk, rho, qqc, theta, vptheta
      REAL(dpk) phi_r1, phi_r2, phi_r4
      !



      qr = q * x

      halfpi = 1.5707963268D+00

      c  = float( l * ( l + 1 )) + a 


!      WRITE(*,*) ' l    = ', l
!      WRITE(*,*) ' a    = ', a
!      WRITE(*,*) ' c    = ', c
!      WRITE(*,*) ' z    = ', z
!      WRITE(*,*) ' q-c  = ', qr**2 - c 

      IF(z.EQ.0.D+00) THEN 


!!!
!!!  potentials of type :   (Z_eff = 0.0D+00) 
!!!
!!!           V(r) = - [ l(l+1) + a ] / r^2  - vp / r^4  + O( r^{-x} )
!!!

         !! ****  NEGATIVE ION CASE  : H- , ...


         xi =  dsqrt( qr**2 - c )

         
         IF(xi.LT.0.0D+00) THEN

            WRITE(*,*) ' (k*r)^2 -c  NEGATIVE '
            WRITE(*,*) '     k = ', q**2
            WRITE(*,*) '     r = ', x**2
            WRITE(*,*) ' l,a,c = ', l, a, c

            STOP
            
         END IF

!!!xxxx   0.208333 ==  5/24
!!!xxxx   0.125    ==  1/8


!!!         SELECT CASE(l)


         IF(c.EQ.0.0D+00) THEN

            !!         CASE(0)

            theta   = 0.125D+00 / qr

            vptheta = vp * 0.5D+00 / ( q * q * x )
 

         ELSE IF(c.GT.0.0D+00) THEN

!!         CASE(1:)

            cc = dsqrt(c)

            theta   = ( ( c + 0.125D+00) / cc) * dasin( cc / qr ) 

            vptheta = vp / ( x * ( xi + q*x) )  

         ELSE

            cc = dsqrt( -c )
            
            dln = DLOG( ( cc  + xi ) / qr )

            theta   = ( ( c + 0.125D+00) / cc) * dln

            vptheta = vp / ( x * ( xi + q*x ) )  


         ENDIF

!!         END SELECT
!!!  2*Z/r             contribution 

         phi_r1 = 0.D+00


!!!   - [ l(l+1) + a ] / r^2          contribution 

         phi_r2 = xi - halfpi * float(l) - 0.125D+00 / xi    &
              & - 0.208333333D+00 * c/ xi**3               &
              & + theta

!!!  -vp/r^4         contribution 


         phi_r4 =  0.5D+00 * vptheta                   
           


         phi = phi_r1 + phi_r2 + phi_r4 


!!!...............  polarization contribution 
!!$         IF(vp.NE.0.0D+00) THEN
!!$
!!$         SELECT CASE(l)
!!$
!!$
!!$         CASE(0)
!!$
!!$            vptheta = vp * x / (6.0D+00 * q )
!!$
!!$         CASE(1:)
!!$
!!$            vptheta = 0.5D+00 * ( ( qr * x)**2 / ( 2.d0 * c ) * &
!!$                 &  vp/ cc * dasin(cc/qr) - xi * x**2 * vp /( 2.0D+00*c) )
!!$
!!$         END SELECT
!!$
!!$            phi = phi + vptheta 
!!$
!!$         ENDIF

!!!...............  


         
      ELSE 


!!!  potentials of type : 
!!!
!!!           V(r) = 2*Z/r - [ l(l+1) + a ] / r^2  - vp / r^4  + O( r^{-x} )
!!!

         !****       v(n-z) potential

         fk  = q / z
         rho = x * z
         xi  = dsqrt( qr**2 + 2.0D+00 * rho - c )
         qqc = fk * fk * c


         IF(xi.LT.0.0D+00) THEN

            WRITE(*,*) ' (k*r)^2 -c  NEGATIVE '
            WRITE(*,*) '     fk  = ', qr**2
            WRITE(*,*) '     rho = ', rho**2
            WRITE(*,*) ' l,a,c   = ', l,a,c

            STOP
            
         END IF


!!!         SELECT CASE(l)

         IF(c.EQ.0.0D+00) THEN

!!         CASE(0)

            theta   = 0.25D+00/ ( xi + qr )

            vptheta = vp *   2 * z * ( xi + 2 * fk * rho )/( 3* rho *(xi + fk * rho )**2 )


         ELSE IF(c.GT.0.0D+00) THEN
!!!         CASE(1:)

            cc = dsqrt( c )

            invcos = dacos(( rho - c + fk * c *xi ) / ( (1.0D+00 + qqc) * rho ) )


            theta = ( ( c + 0.125D+00) / cc ) * invcos 
            
            vptheta = vp * (  - z * ( 2 * rho - c )  / ( c * rho *(xi + fk * rho ))&

            & + (z/ c * cc) * invcos  )

         ELSE

            cc = dsqrt( - c )

            dln = DLOG( ( rho - c + xi * cc ) / ( (1.0D+00 + fk *cc ) * rho ) )

            theta = ( ( c + 0.125D+00) / cc ) * dln 


            vptheta = vp * (  - z * ( 2 * rho - c )  / ( c * rho *(xi + fk * rho ))&

            & + (z/ c * cc) * dln  )
            
         ENDIF



!!         END SELECT


!!!  2*Z/r             contribution 


         phi_r1 =  - ( 1.0D+00 + dlog(fk)) / fk + pscoul( l, z, q)


!!!   - [ l(l+1) + a ] / r^2          contribution 

         phi_r2 = xi                                               &  
              & + dlog( 1.0D+00 + qr * fk + fk * xi )/fk           &
              & -  halfpi * float(l)                               &                                
              & - ( xi * ( 3.0D+00 * qqc + 4.0D+00 )               &
              & + qr *( 3.0D+00 * fk*fk + 2.0D+00) + fk * c ) /    &    
              & ( 24.0D+00 * (1.0D+00 + qqc) * xi * (xi + qr ) )   &  
              & + 0.208333333D+00*( rho - c ) / xi**3              &
              & + theta

!!!  -vp/r^4         contribution 


              phi_r4 = 0.5D+00 * vptheta        


              phi = phi_r1 + phi_r2 + phi_r4 


!!!...... polarization contribution 
!!$         IF(vp.NE.0.0D+00) THEN
!!$
!!$            SELECT CASE(l)
!!$
!!$            CASE(0)
!!$
!!$               vptheta = (2.d0/15.d0) * vp * x**2 * &
!!$                       & (9.d0*qr*x + 11.d0 * qr**2 + 6.d0 *rho) /( xi + qr)**3
!!$
!!$            CASE(1:)
!!$
!!$               vptheta = z**2*vp*x**4/(2.d0*c**2)*(3.d0+qqc)*&
!!$                       & dacos((rho-c+fk*c*xi)/(1.d0+qqc)/rho)/cc&
!!$                       & -z*vp*x**3/(2.d0*c**2)*(3.d0*(2.d0*rho-c)/&
!!$                       & (xi+qr) + c * xi/rho)
!!$
!!$            END SELECT

           ENDIF

      
!           WRITE(*,*)  ' Z/r                    COULOMBIC PHASE = ', phi_r1
!           WRITE(*,*)  ' [l(l+1) + a ]/r^2      CI        PHASE = ', phi_r2
!            WRITE(*,*)  '  b/r^4              POLARIZATION PHASE = ', phi_r4
!           WRITE(*,*)  '                            TOTAL PHASE = ', phi

         END FUNCTION phi
    
!!!................................................

!     this function calculate arg of gamm function
    REAL(dpk) FUNCTION pscoul( l, qz, qk)
      !
      USE precision, only: dpk
      !
      IMPLICIT NONE
      !
      INTEGER l
      REAL(dpk) qz, qk, q1, pi
      COMPLEX(dpk) ak, cc, z, ga, gam
      !
      q1 = l + 1
      ak = dcmplx(0.d0, qz/qk )
      ga = q1 - ak
      pi = dasin(1.d0)*2.0D+00

      cc = cdlog(ga) + cdlog(1.0D+00 + ga) + cdlog(2.0D+00 + ga) + cdlog(3.0D+00 + ga )
      cc = cc + cdlog(4.d0+ga)+cdlog(5.d0+ga)+cdlog(6.d0+ga)+cdlog(7.d0+ga)
      z  = ga + 8.d0
      gam = -cc + (z-0.5d0)*cdlog(z)-z+0.5d0*dlog(2.d0*pi)+1.d0/(12.d0*z)
      gam = gam - 1.d0/ (360.d0*z**3)+1.d0/(1260.d0*z**5)-1.d0/(1680.d0*z**7)
      pscoul = dimag(gam)
      !
    END FUNCTION pscoul
    
!!!...........................
    
  END SUBROUTINE nom
  
END MODULE norm



!!!####################################################################


 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! comment (1)
!!!
!!     essl routine computes the A^{-1} matrix 
!!

!        aa(ncs,nd)  :  matrix A
!        ncs         :  lda       , leading dimension
!        nd         :  N         , order of matrix
!         2         :  inverse & determinant are computed ,
!                                       type of computation   
!        rcond      : |S|^2 
!        det        :
!        aux(naux)  :
!        33*nd      :   size of work area for aux
!
!        on return  : A --> A^-1
!        rcond      : reciprocal of condition number
!        det        : vector with two components  D = det_1(10^det_2)


!                old code
!.................................................................
!      allocate(aux(33*nd))
!
!      call dgeicd(aa, ncs, nd, 2, rcond, det, aux, 33*nd)
!
!      deallocate(aux)
!..................................................................


!!! comment (2)
!  DGETRI computes the inverse of a matrix using the LU factorization
!  computed by DGETRF.
!
!  This method inverts U and then computes inv(A) by solving the system
!  inv(A)*L = inv(U) for inv(A).
!  Arguments
!  =========
!  N       (input) INTEGER   The order of the matrix A.  N >= 0.
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the factors L and U from the factorization
!          A = P*L*U as computed by DGETRF.
!          On exit, if INFO = 0, the inverse of the original matrix A.
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
!          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,N).
!          For optimal performance LWORK >= N*NB, where NB is
!          the optimal blocksize returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!                singular and its inverse could not be computed.
!
!  =====================================================================
!!! comment (3)
!!     dlincg 
!!     imsl routine : computes the inverse of a general complex matrix
!!  
!!     nd    :  N,      order of matrix
!!     s1x   :  A,      matrix 
!!     ncs    :  lda,    leading dimension
!!    s1xinv :  Ainv,   inverted matrix
!!     ncs    :  ldaiinv, leading dimension of the inverted matrix

!                            old code
!.....................................................................
!
!      call dlincg(nd, s1x, ncs, s1xinv, ncs)
!      sx = matmul((ux - uim*kkx), s1xinv)
!
!......................................................................

!     LAPACK ROUTINE V. 3.0
!  
!     SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
!
!     For explanation of the arguments see the double precision routine
!     above,
!
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the factors L and U from the factorization
!          A = P*L*U as computed by ZGETRF.
!          On exit, if INFO = 0, the inverse of the original matrix A.
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
!          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
!
!.... note
! 
!     here i use as working array the original one 's1xinv' as used in 
!     essl-imsl version of the same code 
!
!     ipiv is used the one used in invertin the double precision matrix
!          above
!
!

!!!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXEOF
