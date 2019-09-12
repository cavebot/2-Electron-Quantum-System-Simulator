!!##################################################################
!!!   f = r^L_(<)/r^{L+1}_(>),   r_(<) = min(x,r/2), r_(>) = max(x,r/2)  

PROGRAM h2p

  USE param
  USE set_grid
  USE one_e_matrix
  USE DATA, ONLY: write_mx, write_v_mx, write_v
  !..................
  !
  IMPLICIT NONE
  INTEGER  i
  INTEGER  l, ns                                              !
  INTEGER NOUT                                                ! I/O FILE S
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: T                 ! B - splines grid
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: B                 ! total B-splines overlap
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: H0                ! total diagonal (in l,l') hamiltonian 
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: V                 ! coupling l's potential
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: en, en_b, en_c    ! energy vectors i = 1,2,..,ns
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: ce, ce_b, ce_c    ! coefficients for state en(i)
!
  INTEGER                                :: ndiag
  CHARACTER(len=100)                     ::inputfile
  CHARACTER*100 ARGV


! executable statements !

  ndiag = 9

  NOUT = 16
  OPEN(NOUT,FILE='out/h2p.log')
!
!

!
!
  CALL GETARG(1, ARGV)              !get arguments
  READ(ARGV,*)   l0                 ! l = 0,2,4... or 1,3,5... 
  CALL GETARG(2, ARGV)              !get arguments
  READ(ARGV,*)   nl 
  CALL GETARG(3, ARGV)              !get arguments
  READ(ARGV,*)   inputfile 


  CALL input("inp/"//inputfile)            !get parameters (parameter.f90)

  nldim = nl*ndim
  

  
  !....
  !

  
  WRITE(*, '(a2,1X,a42,2i3,1X,E15.2)') "#"," h2p.f90: Rh2p l0, nl, input_file :", l0, nl, B_F
  WRITE(*, '(a2,1X,a42,i3)')    "#"," Input file: inp/h2p.inp"
  WRITE(*, '(a2,1X,a42,i3)')    "#","!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  WRITE(*, '(a2,1X,a42,i3)')    "#","           partial wave =", l0
  WRITE(*, '(a2,1X,a42,i3)')    "#"," nof partial waves   nl =", nl
  WRITE(*, '(a2,1X,a42,E15.2)') "#","     magnetic field b_f =", b_f

  !
  !..    allocations!
  !
  !


  ALLOCATE(        T(NB + KB) )  !   grid knot points 
  !
  ALLOCATE(  B( nldim, nldim) )  !  
  ALLOCATE( H0( nldim, nldim) )  !  
  ALLOCATE(  V( nldim, nldim) )  !
  !      
  ALLOCATE( EN( nldim ))         !
  ALLOCATE( CE( nldim, nldim) )  !
  !
  !
  !
  ! Construct matrix representations of operators on B-spline basis  
  !
  !      H = H_0 + V       (V: coupling part)
  !
  !      H_0 = -1/2 d^2/dr^2 + l*(l+1)/2r^2  
  !
  ! note: thus the 1-e radial Hamiltonian DOES NOT include the Coulombic potential
  !        
  !  
  
  
  CALL mkgrid(t)       ! make grid 


  b  = 0.0_dpk
  h0 = 0.0_dpk
  v  = 0.0_dpk             ! hydrogen (default)
  
  solve_the_SE:IF(problem=="h2plb") THEN
     
     CALL   MAKE_B_LB( t,  b )                ! construct b    
     CALL  MAKE_H0_LB( t, h0 )                ! construct h0
     CALL  MAKE_V_MOLECULAR_LB( t, v)         ! molecular hydrogen ion (h2+) potential
     
     !
     CALL DIAGONALIZE_DSYGVD(b, h0, v, en, ce)


  ELSE IF(problem=="h2pbl") THEN
     
     CALL            MAKE_B_BL( t,  b )        ! construct b    
     CALL           MAKE_H0_BL( t, h0 )        ! construct h0
     CALL  MAKE_V_MOLECULAR_BL( t, v)          !molecular hydrogen ion (h2+) potential
     !
     CALL DIAGONALIZE_DSBGV(b, h0, v, en, ce) !get en(i) and C(i)

  ENDIF solve_the_SE

     CALL write_v(l0,         en, "en1e-")
     CALL write_v(l0, ce(:,nldim), "ce-" )
     CALL write_v_mx(l0, en, ce,  "h1e-" )    

     CALL   ENERGY_SPECTRUM( L0, en, 'IN')             

     IF(nl.EQ.1) THEN 
        CALL WRITE_MX(L0, H0, "hmx-")
        CALL WRITE_MX(L0,  B, "bmx-")
     ENDIF

!
! solve Schrodinger Equation!  
!
!              H * C = E*B*C  
!      or    A(E)* C = C_0,        A(E) = H - E * B 
!





!!%  IF(method.EQ.'l') THEN              ! linear equations method    ( P(R) free )
!!%
!!%     ns = nbs + ncs
!!%     ALLOCATE( EN(ns) )
!!%     ALLOCATE( CE(ns, nldim) )
!!%
!!%     !     CALL  ENERGY_SPECTRUM( L0, EN, 'OUT')             ! calculate en(i) 
!!%
!!%     CALL LINEARIZE_F07BEF( H0,  B, EN, CE)            ! solve for en>0
!!%
!!%  ELSE IF(method == 'd') THEN       ! diagonalization method ( P(R) fxd  (==0) )
!!%        
!!%     ALLOCATE( EN( nldim ))
!!%     ALLOCATE( CE( nldim, nldim) )
!!%
!!%     !     CALL DIAGONALIZE_F02FHF( H, B )                 !get en(i) and C(i)
!!%
!!%     CALL DIAGONALIZE_DSBGV( B, H0, V, EN, CE)        !get en(i) and C(i)
!!%  ENDIF
!!%
!!%     CALL   ENERGY_SPECTRUM( L0, EN, 'IN')



!!%
!!%  IF(method == 'dl') THEN ! diagonalization (bound states) + linear (continuum states) method
!!%
!!%
!!%     ALLOCATE( EN_B(nldim ))
!!%     ALLOCATE( CE_B(ndim, nldim) )
!!%
!!%     CALL DIAGONALIZE_DSBGV( B, H0, V, EN_B, CE_B)          !get E(n) and C(n)
!!%
!!%     nbs = 0
!!%     DO i = 1, SIZE(en_b)                            ! get nof of b. states by fxd conditions
!!%        IF(en_b(i) < 0.0_dpk) nbs = nbs + 1
!!%     ENDDO
!!%
!!%     WRITE(*,*) '# solve_banded_diag       nbs = ', nbs
!!%
!!%     IF(spectrum=='fxd') THEN  ! note that here we need only eigenvalues of the fxd problem. 
!!%                               ! calculate coe for E_n > 0 by linear-equation solver 
!!%        ns = SIZE(en_b)        ! (results should be the same as in 'd' case
!!%
!!%        WRITE(*,*) '         ns = ', ns
!!%
!!%        ALLOCATE( EN(ns) )
!!%        ALLOCATE( CE(ns, nldim) )
!!%
!!%        en = en_b
!!%
!!%        CALL LINEARIZE_F07BEF(H0, B, EN, CE)       ! solve for E(n) > 0  (P(R) = 0) 
!!%                                                   ! and get   C(n) by LE solver 
!!%        ce(1:nbs,:)    = ce_b(1:nbs,:)
!!%
!!%        DEALLOCATE( EN_B )    !en_b
!!%        DEALLOCATE( CE_B )    !ce_b
!!%
!!%     ELSE IF (spectrum=='mxd') THEN         
!!%
!!%        ALLOCATE( EN_C(NCS) )
!!%        ALLOCATE( CE_C(NCS, NL*(NB-1)) )
!!%
!!%        CALL ENERGY_SPECTRUM(L0, EN_C, 'INOUT')         !calculate en_c(i) (do not store)
!!%        CALL LINEARIZE_F07BEF(H0, B, EN_C, CE_C)        !solve for E > 0
!!%                                                        !and get  C(E) by LE solver
!!%        ns = nbs + ncs
!!%        
!!%        WRITE(*,*) "# h1e::                       nbs = ", nbs
!!%        WRITE(*,*) "# h1e::                       ncs = ", ncs
!!%
!!%        ALLOCATE( EN(ns) )
!!%        ALLOCATE( CE(ns, nldim) )
!!%       
!!%!... assemble bound and continuum results for E,C(E)
!!%
!!%        en(     1 : nbs ) = en_b( 1 : nbs)           ! energy spectrum (b + c)
!!%        en( nbs+1 :  ns ) = en_c( 1 : ncs)
!!%        !
!!%        ce(    1: nbs,  : ) = ce_b( 1 : nbs, : )     ! coefficients    (b + c)
!!%        ce(nbs+1:  ns , : ) = ce_c( 1 : ncs, : )
!!%
!!%!...
!!%        DEALLOCATE( EN_B )    !en_b
!!%        DEALLOCATE( CE_B )    !ce_b
!!%        DEALLOCATE( EN_C )    !en_c
!!%        DEALLOCATE( CE_C )    !ce_c
!!%
!!%     ENDIF



!  DATA_FILE = "h1e-"
!  IF(code == 'bs_fx')   DATA_FILE = "d1e-"      ! P(R) == 0 
!  CALL WRITE_V_MX(L0, EN, CE, DATA_FILE)

  PRINT*, code
  IF(code == 'bs_fx') THEN                 ! P(R) == 0
 
     PRINT*, SIZE(en), SIZE(ce,dim=1), SIZE(ce,dim=2)

     CALL WRITE_V_MX(L0, en, ce, "d1e-")
  ELSE
     CALL WRITE_V_MX(L0, en, ce, "h1e-")
  ENDIF


  !........................  DEALLOCATE  ARRAYS FOR L SYMMETRY
  
  DEALLOCATE( T )
  DEALLOCATE( H0,B,V)
  DEALLOCATE( EN)
  DEALLOCATE( CE)


  WRITE(*,*) '# h1e:: Hamiltonian matrix:       h(ij)   = < B(i) | H | B(j) > : dat/hmx.dat'
  WRITE(*,*) '# h1e:: Hamiltonian matrix:       b(ij)   = < B(i) |B(j) >      : dat/bmx.dat'
!  WRITE(*,*) '# h1e:: energies+coefficients in: ',DATA_FILE

  WRITE(*,*) '# h1e:: end succesfully.'
  WRITE(*, '(a2,1X,a42,i3)') "#","!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  CLOSE(NOUT)

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  B
SUBROUTINE MAKE_B_BL(T, B)
  !
  USE utils,       only : print_mx
  USE one_e_matrix, only: f
  !
  !......................

  IMPLICIT NONE

  !declarations!

  REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: b 
  !loc
  REAL(DPK), DIMENSION(ndim,ndim)          :: mx_r0
  INTEGER                                  :: i,j,il,jl
  REAL(dpk)                                :: r0
  !executable part!....


  CALL  bsp_integral_f( nb, t,  mx_r0,  f, 0, 0, 0 )   ! < B_i | r^0 | B_j  > = <B_i|B_j>   

  WRITE(*,*) '# h1e::make_B:       call.'
!
  !  dim_h_1 = nl* SIZE(mx_r0, 1)

!...
!
  IF(print_out=='yes') THEN
     WRITE(*,*) "# subsolve::solve_banded_lin:          print B:"
     CALL print_mx(ndim, mx_r0,  'B', 'f')
  ENDIF



!#      GET THE SYMMETRIC BANDED FORM OF THE A,B MATRICES

!
! Construct the  nl*ndim x nl*ndim overlap matrix  BLMAX(il;jl) from  B(ndim,ndim)
!


  B = 0.0_dpk
  
  bsp_i:DO i = 1, ndim
     bsp_j:DO j = 1 , ndim

        r0 = mx_r0(i,j)
            
        loop_over_l:DO il = 1, nl
           !           loop_over_lp:DO jl = 1, nl

           !              IF( il == jl ) THEN

           B( (i-1)*nl+il, (j-1)*nl+il) =  r0

                 !              ENDIF


              !              ENDDO loop_over_lp
           ENDDO loop_over_l
           !
     ENDDO bsp_j
  ENDDO bsp_i


  IF(print_out=='yes') THEN
     CALL print_mx(nldim, b, 'b', 'f')
  ENDIF

END SUBROUTINE MAKE_B_BL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MAKE_H0_BL(T, H0)

  USE utils,      only : print_mx
  USE one_e_matrix, only:f
!......................

  IMPLICIT NONE

!declarations!
  REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: h0
! locals !
  REAL(DPK), DIMENSION(ndim,ndim)          :: mx_dr_0
  REAL(DPK), DIMENSION(ndim,ndim)          :: mx_r_1, mx_r_2
  INTEGER                                  :: i,j,il,jl
  REAL(dpk)                                :: la
  REAL(dpk)                                :: p2, r_2, r_1
!  external f

!!!....  
  
 
  CALL  bsp_integral_f( nb, t,  mx_dr_0,  f,  0, 1, 0 )  ! < B'_i| 1      | B'_j > kinetic     operator
  CALL  bsp_integral_f( nb, t,   mx_r_1,  f, -1, 0, 0 ) ! < B_i | x^{-1} | B_j  > coulombing  potential
  CALL  bsp_integral_f( nb, t,   mx_r_2,  f, -2, 0, 0 )  ! < B_i | x^(-2) | B_j  > centrifugal operator

!executable part!....

  WRITE(*,*) '# h1e::make_h:       call.'


!...
!
 WRITE(*,*) '# h1e::make_h:       nb = ', nb, ndim
  IF(print_out=='yes') THEN
     WRITE(*,*) "# subsolve::solve_banded_lin:          print bp,v,c:"
     CALL print_mx(ndim, mx_dr_0, 'db_1_db',   'f')
     CALL print_mx(ndim, mx_r_1,  'b_1_r_b',  'f')
     CALL print_mx(ndim, mx_r_2,  'b_1_r2_b', 'f')
  ENDIF



! Construct the  nl*ndim x nl*ndim hamiltonian         H(il;jl) from h0(ndim,ndim), O,A


  !  WRITE(*,*) '# make_h0::        size = ', dim_h_1


  H0 = 0.0_dpk
  
  bsp_i: DO i = 1, ndim
     bsp_j: DO j = 1 , ndim

        p2   = mx_dr_0(i,j)
        !        r_1  = mx_r_1(i,j)
        r_2  = mx_r_2(i,j)

        la = dble(l0)
        loop_over_li:DO il = 1, nl
           !           loop_over_lp:DO jl = 1, nl

           !                 IF( il == jl ) THEN
                    H0( (i-1)*nl+il, (j-1)*nl+il) = p2  /(2.0_dpk*ma)  + r_2 * la*(la+1)/2.0_dpk   !&
                    !                                                  - r_1 * za  * 2.0_dpk                 
                    !                 ENDIF
                 
                 !              ENDDO loop_over_lp
           la = la + 2.0_dpk
        ENDDO loop_over_li
        !
     ENDDO bsp_j
  ENDDO bsp_i

  IF(print_out=='yes') THEN
     CALL print_mx(nldim, h0, 'h0', 'f') !or change f->e for scientific output
  ENDIF

END SUBROUTINE MAKE_H0_BL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (he atomic problem) 
!SUBROUTINE MAKE_V_2E(T, V)
!REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
!REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: v
!WRITE(*,*) ' h1e:: need to be implemented 31.03.2006'
!END SUBROUTINE MAKE_V_2E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (h2+ molecular problem)
SUBROUTINE MAKE_V_MOLECULAR_BL(T, V)

  USE ang_utils
  USE ang,        ONLY : c3j, cleb
  USE utils,      only : print_mx
  USE units,      only : m_pi
  USE one_e_matrix, only:f
!................................
!

  IMPLICIT NONE

!
!declarations!
!

  REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: v
! locals !
  REAL(DPK), DIMENSION(ndim,ndim,nl) :: mx_r_lm
  REAL(DPK), DIMENSION(nl,nl,nl)     :: f_gaunt
!
  INTEGER                            :: dim_h_1
  REAL(dpk)                          :: la, lb, lm, lab 
  REAL(DPK)                          :: vlm
  INTEGER                            :: i,j,il,jl
  REAL(dpk)                          :: c3j_0, c3j_m
  REAL(dpk)                          :: zero
!executable part!

  WRITE(*,*) '# h1e::make_v_molecular:       call.'
!
!
  !  dim_h_1 = nl* SIZE(mx_r_lm, 1)

!
!...
!

  !
  !   Angular part of the molecular potential on the spherical harmonic basis  Y_(lml)(theta,phi)
  !
  !                            < Y_(lm_l) | Y_(L0)| Y_(l'm_l')> 
  !
  !
  !note 1: below the odd L numbers are excluded (for homonuclear systems)
  !        only the L = 0, 2, 4,... are included 
  !        this affects the couplings for l and l'  |l-l'| = 0,2,4,..
  !
  !note 2: check the sequence (l,l',L) in the wigner symbols below 
  !        it is correct for m_l = 0, but needs to be checked otherwise.
  !


  f_gaunt = 0.0_dpk
    c3j_0 = 0.0_dpk
    c3j_m = 0.0_dpk
  
  zero = 0.0_dpk

  la = DBLE(l0)
  loop_over_l:DO il = 1, nl

     lb = DBLE(l0)
     loop_over_lp:DO jl = 1, nl

        
        lab  = (2.0_dpk * la + 1.0_dpk) * (2.0_dpk * lb + 1.0_dpk)
!
        lm = 0.0_dpk
        loop_over_multipoles_L:DO l = 1, nl     
       
           CALL  threejsymbol(lm, la, lb, zero, zero, zero, c3j_0 )  
           CALL  threejsymbol(lm, la, lb, zero, DBLE(m_l), -DBLE(m_l), c3j_m )  


           f_gaunt(jl,il,l) = SQRT(lab) * c3j_0 * c3j_m 


!
!           WRITE(*,*) '(', int(la), int(lm), int(lb), ') <-','->(', il, l,jl, ')', m_l
!           WRITE(*,*)  c3j_0, c3j_m, f_gaunt(il,jl,l)
!
!

           lm = lm + 2.0_dpk
        ENDDO loop_over_multipoles_L
!
        lb = lb + 2.0_dpk
     ENDDO loop_over_lp
!
     la = la + 2.0_dpk
  ENDDO loop_over_l

!stop

  !
  !         Radial part of the molecular part on the B-spline basis   B_i(r)
  !
  !             U_{ll'L}(i,j) = <B_i| (r_larger)**L/(r_smaller)**(L+1)|B_j>
  !


  mx_r_lm = 0.0_dpk

  lm = 0.0_dpk
  DO l = 1, nl
     CALL  bsp_integral_f(nb, t, mx_r_lm(:,:,l), f, INT(lm), 0, 1)  !multipole expansion on B -splines
     lm = lm + 2.0_dpk
  ENDDO

  !
  ! ready calculate the molecular potential: mix radial multipoles and angular factors
  !
  !
  !   <B_iY_lm_l| V_{ll'}(r,theta,phi) | B_j Y_(l'm_l) = - 2.0 * za * Sum_L  U_{ll'L}(i,j) * <lm_l|L0|l'm_l> 
  !

  V = 0.0_dpk
  
  bsp_i: DO i = 1, ndim
     bsp_j:   DO j = 1 , ndim
        loop_over_ll:DO il = 1, nl 
           loop_over_llp:DO jl = 1, nl

              
                
              !!vlm = 0.0_dpk   
              !!DO l = 1, nl
              !!        vlm = vlm + f_gaunt(il,jl,l) * mx_r_lm(i,j,l) 
              !!ENDDO
              !! V( (i-1)*nl+il, (j-1)*nl+jl) = - 2.0_dpk * za * vlm

                 V( (i-1)*nl+il, (j-1)*nl+jl ) =  - 2.0 * za * DOT_PRODUCT(mx_r_lm(i,j,:), f_gaunt(il,jl,:))
                 
              ENDDO loop_over_llp
        ENDDO loop_over_ll
     ENDDO bsp_j
  ENDDO bsp_i
!
!
  IF(print_out=='yes') THEN
     CALL print_mx(ndim, mx_r_lm(:,:,1), 'b_1_r_b', 'f') !or change f->e for scientific output
     CALL print_mx(nldim, v, 'v', 'f')
  ENDIF

END SUBROUTINE MAKE_V_MOLECULAR_BL
!
!
!
SUBROUTINE MAKE_B_LB(T, B)
  !
  USE utils,       only : print_mx
  USE one_e_matrix, only: f
  !
  !......................

  IMPLICIT NONE

  !declarations!

  REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: b 
  !loc
  REAL(DPK), DIMENSION(ndim,ndim)          :: mx_r0
  INTEGER                                  :: i,j,il

  !executable part!....


  CALL  bsp_integral_f( nb, t,  mx_r0,  f, 0, 0, 0 )   ! < B_i | r^0 | B_j  > = <B_i|B_j>   

  WRITE(*,*) '# h1e::make_B:       call.'
!


!...
!
  IF(print_out=='yes') THEN
     WRITE(*,*) "# subsolve::solve_banded_lin:          print B:"
     CALL print_mx(ndim, mx_r0,  'B', 'f')
  ENDIF



!#      GET THE SYMMETRIC BANDED FORM OF THE A,B MATRICES

!
! Construct the  nl*ndim x nl*ndim overlap matrix  BLMAX(il;jl) from  B(ndim,ndim)
!


  B = 0.0_dpk            
  loop_over_l:DO il = 1, nl
     !     loop_over_lp:DO jl = 1, nl
        !
        bsp_i:DO i = 1, ndim
           bsp_j:DO j = 1, ndim
              
         
              !              IF( il == jl ) THEN
                 B( (il-1)*ndim+i, (il-1)*ndim+j) =  mx_r0(i,j)

                 !              ENDIF
              
           ENDDO bsp_j
        ENDDO bsp_i
        !
        !     ENDDO loop_over_lp
  ENDDO loop_over_l
  

  IF(print_out=='yes') THEN
     CALL print_mx(nldim, b, 'b', 'f')
  ENDIF

END SUBROUTINE MAKE_B_LB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MAKE_H0_LB(T, H0)

  USE utils,      only : print_mx
  USE one_e_matrix, only:f
!......................

  IMPLICIT NONE

!declarations!
  REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: h0
! locals !
  REAL(DPK), DIMENSION(ndim,ndim)          :: mx_dr_0
  REAL(DPK), DIMENSION(ndim,ndim)          :: mx_r_1, mx_r_2
  INTEGER                                  :: i,j,il
  REAL(dpk)                                :: la
  REAL(dpk)                                :: p2, r_2, r_1
!  external f

!!!....  
  
  CALL  bsp_integral_f( nb, t,  mx_dr_0,  f,  0, 1, 0 )  ! < B'_i| 1      | B'_j > kinetic     operator
  CALL  bsp_integral_f( nb, t,   mx_r_1,  f, -1, 0, 0 ) ! < B_i | x^{-1} | B_j  > coulombing  potential
  CALL  bsp_integral_f( nb, t,   mx_r_2,  f, -2, 0, 0 )  ! < B_i | x^(-2) | B_j  > centrifugal operator

!executable part!....

  WRITE(*,*) '# h1e::make_h:       call.'

!


!...
!
 WRITE(*,*) '# h1e::make_h:       nb = ', nb, ndim
  IF(print_out=='yes') THEN
     WRITE(*,*) "# subsolve::solve_banded_lin:          print bp,v,c:"
     CALL print_mx(ndim, mx_dr_0, 'db_1_db',   'f')
     CALL print_mx(ndim, mx_r_1,  'b_1_r_b',  'f')
     CALL print_mx(ndim, mx_r_2,  'b_1_r2_b', 'f')
  ENDIF



! Construct the  nl*ndim x nl*ndim hamiltonian         H(il;jl) from h0(ndim,ndim), O,A


  WRITE(*,*) '# make_h0::          l0 = ', l0

  la = DBLE(l0)
  loop_over_li:DO il = 1, nl
     !              
        bsp_i: DO i = 1, ndim
           bsp_j: DO j = 1 , ndim

              !              r_1  = mx_r_1(i,j)
              p2   = mx_dr_0(i,j)
              r_2  = mx_r_2(i,j)
              
              H0( (il-1)*ndim+i, (il-1)*ndim+j) = p2/(2.0_dpk*ma) + r_2 * la*(la+1)/2.0_dpk
                 !    - r_1 * za  * 2.0_dpk                 

              ENDDO bsp_j
           ENDDO bsp_i
           
        la = la + 2.0_dpk
     ENDDO loop_over_li



  IF(print_out=='yes') THEN
     CALL print_mx(nldim, h0, 'h0', 'f') !or change f->e for scientific output
  ENDIF

END SUBROUTINE MAKE_H0_LB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (he atomic problem) 
!SUBROUTINE MAKE_V_2E(T, V)
!REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
!REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: v
!WRITE(*,*) ' h1e:: need to be implemented 31.03.2006'
!END SUBROUTINE MAKE_V_2E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (h2+ molecular problem)
SUBROUTINE MAKE_V_MOLECULAR_LB(T, V)

  USE ang_utils
  USE ang,        ONLY : c3j, cleb
  USE utils,      only : print_mx
  USE units,      only : m_pi
  USE one_e_matrix, only:f
!................................
!

  IMPLICIT NONE

!
!declarations!
!

  REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: v
! locals !
  REAL(DPK), DIMENSION(ndim,ndim,nl) :: mx_r_lm
  REAL(DPK), DIMENSION(nl,nl,nl)     :: f_gaunt
!
  INTEGER                            :: dim_h_1
  REAL(dpk)                          :: la, lb, lm, lab 
  REAL(DPK)                          :: vlm
  INTEGER                            :: i,j,il,jl
  REAL(dpk)                          :: c3j_0, c3j_m
  REAL(dpk)                          :: zero
!executable part!

  WRITE(*,*) '# h1e::make_v_molecular:       call.'
!
!


!
!...
!

  ! make gaunt factor 

  f_gaunt = 0.0_dpk
    c3j_0 = 0.0_dpk
    c3j_m = 0.0_dpk
  
  zero = 0.0_dpk

  la = DBLE(l0)
  loop_over_l:DO il = 1, nl

     lb = DBLE(l0)
     loop_over_lp:DO jl = 1, nl

        
        lab  = (2.0_dpk * la + 1.0_dpk) * (2.0_dpk * lb + 1.0_dpk)
!
        lm = 0.0_dpk
        loop_over_multipoles_L:DO l = 1, nl     
       
           CALL  threejsymbol(la, lm, lb, zero, zero, zero, c3j_0 )  
           CALL  threejsymbol(la, lm, lb,  DBLE(m_l), zero, -DBLE(m_l), c3j_m )  


           f_gaunt(il,jl,l) = SQRT(lab) * c3j_0 * c3j_m 

!
!
!           WRITE(*,*) '(', int(la), int(lm), int(lb), ') <-','->(', il, l,jl, ')', m_l
!           WRITE(*,*)  c3j_0, c3j_m, f_gaunt(il,jl,l)

!                            * c3j( 2*la,     0, 2*lm, 0, 2*lb,   0   )    &
!                            * c3j( 2*la, 2*m_l, 2*lm, 0, 2*lb, -2*m_l )
!           WRITE(*,*)  c3j( 2*la,     0, 2*lm, 0,  2*lb,    0   )    &
!                    &, c3j( 2*la, 2*m_l, 2*lm, 0,  2*lb, -2*m_l )
!           WRITE(*,*)  cleb( 2*la,     0, 2*lm, 0,  2*lb,    0   )    &
!                    &, cleb( 2*la, 2*m_l, 2*lm, 0,  2*lb, -2*m_l )


           lm = lm + 2.0_dpk
        ENDDO loop_over_multipoles_L
!
        lb = lb + 2.0_dpk
     ENDDO loop_over_lp
!
     la = la + 2.0_dpk
  ENDDO loop_over_l

!stop

  !
  ! make molecular potential
  ! mix radial multipoles and angular factors
  !

  mx_r_lm = 0.0_dpk

  lm = 0.0_dpk
  DO l = 1, nl
     CALL  bsp_integral_f(nb, t, mx_r_lm(:,:,l), f, INT(lm), 0, 1)  !multipole expansion on B -splines
     lm = lm + 2.0_dpk
  ENDDO


  !
  !
  ! diagonal part of the molecular potential
  !
  !


!!%  loop_over_ll_d:DO il = 1, nl 
!!%        bsp_i_d: DO i = 1, ndim
!!%           bsp_j_d:   DO j = 1 , ndim
!!%              DO  l = 1, nl
!!%                 V( (il-1)*ndim+i, (il-1)*ndim+j) = -  2.0_dpk * za *mx_r_lm(i,j,l) * f_gaunt(il,il,l)  
!!%              ENDDO
!!%           ENDDO bsp_j_d
!!%        ENDDO bsp_i_d        
!!%     ENDDO loop_over_ll_d

  !
  ! 
  ! non-diagonal part of the molecular potential
  !  
  !
  loop_over_ll_v:DO il = 1, nl 
     loop_over_llp_v:DO jl = 1, nl

  
        bsp_i_v: DO i = 1, ndim
           bsp_j_v:   DO j = 1 , ndim
              
              vlm = 0.0_dpk
              !              IF( il.NE.jl ) THEN
                 
              DO  l = 1, nl
                 vlm = vlm + f_gaunt(il,jl,l) * mx_r_lm(i,j,l) 
              ENDDO

              V( (il-1)*ndim+i, (jl-1)*ndim+j) = - 2.0_dpk * za * vlm 

                 !              ENDIF
              !
           ENDDO bsp_j_v
        ENDDO bsp_i_v
        
     ENDDO loop_over_llp_v
  ENDDO loop_over_ll_v

  


!
!
  IF(print_out=='yes') THEN
     CALL print_mx(ndim, mx_r_lm(:,:,1), 'b_1_r_b', 'f') !or change f->e for scientific output
     CALL print_mx(nldim, v, 'v', 'f')
  ENDIF

END SUBROUTINE MAKE_V_MOLECULAR_LB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (Magnetized hydrogen) 
SUBROUTINE MAKE_V_MAGNETIC(T, V)

  USE ang,         only : c3j
  USE utils,       only : print_mx
  USE one_e_matrix, only:f
!......................

  IMPLICIT NONE
!
!declarations!
!
  REAL(DPK), DIMENSION(:),    INTENT(IN)   :: T
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT):: v 
!
! locals !
!
  REAL(DPK), DIMENSION(ndim,ndim)    :: mx_r0, mx_r2
  REAL(DPK), DIMENSION(nl,nl)        :: vlm
  INTEGER i,j,il,jl
  INTEGER dim_h_1
  INTEGER la, lb
  REAL(DPK) lab
  REAL(dpk) r2, r0
!  REAL(DPK) h_b_o, h_l_m, h_b_m
 ! EXTERNAL f

!executable part!....
  WRITE(*,*) '# h1e::make_v_magnetic:       call.'
!  stop

  CALL  bsp_integral_f( nb, t,   mx_r0,  f, 0, 0, 0 )   ! < B_i | x^{0} | B_j  >    
  CALL  bsp_integral_f( nb, t,   mx_r2,  f, 2, 0, 0 )    ! < B_i | x^(2) | B_j  >    

!


!...
!
  IF(print_out=='yes') THEN
     WRITE(*,*) "# subsolve::solve_banded_lin:          print O,B:"
     CALL print_mx(ndim, mx_r0,  'B', 'f')
     CALL print_mx(ndim, mx_r2,  'O', 'e')
  ENDIF


! Atoms in strong magnetic fields, Ruder, Wunner, Herold, Geyer, Springer, Berlin, 1994
! magnetic contribution to field-free hydrogenic hamiltonian 
! l's-coupling          : h_l_m
!  field-free hydrogenic states are recovered if  b_f == 0


!
!.... calculate vlm first (independent on B-splines basis)
!

  la = l0
  DO il = 1, nl
     lb = l0
     DO jl = 1, nl
        
        lab = (2.0_dpk * DBLE(la) + 1)*(2.0_dpk * DBLE(lb) + 1)
        
        vlm(il,jl) = SQRT(lab)                                     &
                   * c3j( 2*la,     0, 4, 0, 2*lb,   0   )         &
                   * c3j( 2*la, 2*INT(m_l), 4, 0, 2*lb, 2*INT(m_l) )
               
        lb = lb + 2
     ENDDO
     la = la + 2
  ENDDO

!
!...
!

  V = 0.0_dpk
  

  bsp_i:DO i = 1, ndim
     bsp_j:DO j =1 , ndim

        r0 = mx_r0(i,j)
        r2 = mx_r2(i,j)
                 
        la = INT(l0)
        loop_over_l:DO il = 1, nl
           lb = INT(l0)
           loop_over_lp:DO jl = 1, nl


              IF( il == jl ) THEN
                 V( (i-1)*nl+il, (j-1)*nl+jl) = r2 * b_f**2 /3.0_dpk   + r0 * b_f* (m_l - 1)
              ELSE

                 !  lab = (2.0D+00 * DBLE(la) + 1)*(2.0d+00 * DBLE(lb) + 1)
                 !  h_l_m = sqrt(lab)                                      &
                 !        * cleb( 2*la,     0, 4, 0, 2*lb,   0   )         &
                 !        * cleb( 2*la, 2*m_l, 4, 0, 2*lb, 2*m_l )

                 !V( (i-1)*nl+il, (j-1)*nl+jl) =  - r2 * h_l_m * b_f**2 /3.0D+00 

                 V( (i-1)*nl+il, (j-1)*nl+jl) =  - r2 * vlm(il,jl) * b_f**2 /3.0_dpk 
              ENDIF

              lb = lb + 2
           ENDDO loop_over_lp
           la = la + 2
        ENDDO loop_over_l
        !
     ENDDO bsp_j
  ENDDO bsp_i


  IF(print_out=='yes') THEN
     CALL print_mx(nldim, v, 'v', 'f')
  ENDIF
  !
END SUBROUTINE MAKE_V_MAGNETIC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DIAGONALIZE_F02FHF(H,B)
  !
  USE utils,      only : print_mx
!......................

  IMPLICIT NONE

!declarations!
  REAL(DPK), DIMENSION(:,:), intent(in):: h         ! full matrix
  REAL(DPK), DIMENSION(:,:), intent(in):: b         ! full matrix

! locals !
  INTEGER I, J
  INTEGER N, K

!F02FHF (eigenvalues only)
  REAL(DPK), DIMENSION( NL*KB, NL*(NB-2) )                        :: hb, bb    !banded matrices
  REAL(DPK), DIMENSION( MAX(NL*(nb-2),4*(NL*kB-1)*(2*NL*KB-1) ) ) :: w
  REAL(DPK), DIMENSION( NL*(NB-2) )                               :: en
  INTEGER LWORK, INFO 

!executable part!....

  WRITE(*,*) '# h1e::solve_se_eigensystem_nag:       call.'


!...

  K   = NL *  KB - 1
  N   = NL * (NB - 2)

!...

!
! hb, bb, DIMENSION(  NL*KB, NL*(NB-2)) 
!
! transform the (n-1)x(n-1) 'free' matrices to (n-2)x(n-2) 'fxd' matrices 
! equivalent to setting P(R) = 0 <=> or excluding the B_n spline.
!

  HB = 0.0_dpk
  BB = 0.0_dpk

  DO J = 1, N
     DO I = MAX(1, J - K), J
        HB( K + 1 + I - J, J ) =   H(I, J)
        BB( K + 1 + I - J, J ) =   B(I, J)
     ENDDO
  ENDDO


  !  diagonalize    :
  !  FO2FHF/NAG/M20 :     A * C = e * B * C,   where, A,B banded real symmetric
  !


  INFO = 0
  LWORK = MAX( NL*(nb-2), 4*(NL*kB-1)*(2*NL*KB-1))    

!...

  CALL F02FHF(N, K, HB, K+1, K, BB, K+1, EN, W, LWORK, INFO)

!...
      
  IF(INFO.NE.0) THEN 
     
     WRITE(*,*) '# ERROR IN  F02FHF/NAG/M20 SUBROUTINE,    INFO = ', INFO 
     WRITE(*,*) '#'
     STOP
     
  ELSE  

     OPEN(10, file='out/en-nag.out') 
     DO  I = NL*(NB-2), 1, -1 
        
        WRITE(10,"(I4,3X,3E20.10)") I, EN(I)
        
     ENDDO
     CLOSE(10)
  ENDIF
      
END SUBROUTINE DIAGONALIZE_F02FHF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DIAGONALIZE_DSBGV(B, H0, V, EN, CE)

  USE utils,       ONLY : print_mx
  !declarations!
  IMPLICIT NONE
  !
  REAL(DPK), DIMENSION(:,:),     INTENT(IN) :: H0
  REAL(DPK), DIMENSION(:,:),     INTENT(IN) :: B
  REAL(DPK), DIMENSION(:,:),     INTENT(IN) :: V
  REAL(DPK), DIMENSION(:,:),  INTENT(INOUT) :: CE
  REAL(DPK), DIMENSION(:),    INTENT(INOUT) :: EN

  !locals!
  INTEGER I, J
  INTEGER N, K

  ! DSBGV/LAPACK
  REAL(DPK), DIMENSION(  NL*KB, nldim )   :: HB, BB
  REAL(DPK), DIMENSION(  nldim, nldim)    :: Z      !dsbgv
  REAL(DPK), DIMENSION( 3*nldim )         :: w      !dsbgv
  INTEGER INFO 

!executable part!

!...

  WRITE(*,*) '# h1e::diagonalize_dsbgv:   call.'

!...

  IF(print_out=='yes') THEN
     WRITE(*,*) "# subsolve::solve_banded_lin:          print    b:"
     CALL print_mx(nl*ndim, b,  'b', 'f')
     WRITE(*,*) "# subsolve::solve_banded_lin:          print   h0:"
     CALL print_mx(nl*ndim, h0,  'h0', 'f')
     WRITE(*,*) "# subsolve::solve_banded_lin:          print    v:"
     CALL print_mx(nl*ndim, v,  'v', 'f')
     !     WRITE(*,*) "# subsolve::solve_banded_lin:          print h0_v:"
     !     CALL print_mx(nl*ndim, h0_v,  'h0_v', 'f')
  ENDIF



  K   = NL *  KB - 1



  WRITE(*,*) '# h1e::diagonalize_dsbgv:         k = ', k
  WRITE(*,*) '# h1e::diagonalize_dsbgv:   nl*ndim = ', nldim

!............................

! hb, bb, DIMENSION(  NL*KB, NL*(NB-2)) 
!
! transform the (n-1)x(n-1) 'free' matrices to (n-2)x(n-2) 'fxd' matrices 
! equivalent to setting P(R) = 0 <=> or excluding the B_n spline.
!

      HB = 0.0_dpk
      BB = 0.0_dpk

      DO J = 1, nldim
         DO I = MAX(1, J - K), J
            HB( K + 1 + I - J, J ) =   H0(I, J) + V(I,J)
            BB( K + 1 + I - J, J ) =   B(I, J) 
         ENDDO
      ENDDO

      !
      !  DSBGV/LAPACK :     A * C = e * B * C,    where, A,B banded real symmetric
      !
      !   EV NORMALIZED AS : Z^T * C * Z = 1
      !   EN(I) <====> Z(J,I)
      !

      INFO = 0

      CALL DSBGV('V', 'U', nldim, k, k, hb, k+1, bb, k+1, en, z, nldim, w, info)
      
      !...
      IF(info.NE.0) THEN 
         WRITE(*,*) '# ERROR IN  DSBGV/LAPACK SUBROUTINE,    INFO = ', info 
         WRITE(*,*) '#'
         STOP
      ELSE  

         CE = TRANSPOSE(Z)                 !assign Z(NB-2,NB-2) ---> CE(NB-2,NDIM)


         !         DO i = 1, SIZE(z,dim=1)                ! force solutions with dP(0)/dr > 0
         !  IF(ce(i,2).LT.0.0_dpk) THEN
         !     ce(i,:) = -ce(i,:)
         !  END IF
         !ENDDO
         

         DO i = 1, ndim
            WRITE(11,*), i, ( ce(i,j), j= ndim, nldim,ndim)
         ENDDO
         !         stop
         DO i = 1, 20
            WRITE(*,*) i, en(i)
         ENDDO


         !         DO i = nldim, 1, -1
         !            WRITE(*,*) i, en(i)
         !ENDDO

      ENDIF
     
      
    END SUBROUTINE DIAGONALIZE_DSBGV
!
!
!  A * X = lambda * B * X 
!
!
    SUBROUTINE DIAGONALIZE_DSYGVD(B, H0, V, EN, CE)

      USE utils,       ONLY : print_mx
      !declarations!
      IMPLICIT NONE
      !
      REAL(DPK), DIMENSION(:,:),     INTENT(IN) :: H0
      REAL(DPK), DIMENSION(:,:),     INTENT(IN) :: B
      REAL(DPK), DIMENSION(:,:),     INTENT(IN) :: V

      REAL(DPK), DIMENSION(:,:),  INTENT(INOUT) :: CE
      REAL(DPK), DIMENSION(:),    INTENT(INOUT) :: EN
      
      !locals!
      INTEGER I, J
      INTEGER N

      
      ! DSYGVD/LAPACK
      REAL(DPK), DIMENSION( nldim, nldim )                  :: hb, bb
      REAL(DPK), DIMENSION( 1 + (64 + 6 + 2*nldim) * nldim) :: work      !
      INTEGER,   DIMENSION( 3+5*nldim )                     :: iwork     !
      !      REAL(dpk), DIMENSION( nldim,nldim)                    :: h0_v
      INTEGER                                               :: info, ki, k, itype
      
      

      !executable part!
      
      !...
      
      WRITE(*,*) '# h1e::diagonalize_dsygd:   call.'
      
      !...
      
      !      h0_v = h0 + v
      
      IF(print_out=='yes') THEN
         WRITE(*,*) "# subsolve::solve_banded_lin:          print    b:"
         CALL print_mx(nl*ndim, b,  'b', 'f')
         WRITE(*,*) "# subsolve::solve_banded_lin:          print   h0:"
         CALL print_mx(nl*ndim, h0,  'h0', 'f')
         WRITE(*,*) "# subsolve::solve_banded_lin:          print    v:"
         CALL print_mx(nl*ndim, v,  'v', 'f')
         !         WRITE(*,*) "# subsolve::solve_banded_lin:          print h0_v:"
         !         CALL print_mx(nl*ndim, h0_v,  'h0_v', 'f')
      ENDIF
      
      !      n  = NL * NDIM      
      k  = 1 + ((64 + 6 + 2*nldim) * nldim )
      ki =  3 + 5 * nldim    !liwork      


      
      WRITE(*,*) '# h1e::diagonalize_dsygvd:      ndim = ',  ndim
      WRITE(*,*) '# h1e::diagonalize_dsygvd:   nl*ndim = ', nldim
      
      !............................
      
      ! hb, bb, DIMENSION(  NL*KB, NL*(NB-2)) 
      !
      ! transform the (n-1)x(n-1) 'free' matrices to (n-2)x(n-2) 'fxd' matrices 
      ! equivalent to setting P(R) = 0 <=> or excluding the B_n spline.
      !
      
      HB = 0.0_dpk
      BB = 0.0_dpk
      
      DO i = 1, nldim
         DO j = i, nldim
            HB( i, j ) =   h0(i, j) + v(i,j)
            BB( i, j ) =   b(i, j) 
         ENDDO
      ENDDO
      
      !
      !  DSYGVD/LAPACK :     A * C = e * B * C,    where, A,B banded real symmetric
      !
      !   EV NORMALIZED AS : Z^T * C * Z = 1
      !   EN(I) <====> Z(J,I)
      !
      


      itype = 1 ! a * x = lambda * b * x
      info  = 0
      
      CALL DSYGVD(1,'V','U',nldim, hb, nldim, bb, nldim, en, work, k, iwork, ki, info)

      
      !...
      IF(info.NE.0) THEN 
         !
         WRITE(*,*) '# ERROR IN  DSYGVD/LAPACK SUBROUTINE,    INFO = ', info
         WRITE(*,*) '#'
         STOP
         !
      ELSE
         !

         OPEN(10, file='out/en-dsbgv.out') 
         DO  I = 1 , nldim            
            WRITE(10,"(I4,3X,3E20.10)") i, en(i)            
         ENDDO
         CLOSE(10)


         DO i = 1, SIZE(ce,dim=1)                ! force solutions with dP(0)/dr > 0
            IF(ce(i,2).LT.0.0_dpk) THEN
               ce(i,:) = -ce(i,:)
            END IF
         ENDDO

         
         CE = TRANSPOSE(hb)                 !assign Z(NB-2,NB-2) ---> CE(NB-2,NDIM)


         DO i = 1, ndim
            WRITE(11,*), i, ( ce(i,j)**2, j= ndim, nldim, ndim)
         ENDDO
!         stop
         DO i = 1, 20
            WRITE(*,*) i, en(i)
         ENDDO
         
      ENDIF
      
      
    END SUBROUTINE DIAGONALIZE_DSYGVD
    !
!
!
    SUBROUTINE LINEARIZE_F07BEF(H, B, EN, CE)  
!  USE nag_f77_f_chapter
  USE utils, ONLY : check
!  USE utils, ONLY : print_mx, check

!declarations!
  IMPLICIT NONE
  REAL(DPK), DIMENSION(:,:),           INTENT(IN) :: H
  REAL(DPK), DIMENSION(:,:),           INTENT(IN) :: B
  REAL(DPK), DIMENSION(:),          INTENT(inout) :: en
  REAL(DPK), DIMENSION(:,:),        INTENT(inout) :: ce

! locals !
  INTEGER i, j
!F07BEF!
  INTEGER K, N
  REAL(DPK), DIMENSION( 3*NL*KB-2, nldim) :: A_E  ! linear
  REAL(DPK), DIMENSION( nldim,    1     ) :: C0
  INTEGER,   DIMENSION( nldim )           :: IW
  INTEGER INFO
!
  REAL(dpk) norm
  INTEGER l, ie
  INTEGER nout
  INTEGER dim, dim_1, dim_2  
  integer er

!executable!
!...
  WRITE(*,*) '# h1e::linearize_f07bef:                  call.'

!...

  NOUT     = 12

!...

  K  = NL *  KB - 1          ! number of sub-diagonials (upper diagonals)

!...

  er = 1


!#      GET THE SYMMETRIC BANDED FORM OF THE A,B MATRICES


  OPEN(NOUT, FILE='out/subsolve_lin.out')

  WRITE(nout, *) '#       LINEAR SOLUTION : '

!...


  dim = size(en)
  dim_1 = SIZE(ce,1)
  dim_2 = SIZE(h, 1)

  WRITE(*,*) '# nof states          ns = ', dim_1, nbs + ncs 
  WRITE(*,*) '# nof basis           nc = ', dim_2, nl * ( nb - 1 )

  CALL check(dim_1, dim, er)
  CALL check(dim_2, nl*ndim, er+1)

    
  DO ie = 1, dim                      ! solve for 'ns' states: en(ie), ie = 1, ..., ns

!...           A(E) = H - E*B 
 
     DO J = 1, nldim
        DO I = MAX(1, J - K), MIN(N, J + K)
           
           A_E(2*K + I - J + 1, J) = H(I,J) - EN(IE) * B(I,J)         

        ENDDO
     ENDDO


!... solve     A(E) * C = C0,

     IW = 0
     CALL F07BDF(nldim, nldim, k, k, A_E, 3*nl*kb-2, iw, info)      ! factorize

     IF(INFO.NE.0) THEN 

        WRITE(*,*) '# subsolve:: f07bdf: Error in factorization routine.'
        WRITE(*,*) '# subsolve:: f07bdf:                         info = ', info

        CALL X04CEF(nldim, nldim, k, 2*k, a_e, 3*nl*kb-2,'X04CEF: A=PLU, U:', info)

     ELSE

        C0 = 0.0_dpk                        ! initialize C0
        C0(NDIM, 1) = 1.0_dpk

        IF(ie.EQ.1) THEN
           WRITE(*,*) '# c0 = ', ie,  SIZE(c0, 1), SIZE(c0, 2) 
        ENDIF


        CALL F07BEF('N', nldim, k, k, 1, a_e, 3*nl*kb-2, iw, c0, nldim, info)

          ENDIF


          IF(INFO.NE.0) THEN 
             WRITE(nout,*) '#         ERROR IN   F07BEF          '
             WRITE(nout,*) '#'
             WRITE(nout,*) '#                     INFO = ', INFO
             WRITE(nout,*) '#'
             STOP
          ELSE

             WRITE(nout, *) '#   INPUT     ENERGY      E = ', EN(IE)
             WRITE(nout, *) '#'
             
                          
             NORM = DOT_PRODUCT( C0(:,1), MATMUL( B, C0(:, 1) ))


!             CE(IE,:) = C(:,1) /SQRT(NORM) 
             

             ! transform output vector from c_i(l=0,1,..,nl) ---> c_l(i=1,2,..,ndim)

             DO l = 1, nl
                DO j = 1, ndim 

                   CE(ie,j +(l-1)*ndim) = C0(l +(j-1)*nl,1) /DSQRT(NORM) 
                ENDDO
             ENDDO

          ENDIF

          WRITE(NOUT, *) IE, EN(IE) 

          IF(IE.EQ.30.AND.L.EQ.0) THEN
             DO  J = 1, NDIM/2
               
                WRITE(NOUT, "(I4,1X,E25.14,2X,I4,1X,E25.14)") J, CE(IE,J), J+NB/2, CE(IE, J+NB/2)
             ENDDO
          ENDIF

       ENDDO

       CLOSE(NOUT)
       
     END SUBROUTINE LINEARIZE_F07BEF

     !
     !
     !
     !


    SUBROUTINE ENERGY_SPECTRUM(L, EN, MODE)
 
      USE param
      USE DATA, ONLY: write_v
      
      IMPLICIT NONE

      INTEGER L
      REAL(DPK), DIMENSION(:), INTENT(INOUT) :: en
      CHARACTER(LEN=*)  MODE
!... locals
      INTEGER I

!.....................

      WRITE(*,*) '# energy_spectrum:    size of en = ', size(en)

      
      IF (MODE =='OUT') THEN 

         DO  I = L + 1 , SIZE(en)
            
            IF(I.LE.NBS + L) THEN 

               EN(I) = - 0.5D+00 * (za/I)**2

            ELSE IF(I.LE.SIZE(en)) THEN 
               
               EN(I) =   EN_1E  +    DE_1E * ( I - NBS - 1 - L) 
     
            ENDIF

         ENDDO

         CALL write_v(l, en, "en1e-")

      ELSE IF (MODE =='INOUT') THEN      ! do not store

         DO  I = 1 , SIZE(en)
            EN(I) =   EN_1E  +    DE_1E * ( I - 1) 
         ENDDO

      ELSE IF (MODE == 'IN') THEN 
         CALL write_v(l, en, "en1e-")
      ENDIF
    END SUBROUTINE ENERGY_SPECTRUM



  END PROGRAM H2P
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EOF
!###################################################################
! laan jan/feb 2003 iesl day  : bsp1ef/h1e.f90
!
!
!   Author         :   Lambros Nikolopoulos
!                  :
!!    JAN/FEB 2003 :
!                  :
!   I/0 files      :   inp/h1e.inp 
!                  :   out/h1e.out
!                  :   out/setmat.log
!                  :   out/subsolve.lin.out
!                  :   out/subsolve.inv.out
!                  :
!   I/O data files :   DAT/HD1E-L.LINDAT,   L = 0, 1,  ..., NL
!                  :   DAT/HD1E-L.INV.DAT,  L = 0, 1,  ..., NL
!                  :   DAT/HMX-L.DAT,        L = 0, 1, ..., NL
!                  :   DAT/BMX-L.DAT,        L = 0, 1, ..., NL
!                  :   dat/knot.dat 
!   LIBRARIES      :   NAG/FO2SDF   
!                  :   NAG/F02FHF 
!                  :   NAG/F07BDF, NAG/F07BEF
!                  :
!   SUB/MOD        :  modules.f90, sub_vhf.f, modio.f90
! 
!
!
!    h1e.f90 :
!               1. Calculates energies + coefficients for l = 0, 1, 2,..,nl
!                  for 1-e atomic hydrogenic systems with free boundary 
!                  conditions at the end of the interval [0, R].
!                  Two methods have been employed :
!                  a. inverse iteration              (FO2SDF)
!                  b. solution of linear systems
!
!               2. Calculates and stores the 1-e Hamiltonian on Bsplines 
!                  as well as the B-splines overlap.
!
!                  H_ij = < B_i | -d^2/dr^2 + V(r) | B_j >
!                  B_ij = < B_i | B_j >
!
!  2003.02.05
!
!      to be done: store the matrices in banded form directly
!                  instead of calling full and then transform them into
!                  banded matrices
!
!
! nbs : nof bound states
! ncs : nof continuum states
!



