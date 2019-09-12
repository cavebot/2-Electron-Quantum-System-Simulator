
!
!    MODULE WF_1E : 
!
!    SUBROUTINES :
!                  WRITE_WF1E, READ_WF1E, OVERLAP, SAVE_WF1E_NL  
!
MODULE WF_1E
  !
  USE PRECISION, ONLY: DPK
  PUBLIC read_wf1e, store_radial_psi_bs, evaluate_dipoles, evaluate_dipoles_diatomics,  read_partial_waves_l
  !
CONTAINS

  !.......................
  SUBROUTINE WRITE_WF1E(nplot, l, en, ce)
      !
     USE param
     USE set_grid
     USE one_e_matrix
     USE utils, ONLY: datafile
     
     !......................
     
     IMPLICIT NONE
     ! arguments
     INTEGER, OPTIONAL,INTENT(in) :: nplot
     INTEGER                      :: l         ! angular number
     REAL(DPK), DIMENSION(:)      :: en        ! energy spectrum
     REAL(DPK), DIMENSION(:,:)    :: ce        ! coefficients
     ! integer ns
     ! locals
     INTEGER                       ::  i, j, ir, jj
     INTEGER                       :: file
     REAL(DPK), DIMENSION(NB)      :: bcoeff
     REAL(DPK), DIMENSION(NB+KX)   :: t     !  r_i, knot sequence
     REAL(DPK), DIMENSION(NP)      :: p     !  P(r_i)
     REAL(DPK), DIMENSION(NP)      :: dp    !  dP(r_i)/dr
     REAL(DPK)                     :: a_s
     REAL(DPK)                     :: bvalue
     REAL(dpk)                     :: sum_boundary_value
     !
     INTEGER                       :: n_store
     INTEGER                       :: nout
     !..........................


     IF(PRESENT(nplot)) THEN
        n_store = nplot
        nout = 16
     ELSE
        n_store = 0
     ENDIF

     FILE = 56
     
     !.........................
     
     CALL r_grid
     CALL mkgrid(t)
     CALL datafile(file, l,"w1e-")
     
     WRITE(FILE) SIZE(en)
     WRITE(FILE) no, hh

     WRITE(*,*) '# write_wf1e::    nof wf points    (np,no) =  ',  np,  no
     WRITE(*,*) '# write_wf1e::    nof wf points  (ndim,nb) =  ',  ndim, nb
     WRITE(FILE) ( r(j), j = 1, no )
     WRITE(FILE) (dr(j), j = 1, no )
     


     !...................
     sum_boundary_value = 0.0_dpk

     bcoeff = 0.0_dpk
     loop_eigenstates: do i = 1, size(en)

        bcoeff(1) = 0.0_dpk
        DO jj = 1, ndim !nb - 1 
           bcoeff(jj + 1) = ce(i, jj)
        END DO
        
        P(1) = 0.0_dpk
        DO J = 2, NO
            P(J) =  BVALUE( T, BCOEFF, NB, KB, R(J), 0 )
            DP(J) = BVALUE( T, BCOEFF, NB, KB, R(J), 1 )
        END DO

        sum_boundary_value = sum_boundary_value + p(no)*p(no)
        
        
        !     IF(last_spline==0) THEN     ! f(R) = 0
        A_S = 0.5_dpk * ( KB - 1 ) / ( t(NB + 1 ) - t(NB) )
        DP(NO) = 2 * A_S * ( BCOEFF(NB) - BCOEFF(NB-1) )
        !     ENDIF
        
        !        WRITE(*,*) i, bcoeff(nb-1)/bcoeff(nb)
        !.....................
        
        !  save radial wavefunction    P_nl( r)  [0, R]
        
        WRITE(FILE)  EN(I) 
        WRITE(FILE) ( P(J), J = 1, NO)
        WRITE(FILE)( DP(J), J = 1, NO)

        store_wf:IF(n_store==i) THEN 
           OPEN(nout,file='out/wf1e.out')
           DO ir = 1, np 
              WRITE(nout,'(3E20.10)') r(ir),p(ir),dp(ir)
           ENDDO
           CLOSE(nout)
        ENDIF store_wf
        
     ENDDO loop_eigenstates
     
     WRITE(*,*) '# write_wf1e::  S_n|p(b)|^2 = ', sum_boundary_value

     CLOSE(FILE) 
     
!!!..............................

   END SUBROUTINE WRITE_WF1E
!
!
!
!#################################################
    SUBROUTINE READ_WF1E(L, PL, DPL, EN, R, DR)
  
      USE units, ONLY : dpk  
      USE utils, ONLY : datafile
      
      IMPLICIT NONE
      
      INTEGER L
      INTEGER FILE
      REAL(DPK), DIMENSION(:,:), POINTER :: pl,dpl       !  P_kl(r_j), dP_kl(r_j)/dr
      REAL(DPK), DIMENSION(:),   POINTER :: r,dr         !  r_j
      REAL(DPK), DIMENSION(:),   POINTER :: en           !  E_i
      INTEGER                            :: ns
      INTEGER                            :: np
      REAL(DPK)                          :: hh           !  dr_j
      !
      INTEGER I, J
      !...............


      FILE = 56
      
      CALL DATAFILE(FILE,L,"wf1e-")

      READ(FILE) NS
      READ(FILE) NP, HH
      
!      WRITE(*,*) 'ns = ', ns 
      ALLOCATE( EN( NS ) )
      ALLOCATE(  PL(NS, NP) )
      ALLOCATE( DPL(NS, NP) )
      ALLOCATE(   R(NP) )
      ALLOCATE(   DR(NP) )
      
      !now read
      READ(FILE) (R(J),    J = 1, NP)
      READ(FILE) (DR(J),   J = 1, NP)      
      DO I = 1, NS
         READ(FILE)   EN(I) 
         READ(FILE) ( PL(I,J), J = 1, NP)
         READ(FILE) (DPL(I,J), J = 1, NP)
      ENDDO
      CLOSE(FILE)     
      !
    END SUBROUTINE READ_WF1E
!
!
!
   SUBROUTINE WRITE_ALL_WF1E(nplot, l, en, ce)
     !
     USE param
     USE set_grid
     USE one_e_matrix
     USE utils, ONLY: datafile
     !
     !...........................................................
     !
     IMPLICIT NONE
     ! arguments
     INTEGER, OPTIONAL,INTENT(in) :: nplot
     INTEGER                      :: l         ! angular number
     REAL(DPK), DIMENSION(:)      :: en        ! energy spectrum
     REAL(DPK), DIMENSION(:,:)    :: ce        ! coefficients
     ! integer ns
     ! locals
     INTEGER                       :: i, j, ir, jj
     INTEGER                       :: ifile
     REAL(DPK), DIMENSION(nb)      :: bcoeff
     REAL(DPK), DIMENSION(nb+kx)   :: t        !  r_i, knot sequence
     REAL(DPK), DIMENSION(np)      :: p        !  P(r_i)
     REAL(DPK)                     :: bvalue
     integer                       :: n_store
     integer                       :: nout
     !

     !
     IF(PRESENT(nplot)) THEN
        n_store = nplot
        nout = 16
     ELSE
        n_store = 0
     ENDIF

     ifile = 56
     
     !.........................
     
     CALL r_grid
     CALL mkgrid(t)
     CALL datafile(ifile, l,"wf1e-")

     bcoeff = 0.0_dpk
     p      = 0.0_dpk

     WRITE(ifile) SIZE(en), np
     loop_eigenstates: DO i = 1, SIZE(en)
        !
        bcoeff(1) = 0.0_dpk
        DO jj = 1, ndim !nb-1
           bcoeff(jj + 1) = ce(i, jj)
        END DO
        

        !
        DO j = 2, np
           p(j) = bvalue(t, bcoeff, nb, kb, r(j), 0 )
        END DO
        
        !  save radial wavefunction    P_nl( r)  [0, R]
        
         WRITE(ifile) (p(j), j = 1, np)

         store_wf:IF(n_store==i) THEN 
            OPEN(nout,file='out/wf1e.out')
            DO ir = 1, np
               WRITE(nout,'(3E20.10)') r(ir), p(ir)
            ENDDO
            CLOSE(nout)
         ENDIF store_wf
         
      ENDDO loop_eigenstates
      WRITE(ifile) ( r(ir), ir = 1, np)
     
     CLOSE(ifile)
   END SUBROUTINE WRITE_ALL_WF1E
   !
   !
   !
   !
   SUBROUTINE WRITE_LOPT_WF1E(nplot, l, en, ce)
     !
     USE param
     USE set_grid
     USE one_e_matrix
     USE utils, ONLY: datafile
     USE units, only:m_pi
     !......................
     
     IMPLICIT NONE
     ! arguments
     INTEGER, OPTIONAL,INTENT(in) :: nplot
     INTEGER                      :: l         ! angular number
     REAL(DPK), DIMENSION(:)      :: en        ! energy spectrum
     REAL(DPK), DIMENSION(:,:)    :: ce        ! coefficients
     ! integer ns
     ! locals
     INTEGER                       :: i, j, ir, jj
     INTEGER                       :: ifile
     REAL(DPK), DIMENSION(nb)      :: bcoeff
     REAL(DPK), DIMENSION(nb+kx)   :: t        !  r_i, knot sequence
     REAL(DPK), DIMENSION(np,nl)   :: p        !  P(r_i)
     REAL(DPK), DIMENSION(nl)      :: c_l        !  P(r_i)
     REAL(dpk)                     :: c_0, f_l
     REAL(DPK)                     :: bvalue
     INTEGER                       :: n_store
     INTEGER                       :: nout


     !
     IF(PRESENT(nplot)) THEN
        n_store = nplot
        nout = 16
     ELSE
        n_store = 0
     ENDIF

     ifile = 56
     
     !.........................
     
     CALL r_grid
     CALL mkgrid(t)
     CALL datafile(ifile, l,"w1e-")


     c_0 = 0.5_dpk/SQRT(m_pi)




     
     WRITE(*,*) '# write_lopt_wf1e::    nof wf points    (np,no) =  ',  np,  no
     WRITE(*,*) '# write_lopt_wf1e::                   (ndim,nb) =  ',  ndim, nb
     
     !     bcoeff = 0.0_dpk
     !     p      = 0.0_dpk

     !     WRITE(ifile) SIZE(en), np
     loop_eigenstates: DO i = 1, SIZE(en)
        !
        bcoeff = 0.0_dpk
        p = 0.0_dpk
        loop_partial_waves: DO l = 1, nl !note that for diatomics l= 0,2,4,... or l=1,3,5,..

           c_l(l) = c_0 * SQRT(4.0_dpk*DBLE(l) - 3.0_dpk )    !c_1 = c_0, c_2 = c_0 *sqrt(5), c_3=c_0*sqrt(9)

           bcoeff(1) = 0.0_dpk 
           DO jj = 1, ndim 
              bcoeff(jj + 1) = ce(i, jj + (l-1)*ndim) ! ce(i, jj)

              IF(jj.EQ.ndim) THEN
                 PRINT*, i, bcoeff(nb-1), bcoeff(nb-2), r(np-1),r(np)
              ENDIF
           END DO
           

           !
           DO j = 2, np
              p(j,l) = bvalue(t, bcoeff, nb, kb, r(j), 0 )
           END DO

           !save radial wavefunction    P_nl( r)  [0, R]
           !   WRITE(ifile) (p(j), j = 1, np)
           !   store_wf:IF(n_store==i) THEN 
           !   OPEN(nout,file='out/wf1e.out')
           !   DO ir = 1, np
           !      WRITE(nout,'(2E20.10,5x,i5)') r(ir), p(ir)/r(ir)
           !   ENDDO
           !   CLOSE(nout)
           !ENDIF store_wf

        ENDDO loop_partial_waves

        WRITE(ifile) (p(1:np,l), l = 1, nl) !save radial wavefunction    P_nl( r)  [0, R]



        store_wf:IF(n_store==i) THEN 
           OPEN(nout,file='out/wf1e.out')
           DO ir = 1, np

              f_l = 0.0_dpk              
              DO l = 1, nl
                 f_l = f_l + c_l(l) * p(ir,l) !/ r(ir)
              ENDDO
! fails for the linear grid where               
!              WRITE(nout,'(20E15.5)') r(ir),f_l, (c_l(l)*p(ir,l)/r(ir), l=1,nl)
              WRITE(nout,'(20E15.5)') r(ir),f_l, (c_l(l)*p(ir,l), l=1,nl)
           ENDDO
           CLOSE(nout)
        ENDIF store_wf

        !
     ENDDO loop_eigenstates
      !      WRITE(ifile) ( r(ir), ir = 1, np)
     
     CLOSE(ifile)
   END SUBROUTINE WRITE_LOPT_WF1E




!
! reads r,  p_l(r), n =1,...,ns for fixed l
!  
!
!
   !
   SUBROUTINE read_partial_waves_l(l, pl, r)
     !
     USE utils, ONLY: datafile
     !......................
     
     IMPLICIT NONE
     ! arguments
     INTEGER                   :: l         ! angular number
     REAL(DPK), DIMENSION(:,:) :: pl
     REAL(DPK), DIMENSION(:)   :: r
     INTEGER                   :: ns, np
     INTEGER                   :: i, j
     INTEGER                   :: ifile
     !
     ifile = 1
     
     !.........................
     

     CALL datafile(ifile, l,"wf1e-")    

     READ(ifile) ns, np


     IF(ns.NE.SIZE(pl,dim=1)) THEN 
        PRINT*,"# read_pl_r::  wrong input file                  ns =  ", ns
        PRINT*,"# read_pl_r::                     nof energy states =  ", SIZE(pl,dim=1)
        STOP
     ENDIF
     IF(np.NE.SIZE(pl,dim=2)) THEN 
        PRINT*,"# read_pl_r::  wrong input file                  np =  ", np
        PRINT*,"# read_pl_r::                             grid size =  ", SIZE(pl,dim=2)
        STOP
     ENDIF

     
     pl     = 0.0_dpk
     r      = 0.0_dpk
     loop_eigenstates: DO i = 1, SIZE(pl,dim=1)
        READ(ifile) ( pl(i, j), j = 1, SIZE(pl,dim=2) )    ! read radial wavefunction    P_nl( r_j )  [0, R]
     ENDDO loop_eigenstates
     READ(ifile) (r(i), i = 1, SIZE(r))            ! read r(j)
     
     CLOSE(ifile)
     !
   END SUBROUTINE READ_PARTIAL_WAVES_L

!
!
!
   SUBROUTINE READ_ALL_WF1E(l, p, r)
     !
     USE utils, ONLY: datafile
     !......................
     
     IMPLICIT NONE
     ! arguments
     INTEGER                                :: l         ! angular number
     REAL(DPK), DIMENSION(:,:),     POINTER :: p
     REAL(DPK), DIMENSION(:),       POINTER :: r
     INTEGER                                :: ns, np
     INTEGER                                :: i, j
     INTEGER                                :: ifile
     !
     ifile = 1
     
     !.........................
     

     CALL datafile(ifile, l,"wf1e-")

     READ(ifile) ns, np

     ALLOCATE( p(ns,np) )
     ALLOCATE( r(   np) )
     
     p      = 0.0_dpk
     r      = 0.0_dpk
     loop_eigenstates: DO i = 1, ns
        READ(ifile) ( p(i, j), j = 1, np )    ! read radial wavefunction    P_nl( r_j )  [0, R]
     ENDDO loop_eigenstates
     READ(ifile) (r(i), i = 1, np)            ! read r(j)

     CLOSE(ifile)
     !
   END SUBROUTINE READ_ALL_WF1E
   !

   !
   !
    !#################################################
    SUBROUTINE READ_BASIS_GRID(l,r)
      !
      USE units, ONLY : dpk  
      USE utils, ONLY : datafile
      !
      IMPLICIT NONE
      !
      INTEGER L
      INTEGER FILE
      REAL(DPK), DIMENSION(:),   POINTER :: r       !  r_j
      !local
      INTEGER                            :: ns
      INTEGER                            :: np
      REAL(DPK)                          :: hh      !  dr_j
      !
      INTEGER I
      !...............

      FILE = 56
      
      CALL DATAFILE(FILE,L,"hwf1e-")
      
      READ(FILE) NS
      READ(FILE) NP, HH
      ALLOCATE(  R(NP) ) ;  
      READ(FILE) (R(I),    I = 1, NP)
      RETURN
    END SUBROUTINE READ_BASIS_GRID
!
!
!!$    SUBROUTINE OVERLAP(L, P1, P2, DR, NP, INTEGRAL)
!!$  
!!$      USE units,       ONLY : dpk
!!$      IMPLICIT NONE
!!$      
!!$      INTEGER L
!!$      REAL(DPK), DIMENSION(np) :: p1,p2       ! radial functions to be overlaped 
!!$      REAL(DPK), DIMENSION(np) :: dr        ! r_i, dr_i
!!$      INTEGER                  :: np          ! nof grid points
!!$      REAL(dpk),   INTENT(out) :: integral
!!$      !
!!$      INTEGER FILE
!!$      INTEGER I, J
!!$      !...............
!!$
!!$      integral = 0.0_dpk
!!$      DO i=1, np                              !calculate <P1|P2> 
!!$         integral  = integral + p1(i) * p2(i) * dr(i)
!!$      ENDDO
!!$      !      
!!$      RETURN
!!$    END SUBROUTINE OVERLAP

    !S
    !S evaluates <P_(na,la) | K_ab * T_G(a,b) | P_(nb,kb) >
    !S
    !S  na = 1, ..., na_max           G = l,v,a
    !S  nb = 1,...., nb_max
    !S
    !S corrections are added  
    !S

    SUBROUTINE evaluate_dipoles(state_a, state_b, gauge, bs, interaction)
      !
      USE PRECISION,    ONLY: dpk
      USE DATA,         ONLY: read_mx, read_mx_p, read_v_mx, write_mx
      USE utils,        ONLY: print_mx, banded_mat_v_mul, asciifile, dfile
      USE one_e_matrix, ONLY: bsp_overlap
      USE my_types,     ONLY: state_index, basis, atomic_system, init_atomic_system
      !
      IMPLICIT NONE
      !
      !ARG!      
      TYPE(state_index),         INTENT(in)  :: state_a
      TYPE(state_index),         INTENT(in)  :: state_b
      CHARACTER(len=6),          INTENT(in)  :: gauge
      TYPE(basis),               INTENT(in)  :: bs
      CHARACTER(len=6),          INTENT(in)  :: interaction
      !CHARACTER(len=3), OPTIONAL,intent(in)  :: surface_term !(or make it logical?)
      !LOC!
      REAL(DPK), DIMENSION(:,:), POINTER     :: c_a, c_b    ! < B_i  | P_nl(x) >
      REAL(DPK), DIMENSION(:),   POINTER     :: e_a, e_b    ! 
      REAL(dpk), DIMENSION(:,:), POINTER     :: rb          ! <B_i| r^q | B_j > 
      REAL(dpk), DIMENSION(:,:), POINTER     :: bdr         ! <B_i| d/dr | B_j > 
      !
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: e_ab        ! e_ab = e_a -e_b)
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: dipole      ! <P_nl|T_g  | P_nl>
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: dipole_r    ! <P_nl|r | P_nl>
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: dipole_v    ! <P_nl|d/dr | P_nl>
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: S_r         ! == (W_ab + P_ab)/2.
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: S_pr        ! == S_p + e_ab*S_r
      !
      TYPE(atomic_system)                    :: atomo
      REAL(dpk)                              :: k_ab, de_ab ! 1-e angular factor 
      !
      REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: vr          ! 

      !
      INTEGER                                :: na,la,ma
      INTEGER                                :: nb,lb,mb
      INTEGER                                :: ne_a, ne_b
      INTEGER                                :: nb_a, nb_b
      INTEGER                                :: cl
      REAL(dpk)                              :: rmax, znuc
      INTEGER                                :: lab
      !
      REAL(DPK)                              :: rcn
      INTEGER                                :: icn, i, j 

      !
      INTEGER                                :: nascii, nbin
      integer                                :: mode
      !EXE!
      
      !      IF(PRESENT(surface_term)) THEN         
      !      ELSE         
      !ENDIF

      WRITE(*,*) '@ subroutine evaluate_dipoles in.'
      nascii = 1
     
      call init_atomic_system(atomo)
     
      !linearly polarized light 
      na = state_a%n
      la = state_a%l
      ma = state_a%m
      !
      nb = state_b%n
      lb = state_b%l
      mb = state_b%m
      !
      rmax = bs%t(bs%n+1)
      znuc = atomo%znuc

      !

      ! angular momentum factor;  

      
      lab = MAX(la,lb)
      k_ab = DSQRT( REAL( ( lab**2 - mb**2 ), dpk )/ REAL(  (4*lab**2-1), dpk) )

      !...,  la
     
      CALL READ_V_MX( la, e_a, c_a, "h1e-")     ! read en + coeff for la symmetry
      nb_a = SIZE(c_a(1,:))
      ne_a = SIZE(e_a)

      !,,,, lb = la + 1

      CALL READ_V_MX( lb, e_b, c_b, "h1e-")     ! read en + coeff for lb symmetry
      nb_b = SIZE(c_b(1,:))
      ne_b = SIZE(e_b)
     
      !
      WRITE(*,'(a60,i10)')  'basis dimension   nb = ',   bs%n
      WRITE(*,'(a60,i10)')  'basis order       kb = ',   bs%k
      !
      WRITE(*,'(a60,i10)')           'la = ', la
      WRITE(*,'(a60,i10)')           'lb = ', lb
      WRITE(*,'(a60,i10)')          'lab = ', lab
      WRITE(*,'(a60,2i10)') '(ne_a,nb_a) = ', ne_a, nb_a
      WRITE(*,'(a60,2i10)') '(ne_b,nb_a) = ', ne_b, nb_b
      WRITE(*,'(a60,E25.14)')          'Kab = ',  k_ab
      
      
      ALLOCATE(   dipole( ne_a, ne_b) )
      ALLOCATE( dipole_v( ne_a, ne_b) )
      ALLOCATE( dipole_r( ne_a, ne_b) )
      ALLOCATE(      S_r( ne_a, ne_b) )
      ALLOCATE(     S_pr( ne_a, ne_b) )
      ALLOCATE(     e_ab( ne_a, ne_b) )
      ALLOCATE(       vr( ne_b)       )  ! temporary storage vector



     ! d(a,b) = Kab * Tab


      cl  = 1

      ! gfortran causes segmentation fault here. Simply removing seems to be working
      ! ifort no problem
      !      rb  = 0.0_dpk
      !      bdr = 0.0_dpk



      IF(gauge=='l') THEN                    
         CALL READ_MX_P(-1,   rb, "bb-r" )     ! read (upper part) < B_i|r|B_j>
      ELSE IF(gauge=='v')  THEN              
         cl = - 1
         CALL READ_MX_P(-1,   rb, "bb-1r")     ! read (upper part) < B_i|1/r|B_j>
         CALL READ_MX_P(-1,  bdr, "bb-dr")     ! read (upper part) < B_i|d/dr|B_j>
      ELSE IF(gauge=='a') THEN                  
         CALL READ_MX_P(-1, rb, "bb-1r2")      ! read (upper part) < B_i|1/r2|B_j>
      ENDIF
      
     !
     !
     !

      dipole = 0.0_dpk

      !,,,,,, inner contribution 
      !
      !  T(a,b) = <a|T_g|b> = Int_{0^R} dr P_a(r) * T_g * P_b(r)
      !

      
      dipole_r = 0.0_dpk
      perform_v_mat_v_mul_dipole: DO i = 1, ne_a
         DO j = 1, ne_b
            vr = c_b(j,:)
            CALL banded_mat_v_mul(rb, vr, bs%k)         ! C_b <--- b*c_b 
            dipole_r(i,j) =  DOT_PRODUCT(c_a(i,:),vr) 
         ENDDO
      ENDDO perform_v_mat_v_mul_dipole
     
     !,,,,,, for the velocity gauge <a|d/dr|b> is also needed

      dipole_v = 0.0_dpk
      evaluate_deriv_part_of_vel: IF(gauge.EQ.'v') THEN  ! velocity gauge
         vr = 0.0_dpk
         perform_v_mat_v_mul_ddr: DO i = 1, ne_a
            DO j = 1, ne_b
               vr = c_b(j,:)              
               CALL banded_mat_v_mul(bdr, vr, bs%k, cl)        ! C_b <--- b*c_b 
               dipole_v(i,j) =  DOT_PRODUCT(c_a(i,:),vr)   !&
            ENDDO
         ENDDO perform_v_mat_v_mul_ddr
      ENDIF evaluate_deriv_part_of_vel
     


     !,,,,,,,,  construct de_ab matrix


      e_ab     = 0.0_dpk
      evaluate_e_ab: DO j = 1, ne_b
         DO i = 1, ne_a
            e_ab(i,j) = e_a(i) - e_b(j)
         ENDDO
      ENDDO evaluate_e_ab
      
     !,,,,,,,, surface corrections  



      S_r  = 0.0_dpk
      S_pr = 0.0_dpk       
      icn  = 0
      evaluate_surface_terms:IF(gauge.NE.'a') THEN

         !
         ! P'(R) = 0 , P(R) /= 0
         !

         non_zero_amplitudes:IF((nb_a.EQ.bs%n-1).AND.((nb_b.EQ.bs%n-1))) THEN

            rcn     = (2.0_dpk * znuc - (la*(la+1)+lb*(lb+1))/(2.0_dpk*rmax))/rmax
            evaluate_surface_term_01:DO j = 1, ne_b
               DO i = 1, ne_a
                  S_r(i,j) =   c_a(i,nb_a) * c_b(j,nb_b)
               ENDDO
            ENDDO evaluate_surface_term_01
            
            S_r =   0.5_dpk * S_r
            
            DO j = 1, ne_b
               DO i = 1, ne_a
                  S_pr(i,j) = ( e_a(i) + e_b(j)  + rcn ) * S_r(i,j)
               ENDDO
            ENDDO

            
            IF(icn.EQ.0) THEN 
               WRITE(*,'(a60)') ' surface terms for dP(R)= 0 calculated.'
               WRITE(*,'(a60,E15.8)') 'box  factor = ', rcn
            ENDIF
            icn = icn + 1
            

         ENDIF non_zero_amplitudes
         !
         ! P(R) = 0 , P'(R) /= 0
         !
         zero_amplitudes:IF((nb_a.EQ.bs%n-2).OR.((nb_a.EQ.bs%n-2))) THEN
            
            rcn     = (bs%k-1)/(bs%t(bs%n+1)-bs%t(bs%n))     ! (kb-1)/[t(nb+1)-t(nb)]
                        
            evaluate_surface_term_00:DO j = 1, ne_b
               DO i = 1, ne_a
                  S_pr(i,j) =  c_a(i,nb_a) * c_b(j,nb_b)
               ENDDO
            ENDDO evaluate_surface_term_00
            
            s_pr = 0.5_dpk * S_pr * rcn**2

!            CALL print_mx(ne_a, s_pr,'s_pr','e')

            IF(icn.EQ.0) THEN 
               WRITE(*,'(a60)') ' surface terms for P(R)= 0 calculated.'
               WRITE(*,'(a60,E15.8)') 'box  factor = ', rcn
            ENDIF
            
            icn = icn + 1
         ENDIF zero_amplitudes
         
      ENDIF evaluate_surface_terms
      
     !
     ! add the contributions
     !

      IF(gauge=='l') THEN
         dipole = dipole_r !- s_pr/(e_ab*e_ab) 
      ELSE IF(gauge=='v') THEN
         dipole =  dipole_v  + (lb-la) * lab * dipole_r - s_r !+ s_pr/e_ab
      ELSE
         dipole = dipole_r    
      ENDIF
      

      !,,, multiply by the angular factor     
      !
      
      !        dipole = k_ab * dipole
      
        IF(la==0) THEN 
           DO i = 1, SIZE(e_b-bs%k)
              WRITE(10,'(2E15.8)') e_b(i)-e_a(1), dipole(1,i)
           ENDDO
        ENDIF
        

        !,,,,,,,,,,, end of calculations ,,,,,,,,,,,,,!

        !
        ! store in binary
        !


        icn = 0
        nbin = 1


        save_for_tdse_calulations:IF(interaction=='tdse') THEN    !tdse calculations


           CALL dfile(nbin,la,lb,'dp1e-',gauge)

           WRITE(nbin) gauge
           WRITE(nbin) la, lb, SIZE(e_a)-icn, SIZE(e_b) - icn         ! partial waves, dimensions
           WRITE(nbin) (e_a(i), i=1, SIZE(e_a)-icn)                   ! eigenenergies of |a>
           WRITE(nbin) (e_b(i), i=1, SIZE(e_b)-icn)                   ! eigenenergies of |b>
           DO j = 1, SIZE(e_b) - icn
              WRITE(nbin) (dipole(i,j), i=1, SIZE(e_a)-icn)                      ! <a|d|b>
           ENDDO
           CLOSE(nbin)
           
        ELSE IF(interaction=='lopt') THEN      ! lopt calculations           

           IF(gauge=='v')      THEN
              mode  = 0 
           ELSE IF(gauge=='l') THEN
              mode= 1
           ELSE IF(gauge=='a') THEN
              mode = 2
           ENDIF
           
           CALL dfile(nbin,la,lb,'d1e-',gauge)
           WRITE(nbin) mode
           WRITE(nbin) la, lb, SIZE(e_a)-icn, SIZE(e_b) - icn               ! partial waves, dimensions
           WRITE(nbin) ((dipole(i,j),j=1,SIZE(e_b)-icn), i=1, SIZE(e_a)-icn)! <a|d|b>
           CLOSE(nbin)
        ENDIF save_for_tdse_calulations
        


        IF(la.EQ.0) THEN 
           OPEN(nbin,file='out/dmx1e-wf1e.out')
           DO j = 1, SIZE(e_b)-icn
              WRITE(nbin,'(3E25.14,i10)') e_b(j)-e_a(1), ABS(dipole(1,j))**2, dipole(1,j), j
           ENDDO
           CLOSE(nbin)
        ENDIF
        
        !     PRINT*, gauge, SIZE(e_a), SIZE(e_b), la,lb
        !     PRINT*, (e_a, i= 2,5)
        !     PRINT*, (e_b, i= 8,13)
        !     DO i = 5, 10
        !        PRINT*, (dipole(i,j),j=9,14)
        !     ENDDO
      
        !,,,,,,, output section
        
      IF((nb_a.EQ.bs%n-1).OR.((nb_a.EQ.bs%n-1))) THEN        
         WRITE(*,*) '& rmx = ', SUM(ABS(c_a(1:ne_a-icn,nb_a))**2), SUM(ABS(c_b(1:ne_a-icn,nb_b))**2)
      ENDIF

      !
      ! S(1s,np;q) = 16.125(-2)/6.75(-1)/3(0)/1.5(1)/1(2) :: hydrogen
      !

        DO i = -2, 2
           rcn = 0._dpk
           DO j=1,SIZE(e_a) - icn 
              rcn = rcn + dipole(na,j)**2 * e_ab(na,j)**i
           ENDDO
           WRITE(*,'(a40,2E20.10,I5)') '# sum_rule_1e::  s(n,q) = ', &
    &        rcn, SUM(dipole(na,:)**2*e_ab(na,:)**i), i
        ENDDO
        
        CALL asciifile(nascii,la,'dmx1e-',gauge)

        de_ab = 0.0_dpk
        write_dmx_from_state_a_to_states_b:DO j = na+1, SIZE(e_b)  - icn

           !        PRINT*, e_a(n),e_b(j), e_ab(j

           IF(e_ab(j,na).EQ.0.0_dpk) THEN 

              IF(gauge.EQ.'l')  THEN
                 de_ab = 10.0_dpk
              ELSE 
                 de_ab = 0.1_dpk
              ENDIF
           ELSE 
              de_ab = ABS(e_ab(na,j))           
              IF(gauge=='v') de_ab = 1.0_dpk/de_ab
              IF(gauge=='a') de_ab = 1.0_dpk/de_ab**3
           ENDIF
           
           !        IF(e_b(j).GT.0.0_dpk) THEN
           !        WRITE(nascii,'(3E30.15,1X,I5)') -e_ab(na,j),  de_ab*ABS(  dipole(na,j) )**2  &
           !                                              ,de_ab * ABS( dipole_r(na,j) )**2, j
           
           WRITE(nascii,'(3E30.15,1X,2I5)') -e_ab(na,j),  de_ab*ABS(dipole(na,j) )**2, dipole(na,j), j, na
           !     ENDIF
           
        ENDDO write_dmx_from_state_a_to_states_b
        CLOSE(nascii)
    

     WRITE(*,'(a60,i10)') 'n_initial = ', na
     WRITE(*,'(a60,i10)') 'n_final   = ', SIZE(e_a) - icn
     WRITE(*,'(a60,i10)') 'n_cut = ', icn
     


     WRITE(*,'(a5,1x,4a14)') " n\m","1","2","3","4"
     DO  i = 1, 6
        WRITE(*,'(2X,I3,1X,1P8E16.8)') i, ( dipole( i,j ), j = 1, 6)
     ENDDO


     DEALLOCATE(e_a,e_b,c_a,c_b,dipole,dipole_r,dipole_v,vr)
     WRITE(*,*) '@ subroutine evaluate_dipoles out.'

   END SUBROUTINE evaluate_dipoles
    !S
    !S
    !S
    !
    !     Calculates < P_{nl} | T_G | B_j(r) >   n = 1,2,...,nb-2 
    !                                            j = 1,2,...., nb
    !
    !S
    !S
    !S
    SUBROUTINE evaluate_overlap_dipoles(state_a, state_b, gauge, bs)
      !
      USE PRECISION,    ONLY: dpk
      USE DATA,         ONLY: read_mx, read_mx_p, read_v_mx, write_mx
      USE utils,        ONLY: print_mx, banded_mat_v_mul, asciifile, dfile
      USE one_e_matrix, ONLY: bsp_overlap
      USE my_types,     ONLY: state_index, basis, atomic_system, init_atomic_system

      !
      IMPLICIT NONE
      !
      !ARG!      
      TYPE(state_index),         INTENT(in)  :: state_a
      TYPE(state_index),         INTENT(in)  :: state_b
      CHARACTER(len=6),          INTENT(in)  :: gauge
      TYPE(basis),               INTENT(in)  :: bs
      !CHARACTER(len=3), OPTIONAL,intent(in) :: surface_term !(or make it logical?)
      !LOC!
      REAL(DPK), DIMENSION(:,:), POINTER     :: c_a, c_b    ! < B_i  | P_nl(x) >
      REAL(DPK), DIMENSION(:),   POINTER     :: e_a, e_b    ! 
      REAL(dpk), DIMENSION(:,:), POINTER     :: rb          ! <B_i| r^q | B_j > 
      REAL(dpk), DIMENSION(:,:), POINTER     :: bdr         ! <B_i| d/dr | B_j > 
      !
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: e_ab        ! e_ab = e_a -e_b)
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: dipole      ! <P_nl|T_g  | P_nl>
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: dipole_r    ! <P_nl|r | P_nl>
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: dipole_v    ! <P_nl|d/dr | P_nl>
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: S_r         ! == (W_ab + P_ab)/2.
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: S_pr        ! == S_p + e_ab*S_r
      !
      TYPE(atomic_system)                    :: atomo
      REAL(dpk)                              :: k_ab        ! 1-e angular factor 
      !
      REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: vr          ! 
      INTEGER                                :: ndim
      !
      INTEGER                                :: na,la,ma
      INTEGER                                :: nb,lb,mb
      INTEGER                                :: ne_a, ne_b
      INTEGER                                :: nb_a, nb_b
      INTEGER                                :: cl
      REAL(dpk)                              :: rmax, znuc
      INTEGER                                :: lab
      INTEGER                                :: i, j 
      !
      !      REAL(DPK)                              :: rcn
      !      INTEGER                                :: icn, i, j 
      !
      INTEGER                                :: nascii, nbin
      integer                                :: mode

      !EXE!
      
      !      IF(PRESENT(surface_term)) THEN       
      !      ELSE         
      !ENDIF

      WRITE(*,*) '@ subroutine evaluate_overlap_dipoles in.'
      nascii = 1
     
      call init_atomic_system(atomo)
     
      !linearly polarized light 
      na = state_a%n
      la = state_a%l
      ma = state_a%m
      !
      nb = state_b%n
      lb = state_b%l
      mb = state_b%m
      !
      rmax = bs%t(bs%n+1)
      znuc = atomo%znuc

      ndim = bs%n

      ! angular momentum factor;  
      
      lab = MAX(la,lb)
      k_ab = DSQRT( REAL( ( lab**2 - mb**2 ), dpk )/ REAL(  (4*lab**2-1), dpk) )

      !...,  la
     
      CALL READ_V_MX( la, e_a, c_a, "h1e-")     ! read en + coeff for la symmetry
      nb_a = SIZE(c_a(1,:))
      ne_a = SIZE(e_a)

      !,,,, lb = la + 1

      CALL READ_V_MX( lb, e_b, c_b, "h1e-")     ! read en + coeff for lb symmetry
      nb_b = SIZE(c_b(1,:))
      ne_b = SIZE(e_b)

      !
      WRITE(*,'(a60,i10)')  'basis dimension   ndim = ',   ndim
      WRITE(*,'(a60,i10)')  'basis order         kb = ',   bs%k
      !
      WRITE(*,'(a60,i10)')           'la = ', la
      WRITE(*,'(a60,i10)')           'lb = ', lb
      WRITE(*,'(a60,i10)')          'lab = ', lab
      WRITE(*,'(a60,2i10)') '(ne_a,nb_a) = ', ne_a, nb_a
      WRITE(*,'(a60,2i10)') '(ne_b,nb_a) = ', ne_b, nb_b
      WRITE(*,'(a60,E25.14)')          'Kab = ',  k_ab
            
      ALLOCATE(   dipole( ne_a, ndim) )
      ALLOCATE( dipole_v( ne_a, ndim) )
      ALLOCATE( dipole_r( ne_a, ndim) )
      ALLOCATE(      S_r( ne_a, ndim) )
      ALLOCATE(     S_pr( ne_a, ndim) )
      ALLOCATE(     e_ab( ne_a, ndim) )
      ALLOCATE(       vr( ndim)       )  ! temporary storage vector



      ! d(a,b) = Kab * Tab

      ! rb(nb-2,nb-2)

      cl  = 1
      rb  = 0.0_dpk
      bdr = 0.0_dpk

      IF(gauge=='l') THEN

         mode = 1
         CALL READ_MX_P(-1,   rb, "bb-r" )     ! read (upper part) < B_i|r|B_j>

      ELSE IF(gauge=='v')  THEN              

         mode = 0
           cl = -1

         CALL READ_MX_P(-1,   rb, "bb-1r")     ! read (upper part) < B_i|1/r|B_j>
         CALL READ_MX_P(-1,  bdr, "bb-dr")     ! read (upper part) < B_i|d/dr|B_j>

      ELSE IF(gauge=='a') THEN
         mode = 2
         CALL READ_MX_P(-1, rb, "bb-1r2")      ! read (upper part) < B_i|1/r2|B_j>
      ENDIF
      
      !
      ! dipole(nb,nb)
      !

      dipole = 0.0_dpk

      !,,,,,, inner contribution 
      !
      !  T(a,j) = <a(n)|T_g|B_j> = Int_{0^R} dr P_(n_a,l_a)(r) * T_g * B_j(r)
      !

      ! dipole_r     (nb-2,nb)
      ! vr           (nb)
      ! c_a          (nb-2,nb-2)
      ! rb           (nb,kb)


      dipole_r = 0.0_dpk           !(nb-2 x nb)
      perform_v_mat_v_mul_dipole: DO i = 1, ne_a
         vr(1)        = 0.0_dpk                          
         vr(2:ndim-1) = c_a(i,1:nb_a)          ! nb_a = ndim - 2 = nb - 2 
         vr(ndim)     = 0.0_dpk
         CALL banded_mat_v_mul(rb, vr, bs%k)   ! vr(nb) <--- T(nb,k)*vr(nb) 
         dipole_r(i,1:ndim) = vr(1:ndim)
      ENDDO perform_v_mat_v_mul_dipole
     
     !,,,,,, for the velocity gauge <a|d/dr|b_j> is also needed

      evaluate_deriv_part_of_vel: IF(gauge.EQ.'v') THEN  ! velocity gauge

         dipole_v = 0.0_dpk
         perform_v_mat_v_mul_ddr: DO i = 1, ne_a
            vr(1)        = 0.0_dpk                          
            vr(2:ndim-1) = c_a(i,1:nb_a)          ! nb_a = ndim - 2 = nb - 2 
            vr(ndim)     = 0.0_dpk         
            CALL banded_mat_v_mul(bdr, vr, bs%k, cl)   ! vr(nb) <--- T_v(nb,k)*vr(nb) !! 
            dipole_v(i,1:ndim) = vr(1:ndim)
         ENDDO perform_v_mat_v_mul_ddr

      ENDIF evaluate_deriv_part_of_vel
     
      IF((gauge=='l').OR.(gauge=='a')) THEN
         dipole = dipole_r 
      ELSE IF(gauge=='v') THEN
         dipole =  -dipole_v  + (lb-la) * lab * dipole_r !- s_r !+ s_pr/e_ab    ! no 
      ENDIF
      

!!%,,,,,,,,  construct de_ab matrix
!!%      e_ab     = 0.0_dpk
!!%      evaluate_e_ab: DO j = 1, ne_b
!!%         DO i = 1, ne_a
!!%            e_ab(i,j) = e_a(i) - e_b(j)
!!%         ENDDO
!!%      ENDDO evaluate_e_ab
!!%
!!%
!!%      !
!!%      !,,,,,,,, surface corrections
!!%      !
!!%      non_zero_amplitudes:IF((nb_a.EQ.bs%n-1)) THEN
!!%
!!%         ! P'(R) = 0.0 , P(R)/= 0  ! non-zero BC for the wf
!!%         
!!%         S_r  = 0.0_dpk
!!%         S_pr = 0.0_dpk 
!!%
!!%
!!%         rcn     = (2.0_dpk * znuc - (la*(la+1)+lb*(lb+1))/(2.0_dpk*rmax))/rmax
!!%         s_r     = 0.0_dpk
!!%         s_pr    = 0.0_dpk      
!!%         evaluate_surface_terms:IF(gauge.NE.'a') THEN
!!%            
!!%            DO j = 1, ne_b
!!%               DO i = 1, ne_a
!!%                  S_r(i,j) =   0.5_dpk * c_a(i,nb_a) * c_b(j,nb_b)
!!%                  S_pr(i,j) = ( e_a(i) + e_b(j)  + rcn ) * S_r(i,j)
!!%               ENDDO
!!%            ENDDO
!!%            
!!%         ENDIF evaluate_surface_terms
!!%
!!%
!!%         IF(gauge=='l') THEN
!!%            dipole = dipole_r - s_pr/(e_ab*e_ab)
!!%         ELSE IF(gauge=='v') THEN        
!!%            dipole =  dipole_v  + (lb-la) * lab * dipole_r - s_r !+ s_pr/e_ab    ! no            
!!%         ELSE
!!%            dipole = dipole_r    ! + dipole_r (outer) 
!!%         ENDIF
!!%         WRITE(*,'(a60)') ' boundary values are included.'
!!%      ELSE     ! P(R) = 0. ! zero BC for the wf
!!%         

         
      WRITE(*,'(a60)') ' no boundary values are included.'

!!%      ENDIF non_zero_amplitudes
      

      !,,, multiply by the angular factor     
      !
      
      !        dipole = k_ab * dipole

        IF(la==0) THEN 
           DO i = 1, SIZE(e_b-bs%k)
              WRITE(10,'(2E15.8)') e_b(i)-e_a(1), dipole(1,i)
           ENDDO
        ENDIF
        
        WRITE(*,'(a60,2i10)')  ' T(n,i) = <P_nl|T|B_i>: '
        WRITE(*,'(a5,1x,6a14)')" i\n","1","2","3","4","5","6"
        DO  j = 1, 6
           WRITE(*,'(2X,I3,1X,1P8E16.8)') j, ( dipole( i,j ), i = 1, 6)
        ENDDO
        WRITE(*,*)'.................................................'
        DO  j = ndim-5, ndim
           WRITE(*,'(2X,I3,1X,1P8E16.8)') j, ( dipole( i,j ), i = 1, 6)
        ENDDO


        !eof calculations !

        
        !i/o

        nbin = 1
           
        CALL dfile(nbin,la,lb,'overlap-d-',gauge)
        WRITE(nbin) mode
        WRITE(nbin) la, lb, SIZE(dipole,dim=1), SIZE(dipole,dim=2)        ! partial waves, dimensions
        WRITE(nbin) ((dipole(i,j),j=1,  SIZE(dipole,dim=2)  ), i=1,SIZE(dipole,dim=1) ) ! <a|d|b>
        CLOSE(nbin)

        

        IF(la.EQ.0) THEN 
           OPEN(nbin,file='out/rwfb-pb.out')
           DO i = 1, SIZE(e_a)
              WRITE(nbin,'(3E25.14,i10)') e_a(i), ABS(dipole(i,ndim))**2, dipole(i,ndim), i
           ENDDO
           CLOSE(nbin)
           OPEN(nbin,file='out/rwfb-bp.out')
           DO j = 1, ndim
              WRITE(nbin,'(i10,2E25.14,2i5)') j, ABS(dipole(na,j))**2, dipole(na,j),na,la
           ENDDO
           CLOSE(nbin)
        ENDIF
                
!!%      CALL asciifile(nascii,la,'d1e-',gauge)
!!%     de_ab = 0.0_dpk
!!%     write_dmx_from_state_a:DO j = na, SIZE(e_a)  - icn
!!%
!!%!        PRINT*, e_a(n),e_b(j), e_ab(j
!!%
!!%        IF(e_ab(j,na).EQ.0.0_dpk) THEN 
!!%
!!%           IF(gauge.EQ.'l')  THEN
!!%              de_ab = 10.0_dpk
!!%           ELSE 
!!%              de_ab = 0.1_dpk
!!%           ENDIF
!!%        ELSE 
!!%           de_ab = ABS(e_ab(na,j))           
!!%           IF(gauge=='v') de_ab = 1.0_dpk/de_ab
!!%           IF(gauge=='a') de_ab = 1.0_dpk/de_ab**3
!!%        ENDIF
!!%        
!!%
!!%        WRITE(nascii,'(3E30.15,1X,I5)') -e_ab(na,j),  de_ab*ABS(dipole(na,j) )**2, dipole(na,j), j
!!%     ENDDO write_dmx_from_state_a
!!%     CLOSE(nascii)
    



     DEALLOCATE(e_a,e_b,c_a,c_b,dipole,dipole_r,dipole_v,vr)
      WRITE(*,*) '@ subroutine evaluate_dipoles out.'
    END SUBROUTINE evaluate_overlap_dipoles
    !
    !
    !     Calculates <P_{nl} | B_j(r) >   n = 1,2,...,nb-2 
    !                                     j = 1,2,...., nb
    !
    SUBROUTINE evaluate_overlaps(state_a, bs)
      !
      USE PRECISION,    ONLY: dpk
      USE DATA,         ONLY: read_v_mx, read_mx_p 
      USE utils,        ONLY: banded_mat_v_mul, datafile
      USE one_e_matrix, ONLY: bsp_overlap
      USE my_types,     ONLY: state_index, basis,atomic_system, init_atomic_system
      !
      IMPLICIT NONE
      !
      !ARG!      
      TYPE(state_index),         INTENT(in)  :: state_a
      TYPE(basis),               INTENT(in)  :: bs
      !CHARACTER(len=3), OPTIONAL,intent(in)  :: surface_term !(or make it logical?)
      !LOC!
      REAL(DPK), DIMENSION(:,:), POINTER     :: c_a         ! < B_i  | P_nl(x) >
      REAL(DPK), DIMENSION(:),   POINTER     :: e_a         ! 
      REAL(dpk), DIMENSION(:,:), POINTER     :: rb          ! <B_i|B_j > 
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: overlap     ! <P_nl| B_j>
      !
      TYPE(atomic_system)                    :: atomo
      !
      REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: vr          ! 
      INTEGER                                :: ndim
      !
      INTEGER                                :: na,la,ma
      INTEGER                                :: ne_a
      INTEGER                                :: nb_a
      REAL(dpk)                              :: rmax, znuc
      !
      INTEGER                                :: nascii, nbin
      INTEGER                                :: i,j

      !EXE!
      
      !      IF(PRESENT(surface_term)) THEN         
      !      ELSE         
      !ENDIF

      WRITE(*,*) '@ subroutine evaluate_overlaps in.'
      nascii = 1
     
      call init_atomic_system(atomo)
     
      !linearly polarized light 
      na = state_a%n
      la = state_a%l
      ma = state_a%m
      !
      rmax = bs%t(bs%n+1)
      znuc = atomo%znuc
      ndim = bs%n



      !,,,,,, inner contribution 

      !
      !  T(a,j) = <P_{na,la}|B_j> = Int_{0^R} dr P_(n_a,l_a)(r) * B_j(r)
      !

      !dims...
      ! overlap(nb-2,nb),  vr(nb), c_a(nb-2,nb-2), rb(nb,kb)
      !


      !...,  la
     
      CALL READ_V_MX( la, e_a, c_a, "h1e-")     ! read en + coeff for la symmetry
      nb_a = SIZE(c_a(1,:))
      ne_a = SIZE(e_a)


      ALLOCATE(   overlap( ne_a, ndim) )
      ALLOCATE(       vr( ndim)       )          ! temporary storage vector

      rb  = 0.0_dpk
      CALL READ_MX_P(-1, rb, "bb-1" )            ! read (upper part) < B_i|B_j>
      
      overlap = 0.0_dpk           
      perform_v_mat_v_mul_dipole: DO i = 1, ne_a

         vr(1)        = 0.0_dpk                          
         vr(2:ndim-1) = c_a(i,1:nb_a)          ! nb_a = ndim - 2 = nb - 2 
         vr(ndim)     = 0.0_dpk

         CALL banded_mat_v_mul(rb, vr, bs%k)   ! vr(nb) <--- T(nb,k)*vr(nb) 
         overlap(i,1:ndim) = vr(1:ndim)

      ENDDO perform_v_mat_v_mul_dipole       !eof main calculation for overlap

      !i/o!  store in binary

      nbin = 1           
      CALL datafile(nbin,la,'overlap-1-') 
      WRITE(nbin) la, SIZE(overlap,dim=1), SIZE(overlap,dim=2)        ! partial waves, dimensions
      WRITE(nbin) ((overlap(i,j),j=1,  SIZE(overlap,dim=2)  ), i=1,SIZE(overlap,dim=1) ) ! <a|d|b>
      CLOSE(nbin)

      !i/o check!

      
      WRITE(*,'(a60,i10)')  'basis dimension   ndim = ',   ndim
      WRITE(*,'(a60,i10)')  'basis order         kb = ',   bs%k
      WRITE(*,'(a60,i10)')                      'la = ', la
      WRITE(*,'(a60,2i10)')            '(ne_a,nb_a) = ', ne_a, nb_a
     
      WRITE(*,'(a60,2i10)')  ' O(n,i) = <P_nl|B_i>: '
      WRITE(*,'(a5,1x,6a14)')" i\n","1","2","3","4","5","6"
      DO  j = 1, 6
         WRITE(*,'(2X,I3,1X,1P8E16.8)') j, ( overlap( i,j ), i = 1, 6)
      ENDDO
      WRITE(*,*)'.................................................'
      DO  j = ndim-5, ndim
         WRITE(*,'(2X,I3,1X,1P8E16.8)') j, ( overlap( i,j ), i = 1, 6)
      ENDDO

      
      IF(la.EQ.0) THEN 
         OPEN(nbin,file='out/overlap.out')
         DO i = 1, SIZE(e_a)
            WRITE(nbin,'(3E25.14,i10)') e_a(i), ABS(overlap(i,ndim))**2, overlap(i,ndim), i
         ENDDO
         CLOSE(nbin)
      ENDIF

      DEALLOCATE(e_a,c_a,overlap,vr)
      WRITE(*,*) '@ subroutine evaluate_overlaps out.'
      
      END SUBROUTINE evaluate_overlaps
    !
    !


!
! stores  in ascii P_nl(r), dP_nl(r), r
!
!

    SUBROUTINE save_wf1e_nl(n, l, e, p, dp, r, np)
      !  
      USE units, ONLY : dpk 
      !
      IMPLICIT NONE
      !
      INTEGER                  :: N
      INTEGER                  :: L
      REAL(DPK)                :: e
      REAL(DPK), DIMENSION(np) :: p,dp    ! radial functions to be stored
      REAL(DPK), DIMENSION(np) :: r       ! r_i
      INTEGER                  :: np      ! nof grid points
      !
      INTEGER nfile
      INTEGER I
      !      REAL(dpk)               ke, nrm
      CHARACTER(LEN=6)        sn, sl
      CHARACTER(len=30)       filename
      !...............

      nfile = 1

      WRITE(sl,'(i6)') l
      WRITE(sn,'(i6)') n

      filename = "out/wf1e-"//TRIM(ADJUSTL(SN))//TRIM(ADJUSTL(SL))//".out"

      OPEN(nfile, file=filename)
      DO i = 1, np
         WRITE(nfile,'(3E25.14)')  r(i), p(i), dp(i)
      ENDDO      
      CLOSE(nfile)

      WRITE(*,*) '#'
      WRITE(*,'(a45,1X,i3 )') '# save_wf1e_nl::                        l = ', l
      WRITE(*,'(a45,1X,i5 )') '# save_wf1e_nl::                        n = ', n
      WRITE(*,'(a45,1X,E25.14 )') '# save_wf1e_nl::                       ek = ', e
!      IF(e.GE.0) THEN
!         WRITE(*,'(a45,1X,E20.14)') '# save_wf1e_nl::                       ke = ', SQRT(2*e)
!         WRITE(*,'(a45,1X,E20.14)') '# save_wf1e_nl::                       ak = ', 1./SQRT(nrm)
!      ENDIF
      WRITE(*,'(a45,1X,a20)') '# save_wf1e_nl::                 filename = ', filename
      WRITE(*,*) '#'

      RETURN
    END SUBROUTINE SAVE_WF1E_NL
    !
    !    f(x) = Sum_n c_n * P_nl(x),  c_n = < P_nl(x)| f(x) >
    !
    SUBROUTINE WRITE_FX(L, CE, CEF)

      USE param
      USE set_grid
      USE one_e_matrix
      USE utils,       ONLY: asciifile
      !      
      IMPLICIT NONE
      ! arguments
      INTEGER                       :: L
      REAL(DPK), DIMENSION(:,:)     :: CE    ! coefficients
      REAL(DPK), DIMENSION(:)       :: CEF   ! coefficients
      ! locals
      REAL(DPK), DIMENSION(NB)      :: cb
      REAL(DPK), DIMENSION(NB+KX)   :: T     !  r_i, knot sequence
      REAL(DPK), DIMENSION(NP)      :: P     !  P(r_i)
      REAL(DPK)  BVALUE
      INTEGER    DERIV
      INTEGER    FILE
      INTEGER    i, j, jj
      !
      !!
      !
      
      FILE = 1
      
      CALL R_GRID          ! wf grid
      CALL MKGRID(T)       ! b-splines grid
      
      DERIV = 0            !
      
      !.....!

      DO J = 1, NO

         p(j) = 0.0d+00
         sum_eigenstates: DO i = 1, SIZE(ce,1)
            
            cb(1) = 0.0D+00
            DO jj = 1, SIZE( ce(i,:) )
               cb(jj + 1) = ce(i, jj)
            END DO
            
            p(j) = p(j) + cef(i) * bvalue( t, cb, nb, kb, r(j), 0 )

         END DO sum_eigenstates
      ENDDO

      CALL asciifile(FILE, L,"fx-")
      DO J = 1, NO
         WRITE(FILE, *)  R(J), P(J) !, fx_gauss_real(r(j))
      ENDDO
      
      CLOSE(FILE) 
      

    END SUBROUTINE WRITE_FX
    !
    !
    !
    !
    !    f(x) = Sum_n c_n * P_nl(x),  c_n = < P_nl(x)| f(x) >
    !
    !
    ! save real(f(x)), imag(f(x)) and norm(f(x)) of a complex valued wf
    !
    !
    SUBROUTINE WRITE_FX_COMPLEX(l, tm, en, cet, ce)

      USE param
      USE set_grid
      USE one_e_matrix
      USE utils,                 ONLY: asciifile
      USE units,                 ONLY: zi
      !      
      IMPLICIT NONE
      ! arguments
      INTEGER                                :: L
      REAL(dpk)                              :: TM
      REAL(DPK),    DIMENSION(:), INTENT(IN) :: EN   ! en(i) 
      REAL(DPK),  DIMENSION(:,:), INTENT(IN) :: CE   ! b-splines  coeff : P_nl = S_j ce(i) B_j(x)
      COMPLEX(DPK), DIMENSION(:), INTENT(IN) :: CET  ! eigenstate coeff : psi = S_n ce_t *Pnl(x)
      ! locals
      COMPLEX(DPK),       DIMENSION(NP)      :: PIN   !  P(r_i) interaction picture
      COMPLEX(DPK),       DIMENSION(NP)      :: PSC   !  P(r_i) schrodinger picture
      REAL(DPK),          DIMENSION(NB)      :: cb
      REAL(DPK),          DIMENSION(NB+KX)   :: T     !  r_i, knot sequence
      REAL(DPK)                              :: BVALUE
      COMPLEX(dpk)                           :: bx, exp_h0_t
      INTEGER             DERIV
      INTEGER             FILE
      INTEGER             i, j, jj
      !
      !!
      !
      
      FILE = 1
      
      CALL R_GRID          ! wf grid
      CALL MKGRID(T)       ! b-splines grid
      
      DERIV = 0            ! only the function values B(x) are needed
      
      !.....!

      DO J = 1, NO
         pin(j) = 0.0D+00
         psc(j) = 0.0D+00

         sum_eigenstates: DO i = 1, SIZE(ce,1)
            
            cb(1) = 0.0D+00
            DO jj = 1, SIZE( ce(i,:) )
               cb(jj + 1) = ce(i, jj)
            END DO

            IF(tm == 0.0D+00) THEN 
               exp_h0_t = 1.0D+00
            ELSE
               exp_h0_t = EXP(-zi*en(i)*tm)
            ENDIF

            bx = bvalue( t, cb, nb, kb, r(j), 0 ) * cet(i)
            
            psc(j) = psc(j) + bx
            pin(j) = pin(j) + bx * exp_h0_t 
            
         END DO sum_eigenstates
      ENDDO

      CALL asciifile(FILE, l,"fx-bs-")
      DO I = 1, NO
         WRITE(FILE, '(4E25.12)')  r(i), ABS(psc(i))**2, ABS(pin(i))**2
      ENDDO
      
      CLOSE(FILE) 
      
    END SUBROUTINE WRITE_FX_COMPLEX
    !
    !  save radial parts of psi(t) 
    !
    !    f(x) = Sum_n c_n * P_nl(x),  c_n = < P_nl(x)| f(x) >
    !
    !
!!%    SUBROUTINE store_radial_psi_bs(it, ce, cet)
!!%
!!%      USE param
!!%      USE set_grid
!!%      USE one_e_matrix
!!%      USE utils,                 ONLY: asciifile
!!%
!!%      !      
!!%      IMPLICIT NONE
!!%      ! arguments
!!%
!!%      INTEGER                                  :: IT
!!%      REAL(DPK),  DIMENSION(:,:,:), INTENT(IN) :: CE   ! b-splines  coeff : P_nl = S_j ce(i) B_j(x)
!!%      COMPLEX(DPK), DIMENSION(:,:), INTENT(IN) :: CET  ! eigenstate coeff : psi = S_n ce_t *Pnl(x)
!!%      ! locals
!!%      COMPLEX(DPK), ALLOCATABLE, DIMENSION(:,:):: P     !  P(r_i) schrodinger picture
!!%      REAL(DPK),          DIMENSION(nb)        :: cb    !  b-spline coeff
!!%      REAL(DPK),          DIMENSION(nb+kx)     :: T     !  r_i, knot sequence
!!%      REAL(DPK)                                :: bvalue
!!%      INTEGER             DERIV
!!%      INTEGER             FILE
!!%      INTEGER             i, j, jj
!!%      INTEGER             l,lmn,lmx
!!%      INTEGER             n_1st, n_last
!!%      !
!!%      !!
!!%      !
!!%      
!!%      FILE = 1
!!%      
!!%      CALL R_GRID          ! wf grid
!!%      CALL MKGRID(T)       ! b-splines grid
!!%      
!!%      DERIV = 0            ! only the function values B(x) are needed
!!%      
!!%      
!!%      lmn = 0
!!%      lmx = SIZE(cet,2)
!!%
!!%      n_1st  = 1
!!%      n_last = SIZE(cet,1)
!!%
!!%      ALLOCATE( p(no, lmin:lmax) )
!!%
!!%      !.....!
!!%
!!%
!!%      CALL asciifile(file, it, "psi_bs_t")
!!%
!!%      radial_grid: DO j = 1, no
!!%
!!%         get_partial_wave: DO l = lmn, lmx
!!%
!!%            p(j,l) = 0.0D+00
!!%            sum_eigenstates: DO i = n_1st, n_last ! i eigenstates index    
!!%            
!!%               cb(1) = 0.0D+00 
!!%               DO jj = 1, SIZE( ce(i,:,l) )        !jj B-splines index
!!%                  cb(jj + 1) = ce(i, jj, l)
!!%               END DO
!!%            
!!%               p(j,l) = p(j,l) +  cet(i,l)  * bvalue( t, cb, nb, kb, r(j), 0 ) 
!!%            
!!%            END DO sum_eigenstates
!!%         ENDDO get_partial_wave
!!%
!!%         WRITE(FILE, '(7E25.12)')  r(i), ( ABS(p(i,l))**2, l = lmn, lmx)
!!%         
!!%      ENDDO radial_grid
!!%      
!!%      CLOSE(FILE) 
!!%      
!!%    END SUBROUTINE STORE_RADIAL_PSI_BS
!!%    !
    ! gets the boundary value
    !
    SUBROUTINE get_fn_boundary(ix, ce, fn_b)

      USE param
      USE set_grid
      USE one_e_matrix
      USE utils,                 ONLY: asciifile
      !      
      IMPLICIT NONE
      ! arguments
      INTEGER                               :: ix
      REAL(DPK), DIMENSION(:,:), INTENT(IN) :: ce        ! eigen-coefficients
      REAL(DPK), DIMENSION(:),  INTENT(OUT) :: fn_b      ! f(b) boundary value
      ! locals
      REAL(DPK), DIMENSION(NB)              :: cb        ! bsp-coefficients
      REAL(DPK), DIMENSION(NB+KX)           :: t         ! r_i, knot sequence
      REAL(DPK)  BVALUE
      INTEGER    i, j
      !
      !
      !

      
      CALL R_GRID          ! wf grid
      CALL MKGRID(T)       ! physical grid

      
      !      f(ix) = S_n c_n * Pn(ix)
      
      
      fn_b = 0.0D+00
      sum_eigenstates: DO i = 1, SIZE(ce,1)

              
         cb(1) = 0.0D+00
         DO j = 1, SIZE( ce(i,:) )
            cb(j + 1) = ce(i, j)
         END DO
         
         fn_b(i) = bvalue( t, cb, nb, kb, r(ix), 0 )

      END DO sum_eigenstates


      
    END SUBROUTINE GET_FN_BOUNDARY
    !
    !
    !
    SUBROUTINE evaluate_fn_at(ri, t, ce, fn_b)

      USE param
      USE set_grid
      USE one_e_matrix
      USE utils,                 ONLY: asciifile
      !      
      IMPLICIT NONE
      !ARG!
      REAL(DPK)                             :: ri        ! point that the eigenstate value is asked
      REAL(DPK), DIMENSION(:)               :: t         ! r_i, B-splines knot sequence
      REAL(DPK), DIMENSION(:,:), INTENT(IN) :: ce        ! eigen-coefficients
      REAL(DPK), DIMENSION(:),  INTENT(OUT) :: fn_b      ! f(b) boundary value
      !LOC!
      REAL(DPK), DIMENSION(NB)              :: cb        ! bsp-coefficients

      REAL(DPK)  BVALUE
      INTEGER    i, j
      !
      !
      !

      

      !      CALL MKGRID(T)       ! physical grid

      
      !      f(ix) = S_n c_n * Pn(ix)
      
      
      fn_b = 0.0D+00
      sum_eigenstates: DO i = 1, SIZE(ce,1)

              
         cb(1) = 0.0D+00
         DO j = 1, SIZE( ce(i,:) )
            cb(j + 1) = ce(i, j)
         END DO
         
         fn_b(i) = bvalue( t, cb, nb, kb, ri, 0 )

      END DO sum_eigenstates


      
    END SUBROUTINE EVALUATE_FN_AT
!
!
!    Lb(n,m) = p_n(b) * dp_n(b)
!
!
    SUBROUTINE get_fn_boundary_matrix(ce, mx_lb)

      USE param
      USE set_grid
      USE one_e_matrix
      USE utils,                 ONLY: asciifile
      !      
      IMPLICIT NONE
      ! arguments
      REAL(DPK), DIMENSION(:,:), INTENT(IN)  :: ce        ! eigen-coefficients
      REAL(DPK), DIMENSION(:,:), INTENT(OUT) :: mx_lb     ! Lb(n,n'; x=R) boundary value
      ! locals
      REAL(DPK), DIMENSION(NB)               :: cbi,cbj   ! bsp-coefficients
      REAL(DPK), DIMENSION(NB+KX)            :: t         ! r_i, knot sequence
      REAL(DPK)  BVALUE
      INTEGER    ie, je, j
      REAL(DPK) a_s
      !
      !!
      !
      
      CALL R_GRID          ! wf grid
     CALL MKGRID(T)       ! physical grid

!      f(ix) = S_n c_n * Pn(ix)

      a_s = 0.5D+00 * ( KB - 1 ) / ( t(NB + 1 ) - t(NB) )
      
      mx_lb = 0.0D+00
      sum_eigenstates: DO ie = 1, SIZE(ce,1)

         cbi(1) = 0.0D+00
         DO j = 1, SIZE( ce(ie,:) ) ! nb - 1
            cbi(j + 1) = ce(ie, j)
         END DO

         DO je = 1, SIZE(ce,1)

            cbj(1) = 0.0D+00
            DO j = 1, SIZE( ce(je,:) ) ! nb - 1
               cbj(j + 1) = ce(je, j)
            END DO
         
!            mx_lb(ie,je) = bvalue( t, cbi, nb, kb, r(no), 0 )* bvalue( t, cbj, nb, kb, r(no), 1 )

!            if(ie==je)  WRITE(*,*) ie,je,mx_lb(ie,je), a_s * ce(ie,nb-1) * ( ce(je,nb-1) - ce(je,nb-2))
            IF(ie==je) THEN
               WRITE(15,'(2I6,1X,6E25.14)') ie,je, ce(ie,nb-1), ce(je,nb-2),            &
                                                 & bvalue( t, cbi, nb, kb, r(no), 0 ) ,& 
                                                 & bvalue( t, cbj, nb, kb, r(no), 1 ) ,&
                                                 & a_s*( ce(je,nb-1) - ce(je,nb-2))   ,&
                                                 & ce(je,nb-1) 

            ENDIF
         ENDDO
      END DO sum_eigenstates
      
    END SUBROUTINE GET_FN_BOUNDARY_MATRIX


    SUBROUTINE WRITE_FX_BSP(CEF)

     USE param
     USE set_grid
     USE one_e_matrix
     USE functions,   ONLY:fx_gauss_real
     USE utils,       ONLY: asciifile
     !      
     IMPLICIT NONE
     ! arguments
     REAL(DPK), DIMENSION(:)       :: CEF   ! coefficients
     ! locals
     REAL(DPK), DIMENSION(NB)      :: cb
     REAL(DPK), DIMENSION(NB+KX)   :: T     !  r_i, knot sequence
     REAL(DPK), DIMENSION(NP)      :: P     !  P(r_i)
     REAL(DPK)  BVALUE
     INTEGER    DERIV
     INTEGER    FILE
     INTEGER    j, jj
     integer    l
      !
     !!
     !
     
     FILE = 1
     
     CALL R_GRID          ! wf grid
     CALL MKGRID(T)       ! b-splines grid
      
     DERIV = 0            !
     l = 0
     
     !.....!

     cb(1) = 0.0D+00
     DO jj = 1, SIZE( cef)
        cb(jj + 1) = cef(jj)
     END DO
     
     p = 0.0D+00
     DO J = 2, NO
        p(j) = bvalue( t, cb, nb, kb, r(j), 0 )            
      ENDDO
      
      CALL asciifile(FILE, L,"fx-bsp-")
      DO J = 1, NO
         WRITE(FILE, *)  R(J), P(J), fx_gauss_real(r(j))
      ENDDO
      
      CLOSE(FILE) 
      
      
    END SUBROUTINE WRITE_FX_BSP
    !
    ! propability density in the interaction and Schrodinger picture
    !
    !  Psi(I) = exp(ih_0*t) Psi(S)
    !
    SUBROUTINE WRITE_FX_COMPLEX_LEAPFROG(IT, TM, EN, CE, CE_T)
      !
      USE param
      USE set_grid
      USE one_e_matrix
      USE utils,        ONLY: asciifile
      USE units,        ONLY: zi
      !      
      IMPLICIT NONE
      ! arguments
      INTEGER                                   :: IT
      REAL(DPK)                                 :: TM
      REAL(DPK),    DIMENSION (:),   INTENT(IN) :: EN      ! e(i)
      REAL(DPK),    DIMENSION (:,:), INTENT(IN) :: CE      ! P_nl(r) = S_j ce(i,j)*B_j(r)
      COMPLEX(DPK), DIMENSION (:,:), INTENT(IN) :: CE_T    ! psi(r,t) = S_n ce_t(n)*P_nl(r)*exp(-i*en(i)*t) 
      ! locals
      REAL(DPK),    DIMENSION(NB)               :: cb      
      REAL(DPK),    DIMENSION(NB+KX)            :: T                 !  r_i, knot sequence
      COMPLEX(DPK), DIMENSION(NP)               :: fx2_sc, fx3_sc  !  P(r_i) schrodinger picture
      COMPLEX(DPK), DIMENSION(NP)               :: fx2_in, fx3_in  !  P(r_i) interaction picture
      REAL(DPK)                                 :: bvalue, bx
      COMPLEX(DPK)                              :: bx2,bx3
      COMPLEX(DPK)                              :: exp_h0_t 
      REAL(DPK)                                 :: prob_s, prob_i 
      INTEGER                                   :: DERIV
      INTEGER                                   :: FILE
      INTEGER                                   :: i, j, jj
      !
      !!
      !
     
      
      FILE = 1
      
      CALL R_GRID          ! wf grid
      CALL MKGRID(T)       ! b-splines grid
      
      DERIV = 0            ! only the function values B(x) are needed
      
      !.....!

      DO j = 1, NO

         fx2_sc(j) = 0.0D+00
         fx3_sc(j) = 0.0D+00
         fx2_in(j) = 0.0D+00
         fx3_in(j) = 0.0D+00
         sum_eigenstates: DO i = 1, SIZE(ce,1)
            
            cb(1) = 0.0D+00
            DO jj = 1, SIZE( ce(i,:) )
               cb(jj + 1) = ce(i, jj)
            END DO

            bx = bvalue( t, cb, nb, kb, r(j), 0 )      ! P_nl(x)

            IF(tm == 0.0D+00) THEN 
               exp_h0_t = 1.0D+00
            ELSE
               exp_h0_t = EXP(-zi*en(i)*tm)
            ENDIF
            
            bx2 = ce_t(i,2) *  bx
            bx3 = ce_t(i,3) *  bx
            
            fx2_in(j) = fx2_in(j) + bx2
            fx3_in(j) = fx3_in(j) + bx3 

            fx2_sc(j) = fx2_sc(j) + bx2 * exp_h0_t
            fx3_sc(j) = fx3_sc(j) + bx3 * exp_h0_t

         END DO sum_eigenstates
      ENDDO
      

      CALL asciifile(FILE, IT,"fx-bs-")
      DO J = 1, NO
         prob_s = ( CONJG(fx3_sc(j))*fx2_sc(j) + CONJG(fx2_sc(j))*fx3_sc(j))/2.0_dpk
         prob_i = ( CONJG(fx3_in(j))*fx2_in(j) + CONJG(fx2_in(j))*fx3_in(j))/2.0_dpk
         WRITE(FILE, '(4E25.12)')  r(j), prob_i, prob_s
      ENDDO
      
      CLOSE(FILE) 
      
    END SUBROUTINE WRITE_FX_COMPLEX_LEAPFROG
    !
    !
    !
    !S
    !S evaluates < P_(na,la) | K_ab * T_G(a,b) | P_(nb,kb) >
    !S
    !S  na = 1, ..., na_max           G = l,v,a
    !S  nb = 1,...., nb_max
    !S
    !S corrections are added  
    !S

    SUBROUTINE evaluate_dipoles_diatomics(state_a, state_b, gauge, bs, interaction)
      !
      USE PRECISION,    ONLY: dpk
      USE DATA,         ONLY: read_mx, read_mx_p, read_v_mx, write_mx
      USE utils,        ONLY: print_mx, banded_mat_v_mul, asciifile, dfile
      USE one_e_matrix, ONLY: bsp_overlap
      USE my_types,     ONLY: state_index, basis, molecular_system, init_molecular_system

      !
      IMPLICIT NONE
      !
      !ARG!      
      TYPE(state_index),         INTENT(in)  :: state_a
      TYPE(state_index),         INTENT(in)  :: state_b
      CHARACTER(len=6),          INTENT(in)  :: gauge
      TYPE(basis),               INTENT(in)  :: bs
      CHARACTER(len=6),          INTENT(in)  :: interaction
      !CHARACTER(len=3), OPTIONAL,intent(in)  :: surface_term !(or make it logical?)
      !LOC!
      REAL(DPK), DIMENSION(:,:), POINTER     :: c_a, c_b    ! < B_i  | P_nl(x) >
      REAL(DPK), DIMENSION(:),   POINTER     :: e_a, e_b    ! 
      REAL(dpk), DIMENSION(:,:), POINTER     :: rb          ! <B_i| r^q | B_j > 
      REAL(dpk), DIMENSION(:,:), POINTER     :: bdr         ! <B_i| d/dr | B_j > 
      !
         
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: e_ab        ! e_ab = e_a -e_b)
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: dipole      ! <P_nl|T_g  | P_nl>
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: S_r         ! == (W_ab + P_ab)/2.
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:) :: S_pr        ! == S_p + e_ab*S_r
      !
      TYPE(molecular_system)                 :: diatomic
      REAL(dpk)                              :: de_ab ! 1-e angular factor 
      REAL(dpk), ALLOCATABLE, dimension(:)   :: k_ab_l 
      REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: c_al, c_bl
      INTEGER                                :: i_bs, j_bs
      REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: vr0, vr1   ! 
      REAL(dpk)                              :: vr
      !
      INTEGER                                :: ndim
      INTEGER                                :: na, la, ma, pa
      INTEGER                                :: nb, lb, mb, pb
      INTEGER                                :: ne_a, ne_b
      INTEGER                                :: nb_a, nb_b
      INTEGER                                :: nl
      INTEGER                                :: cl
      REAL(dpk)                              :: rmax, znuc
      INTEGER                                :: lab, l
      !
      REAL(DPK)                              :: rcn
      INTEGER                                :: icn, i, j 
      !
      INTEGER                                :: nascii, nbin
      integer                                :: mode
      !EXE!
      
      !      IF(PRESENT(surface_term)) THEN         
      !      ELSE         
      !ENDIF

      WRITE(*,*) '@ subroutine evaluate_diatomic_dipoles in.'
      nascii = 1
     
      call init_molecular_system(diatomic)
     
      !linearly polarized light 
      na  = state_a%n
      la  = state_a%l
      ma  = state_a%m
      pa  = state_a%p
      nl = state_a%n_l
      !
      nb  = state_b%n
      lb  = state_b%l
      mb  = state_b%m
      pb  = state_b%p

      IF((state_a%n_l).NE.(state_b%n_l)) THEN
         PRINT*, "@  evaluate_dipoles_diatomics::"
         PRINT*, "@  state a and state span couple different partial wave basis set"
         PRINT*, "@                                           nl_a = ", state_a%n_l
         PRINT*, "@                                           nl_b = ", state_b%n_l
         stop
      ELSE
         nl = state_a%n_l
      ENDIF


      !

      rmax = bs%t(bs%n+1)
      znuc = diatomic%znuc(1)  ! implemented for homocharged diatomic


      !

      ! angular momentum factor;  

      
      ! l = 1 --> < 0,m_a | Y_10 | 1,0>       = k_ab_l(1)
      ! l = 2 --> < 1,m_a | Y_10 | 2,0>       = k_ab_l(2)
      ! l = 3 --> < 2,m_a | Y_10 | 3,0>       = k_ab_l(3)
      !....
      !
      ! l = nl --> < nl-1,m_a | Y_10 | n,0>   = k_ab_l(nl)

      ALLOCATE(k_ab_l(2*nl-1))

      DO l = 1, 2*nl-1
         k_ab_l(l) = DSQRT( REAL( (l**2 - mb**2 ), dpk )/ REAL((4*l**2-1), dpk) ) 
      ENDDO

      !...,  la
     
      CALL READ_V_MX( la, e_a, c_a, "h1e-")     ! read en + coeff for la symmetry
      nb_a = SIZE(c_a(1,:))                     ! info for sizes in h1e-la.dat
      ne_a = SIZE(e_a)

      !,,,, lb = la + 1

      CALL READ_V_MX( lb, e_b, c_b, "h1e-")     ! read en + coeff for lb symmetry
      nb_b = SIZE(c_b(1,:))                     ! info for sizes in h1e-lb.dat 
      ne_b = SIZE(e_b)
     
      IF(ne_a.NE.ne_b) THEN
         PRINT*, "@  evaluate_dipoles_diatomics::"
         PRINT*, "@  state a and state spans different b-splines set"
         PRINT*, "@                                          nb_a = ", nb_a
         PRINT*, "@                                          nb_a = ", nb_b
      ENDIF






      !
      WRITE(*,'(a60,i10)')  'basis dimension   nb = ',   bs%n
      WRITE(*,'(a60,i10)')  'basis order       kb = ',   bs%k
      !
      WRITE(*,'(a60,i10)')           'la = ', la
      WRITE(*,'(a60,i10)')           'lb = ', lb
      WRITE(*,'(a60,i10)')          'lab = ', lab
      WRITE(*,'(a60,2i10)') '(ne_a,nb_a) = ', ne_a, nb_a
      WRITE(*,'(a60,2i10)') '(ne_b,nb_a) = ', ne_b, nb_b
      WRITE(*,'(a60,20G10.5)')   'Kab_l = ',  k_ab_l
      
      
      ALLOCATE(   dipole( ne_a, ne_b) )
      ALLOCATE(      S_r( ne_a, ne_b) )
      ALLOCATE(     S_pr( ne_a, ne_b) )
      ALLOCATE(     e_ab( ne_a, ne_b) )
      !

      !
      ! d(a,b) = Sum_{l=1,..,nl} [ Kab(l,l+1) * T(a,l/l+1,b) + Kab(l,l-1) * T(a,l/l-1,b)]
      !


      cl  = 1

      ! gfortran causes segmentation fault here. Simply removing seems to be working
      ! ifort no problem
      !      rb  = 0.0_dpk
      !      bdr = 0.0_dpk
      


      IF(gauge=='l') THEN                    
         CALL READ_MX_P(-1,   rb, "bb-r" )     ! read (upper part) < B_i|r|B_j>
      ELSE IF(gauge=='v')  THEN              
         cl = - 1
         CALL READ_MX_P(-1,   rb, "bb-1r")     ! read (upper part) < B_i|1/r|B_j>
         CALL READ_MX_P(-1,  bdr, "bb-dr")     ! read (upper part) < B_i|d/dr|B_j>
      ELSE IF(gauge=='a') THEN                  
         CALL READ_MX_P(-1, rb, "bb-1r2")      ! read (upper part) < B_i|1/r2|B_j>
      ENDIF
      

      ndim = SIZE(rb(1,:))                     ! info for sizes in h1e-lb.dat 


      
     !
     !
     !

      PRINT*, "@ evaluate_dipoles_diatomics::                         ndim = ", ndim

      !,,,,,, inner contribution 
      !
      !  T(a,b) = <a|T_g|b> = Int_{0^R} dr P_a(r) * T_g * P_b(r)
      !


 
      ALLOCATE( c_al( nl, ndim))
      ALLOCATE( c_bl( nl, ndim))
      ALLOCATE(       vr0(ndim))
      ALLOCATE(       vr1(ndim))
      !
      vr0      = 0.0_dpk
      vr1      = 0.0_dpk
      c_al     = 0.0_dpk
      c_bl     = 0.0_dpk
      dipole   = 0.0_dpk
      !
      dipole_initial_state_a: DO i = 1, ne_a
         !
         DO i_bs = 1, ndim 
            DO l = 1, nl
               c_al(l, i_bs) = c_a(i, i_bs + (l-1)*ndim)       ! c_il(n_a)
            END DO
         ENDDO

!         PRINT*, "ia = ", i

         dipole_final_state_b:DO j = 1, ne_b

            DO j_bs = 1, ndim 
               loop_over_partial_waves: DO l = 1, nl
                  c_bl(l, j_bs) = c_b(j, j_bs + (l-1)*ndim)    ! c_jl(n_b)
               ENDDO loop_over_partial_waves
            END DO


!            PRINT*, "ib = ", i
            !
            ! <c_al(l,:)|r|c_bl(l,  :> =  <0|r|1>, <2|r|3>, <4|r|5>, .... <a;l_max|r|b;l_max+1> 
            ! <c_bl(l,: |r|c_al(l+1,:> =  <1|r|2>, <3|r|4>, <5|r|6>, .... <b;l_max-1|r|a;l_max> 
            !

            
            vr0 = c_bl(1,:)
            CALL banded_mat_v_mul(rb, vr0, bs%k )                ! c_b <--- r_b *c_b             
            vr =  DOT_PRODUCT( c_al(1,:), vr0 ) * k_ab_l(1)      ! <0| r * Y_10 |1>  

            sum_over_partial_waves:DO l = 2, nl                  !

               vr0 = c_bl(l,  :)
               vr1 = c_bl(l-1,:)
               CALL banded_mat_v_mul(rb, vr0, bs%k )    ! vr0 <--- rb*vr0   
               CALL banded_mat_v_mul(rb, vr1, bs%k )    ! vr1 <--- rb*vr1 

               
               vr =  vr + DOT_PRODUCT( c_al(l,:), vr0 ) * k_ab_l(2*l-1)  &
                    &   + DOT_PRODUCT( c_al(l,:), vr1 ) * k_ab_l(l)       !<(2*l)|r*Y_l0|(2*l-1)>

            ENDDO sum_over_partial_waves

            dipole(i,j) = vr


         ENDDO dipole_final_state_b
      ENDDO dipole_initial_state_a
     

      DEALLOCATE(c_al,c_bl)

        !,,,,,,,,,,, end of calculations ,,,,,,,,,,,,,!

        !
        ! store in binary
        !


        icn = 0
        nbin = 1


        save_for_tdse_calulations:IF(interaction=='tdse') THEN    !tdse calculations


           CALL dfile(nbin,la,lb,'dp1e-',gauge)

           WRITE(nbin) gauge
           WRITE(nbin) la, lb, SIZE(e_a)-icn, SIZE(e_b) - icn         ! partial waves, dimensions
           WRITE(nbin) (e_a(i), i=1, SIZE(e_a)-icn)                   ! eigenenergies of |a>
           WRITE(nbin) (e_b(i), i=1, SIZE(e_b)-icn)                   ! eigenenergies of |b>
           DO j = 1, SIZE(e_b) - icn
              WRITE(nbin) (dipole(i,j), i=1, SIZE(e_a)-icn)                      ! <a|d|b>
           ENDDO
           CLOSE(nbin)
           
        ELSE IF(interaction=='lopt') THEN      ! lopt calculations           

           IF(gauge=='v')      THEN
              mode  = 0 
           ELSE IF(gauge=='l') THEN
              mode= 1
           ELSE IF(gauge=='a') THEN
              mode = 2
           ENDIF
           
           CALL dfile(nbin,la,lb,'d1e-',gauge)
           WRITE(nbin) mode
           WRITE(nbin) la, lb, SIZE(e_a)-icn, SIZE(e_b) - icn               ! partial waves, dimensions
           WRITE(nbin) ((dipole(i,j),j=1,SIZE(e_b)-icn), i=1, SIZE(e_a)-icn)! <a|d|b>
           CLOSE(nbin)
        ENDIF save_for_tdse_calulations
        


        IF(la.EQ.0) THEN 
           OPEN(nbin,file='out/dmx1e-wf1e.out')
           DO j = 1, SIZE(e_b)-icn
              WRITE(nbin,'(3E25.14,i10)') e_b(j)-e_a(1), ABS(dipole(1,j))**2, dipole(1,j), j
           ENDDO
           CLOSE(nbin)
        ENDIF
        
        !     PRINT*, gauge, SIZE(e_a), SIZE(e_b), la,lb
        !     PRINT*, (e_a, i= 2,5)
        !     PRINT*, (e_b, i= 8,13)
        !     DO i = 5, 10
        !        PRINT*, (dipole(i,j),j=9,14)
        !     ENDDO
      
        !,,,,,,, output section
        
      IF((nb_a.EQ.bs%n-1).OR.((nb_a.EQ.bs%n-1))) THEN        
         WRITE(*,*) '& rmx = ', SUM(ABS(c_a(1:ne_a-icn,nb_a))**2), SUM(ABS(c_b(1:ne_a-icn,nb_b))**2)
      ENDIF

      !
      ! S(1s,np;q) = 16.125(-2)/6.75(-1)/3(0)/1.5(1)/1(2) :: hydrogen
      !

        DO i = -2, 2
           rcn = 0._dpk
           DO j=1,SIZE(e_a) - icn 
              rcn = rcn + dipole(na,j)**2 * e_ab(na,j)**i
           ENDDO
           WRITE(*,'(a40,2E20.10,I5)') '# sum_rule_1e::  s(n,q) = ', &
    &        rcn, SUM(dipole(na,:)**2*e_ab(na,:)**i), i
        ENDDO
        
        CALL asciifile(nascii,la,'dmx1e-',gauge)

        de_ab = 0.0_dpk
        write_dmx_from_state_a_to_states_b:DO j = na+1, SIZE(e_b)  - icn

           !        PRINT*, e_a(n),e_b(j), e_ab(j

           IF(e_ab(j,na).EQ.0.0_dpk) THEN 
              IF(gauge.EQ.'l')  THEN
                 de_ab = 10.0_dpk
              ELSE 
                 de_ab = 0.1_dpk
              ENDIF
           ELSE 
              de_ab = ABS(e_ab(na,j))           
              IF(gauge=='v') de_ab = 1.0_dpk/de_ab
              IF(gauge=='a') de_ab = 1.0_dpk/de_ab**3
           ENDIF
           
           !        IF(e_b(j).GT.0.0_dpk) THEN
           !        WRITE(nascii,'(3E30.15,1X,I5)') -e_ab(na,j),  de_ab*ABS(  dipole(na,j) )**2  &
           !                                              ,de_ab * ABS( dipole_r(na,j) )**2, j
           
           WRITE(nascii,'(3E30.15,1X,2I5)') -e_ab(na,j),  de_ab*ABS(dipole(na,j) )**2, dipole(na,j), j, na
           !     ENDIF
           
        ENDDO write_dmx_from_state_a_to_states_b
        CLOSE(nascii)
    

     WRITE(*,'(a60,i10)') 'n_initial = ', na
     WRITE(*,'(a60,i10)') 'n_final   = ', SIZE(e_a) - icn
     WRITE(*,'(a60,i10)') 'n_cut = ', icn
     


     WRITE(*,'(a5,1x,4a14)') " n\m","1","2","3","4"
     DO  i = 1, 6
        WRITE(*,'(2X,I3,1X,1P8E16.8)') i, ( dipole( i,j ), j = 1, 6)
     ENDDO


     DEALLOCATE(e_a,e_b,c_a,c_b)
     DEALLOCATE(dipole)

     WRITE(*,*) '@ subroutine evaluate_dipoles out.'

   END SUBROUTINE evaluate_dipoles_diatomics
   !
   !
   ! compact evaluation
   !
   !
    SUBROUTINE evaluate_dipoles_diatomics_2(state_a, state_b, gauge, bs, interaction)
      !
      USE PRECISION,    ONLY: dpk
      USE DATA,         ONLY: read_mx, read_mx_p, read_v_mx, write_mx
      USE utils,        ONLY: print_mx, banded_mat_v_mul, asciifile, dfile
      USE one_e_matrix, ONLY: bsp_overlap
      USE my_types,     ONLY: state_index, basis, molecular_system, init_molecular_system

      !
      IMPLICIT NONE
      !
      !ARG!      
      TYPE(state_index),         INTENT(in)    :: state_a
      TYPE(state_index),         INTENT(in)    :: state_b
      CHARACTER(len=6),          INTENT(in)    :: gauge
      TYPE(basis),               INTENT(in)    :: bs
      CHARACTER(len=6),          INTENT(in)    :: interaction
      !LOC!
      REAL(DPK), DIMENSION(:,:), POINTER       :: c_a, c_b    ! < B_i  | P_nl(x) >
      REAL(DPK), DIMENSION(:),   POINTER       :: e_a, e_b    ! 
      REAL(dpk), DIMENSION(:,:), POINTER       :: rb          ! <B_i| r^q | B_j > 
      !
         
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:)   :: e_ab        ! e_ab = e_a -e_b)
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:)   :: dipole      ! <P_nl|T_g  | P_nl>
      !
      TYPE(molecular_system)                   :: diatomic
      REAL(dpk)                                :: de_ab ! 1-e angular factor 
      REAL(dpk), ALLOCATABLE, dimension(:)     :: k_ab_l 
      REAL(dpk), ALLOCATABLE, DIMENSION(:,:)   :: c_l
      REAL(DPK), ALLOCATABLE, DIMENSION(:,:,:) :: dipole_l      ! <P_nl |T_g  | P_nl+1>
      INTEGER                                  :: i_bs, j_bs
      !
      INTEGER                                  :: ndim, nl
      INTEGER                                  :: na, la, ma, pa
      INTEGER                                  :: nb, lb, mb, pb
      INTEGER                                  :: ne_a, ne_b
      INTEGER                                  :: nb_a, nb_b
      REAL(dpk)                                :: rmax, znuc
      INTEGER                                  :: lab, l
      INTEGER                                  :: i, j 
      !
      INTEGER                                  :: nbin
      integer                                  :: mode
      !EXE!
      
      !      IF(PRESENT(surface_term)) THEN         
      !      ELSE         
      !ENDIF

      WRITE(*,*) '@ subroutine evaluate_diatomic_dipoles in.'
     
      call init_molecular_system(diatomic)
     
      !linearly polarized light 
      na  = state_a%n
      la  = state_a%l
      ma  = state_a%m
      pa  = state_a%p
      nl = state_a%n_l
      !
      nb  = state_b%n
      lb  = state_b%l
      mb  = state_b%m
      pb  = state_b%p

      IF((state_a%n_l).NE.(state_b%n_l)) THEN
         PRINT*, "@  evaluate_dipoles_diatomics::"
         PRINT*, "@  state a and state span couple different partial wave basis set"
         PRINT*, "@                                           nl_a = ", state_a%n_l
         PRINT*, "@                                           nl_b = ", state_b%n_l
         stop
      ELSE
         nl = state_a%n_l
      ENDIF


      !

      rmax = bs%t(bs%n+1)
      znuc = diatomic%znuc(1)  ! implemented for homocharged diatomic


      !

      ! angular momentum factor;  

      
      ! l = 1 --> < 0,m_a | Y_10 | 1,0>       = k_ab_l(1)
      ! l = 2 --> < 1,m_a | Y_10 | 2,0>       = k_ab_l(2)
      ! l = 3 --> < 2,m_a | Y_10 | 3,0>       = k_ab_l(3)
      !....
      !
      ! l = nl --> < nl-1,m_a | Y_10 | n,0>   = k_ab_l(nl)


      !...,  la
     
      CALL READ_V_MX( la, e_a, c_a, "h1e-")     ! read en + coeff for la symmetry
      nb_a = SIZE(c_a(1,:))                     ! info for sizes in h1e-la.dat
      ne_a = SIZE(e_a)

      !,,,, lb = la + 1

      CALL READ_V_MX( lb, e_b, c_b, "h1e-")     ! read en + coeff for lb symmetry
      nb_b = SIZE(c_b(1,:))                     ! info for sizes in h1e-lb.dat 
      ne_b = SIZE(e_b)
     
      IF(ne_a.NE.ne_b) THEN
         PRINT*, "@  evaluate_dipoles_diatomics_2::"
         PRINT*, "@  state a and state spans different b-splines set"
         PRINT*, "@                                          nb_a = ", nb_a
         PRINT*, "@                                          nb_b = ", nb_b
      ENDIF
      WRITE(*,'(a60,i10)')  '@  basis dimension   nb = ', bs%n
      WRITE(*,'(a60,i10)')  '@  basis order       kb = ', bs%k
      WRITE(*,'(a60,i10)')                    '@  la = ', la
      WRITE(*,'(a60,i10)')                    '@  lb = ', lb
      WRITE(*,'(a60,i10)')                    '@ lab = ', lab
      WRITE(*,'(a60,2i10)')          '@  (ne_a,nb_a) = ', ne_a, nb_a
      WRITE(*,'(a60,2i10)')          '@  (ne_b,nb_a) = ', ne_b, nb_b
      WRITE(*,'(a60,20G10.5)')             '@  Kab_l = ',  k_ab_l





      ALLOCATE(   dipole( ne_a, ne_b) )
      ALLOCATE(     e_ab( ne_a, ne_b) )

      ALLOCATE( dipole_l( 0:2*nl-1, ne_a, ne_b) )      
      ALLOCATE(   k_ab_l( 0:2*nl-1)             )



      !
      ! d(a,b) = Sum_{l=1,..,nl} [ Kab(l,l+1) * T(a,l/l+1,b) + Kab(l,l-1) * T(a,l/l-1,b)]
      !

      !
      ! angular matrix         Kab(l,l+-1,m_l) = <l m_l| Y_10 | l+-1, m_l> 
      !




      DO l = 0, 2*nl-1
         k_ab_l(l) = DSQRT( REAL( (l**2 - mb**2 ), dpk )/ REAL((4*l**2-1), dpk) ) 
      ENDDO



      PRINT*, "@ evaluate_dipoles_diatomics_2::      angular part set  l_max = ", 2*nl+1 
      !
      DO l = 1, 2*nl+1
         WRITE(*,'(a6,i3,a4, e15.4))') 'K_ab_l(',l,') = ',k_ab_l(l)
      ENDDO

      !
      !             r_ij = <B_i|r|B_j>
      !


      CALL READ_MX_P(-1,   rb, "bb-r" )     ! read upper part 


      ndim = SIZE(rb(1,:))                  ! info for sizes in h1e-lb.dat


      PRINT*, "@ evaluate_dipoles_diatomics_2::       B-splines matrix set  ndim = ", ndim
 



      ALLOCATE(     c_l(  ndim,       0:2*nl) )      




      !,,,,,, inner contribution 
      !
      !  T(a,b) = <a|T_g|b> = Int_{0^R} dr P_a(r) * T_g * P_b(r)
      !
      

 
      dipole_l = 0.0_dpk
      dipole   = 0.0_dpk
      !

      dipole_final_state_b:DO j = 1, ne_b


         assign_odd_coeff_from_b: DO l = 0, nl-1
            bsp_basis_jb:DO j_bs = 1, ndim 

               c_l(j_bs, 2*l+1) = c_b(j, j_bs +  l * ndim )        ! cl_{j(2l+1)}(n_b)
 
           END DO bsp_basis_jb
         ENDDO assign_odd_coeff_from_b
            

         dipole_initial_state_a: DO i = 1, ne_a


               assign_even_coeff_a:DO l = 0, nl-1               
                  bsp_basis_ia:DO i_bs = 1, ndim                      

                     c_l(i_bs, 2*l) = c_a(i, i_bs + l * ndim )       ! cl_{il}(n_a)

                  END DO bsp_basis_ia
               ENDDO assign_even_coeff_a
               

!            PRINT*, "ib = ", i
            !
            ! <c_al(l,:)|r|c_bl(l,  :> =  <0|r|1>, <2|r|3>, <4|r|5>, .... <a;l_max|r|b;l_max+1> 
            ! <c_bl(l,: |r|c_al(l+1,:> =  <1|r|2>, <3|r|4>, <5|r|6>, .... <b;l_max-1|r|a;l_max> 
            !

               
               CALL banded_mat_v_mul(rb, c_l(:,1), bs%k )                           ! cl <--- r_b *cl 
               dipole_l(0, i, j ) = DOT_PRODUCT( c_l(:,0),  c_l(:,1) ) * k_ab_l(1)  ! <0| r * Y_10 |1>  


               all_partial_waves:DO l = 1, nl-1     !
                  CALL banded_mat_v_mul(rb, c_l(:,2*l-1), bs%k )                                   ! cl <--- r_b *cl 
                  CALL banded_mat_v_mul(rb, c_l(:,2*l+1), bs%k )                                   ! cl <--- r_b *cl 
                  dipole_l(2*l,   i, j ) = DOT_PRODUCT( c_l(:,2*l),  c_l(:,2*l+1) ) * k_ab_l(2*l+1)  ! <0| r * Y_10 |1>  
                  dipole_l(2*l-1, i, j ) = DOT_PRODUCT( c_l(:,2*l),  c_l(:,2*l-1) ) * k_ab_l(2*l)    ! <0| r * Y_10 |1>  
               ENDDO all_partial_waves
               dipole(i,j) = SUM(dipole_l(:,i,j))            


               
            ENDDO dipole_initial_state_a
         ENDDO dipole_final_state_b
     

         DEALLOCATE(c_l)

        !,,,,,,,,,,, end of calculations ,,,,,,,,,,,,,!

        !
        ! store in binary
        !



        nbin = 1


        save_for_tdse_calulations:IF(interaction=='tdse') THEN    !tdse calculations

           CALL dfile(nbin,la,lb,'dp1e-',gauge)

           WRITE(nbin) gauge
           WRITE(nbin) la, lb, SIZE(e_a), SIZE(e_b)          ! partial waves, dimensions
           WRITE(nbin) (e_a(i), i=1, SIZE(e_a))                   ! eigenenergies of |a>
           WRITE(nbin) (e_b(i), i=1, SIZE(e_b))                   ! eigenenergies of |b>
           DO j = 1, SIZE(e_b)
              WRITE(nbin) (dipole(i,j), i=1, SIZE(e_a))                      ! <a|d|b>
           ENDDO
           CLOSE(nbin)
           
        ELSE IF(interaction=='lopt') THEN      ! lopt calculations           

           mode = 1
           CALL dfile(nbin,la,lb,'d1e-',gauge)
           WRITE(nbin) mode
           WRITE(nbin) la, lb, SIZE(e_a), SIZE(e_b)               ! partial waves, dimensions
           WRITE(nbin) ((dipole(i,j),j=1,SIZE(e_b)), i=1, SIZE(e_a))! <a|d|b>
           CLOSE(nbin)
        ENDIF save_for_tdse_calulations
        

        IF(la.EQ.0) THEN 
           OPEN(nbin,file='out/dmx1e.out')
           DO j = 1, SIZE(e_b)
              WRITE(nbin,'(3E25.14,i10)') e_b(j)-e_a(1), ABS(dipole(1,j))**2, dipole(1,j), j
           ENDDO
           CLOSE(nbin)
        ENDIF
    

     WRITE(*,'(a60,i10)') 'n_initial = ', na
     WRITE(*,'(a60,i10)') 'n_final   = ', SIZE(e_a) 
     


     WRITE(*,'(a5,1x,4a14)') " n\m","1","2","3","4"
     DO  i = 1, 6
        WRITE(*,'(2X,I3,1X,1P8E16.8)') i, ( dipole( i,j ), j = 1, 6)
     ENDDO


     DEALLOCATE(e_a,e_b,c_a,c_b)
     DEALLOCATE(dipole,dipole_l)

     WRITE(*,*) '@ subroutine evaluate_dipoles out.'

   END SUBROUTINE evaluate_dipoles_diatomics_2




!!%    SUBROUTINE store_radial_psi_bs(ce, psi,r)
!!%      !
!!%      USE utils,                       ONLY: asciifile
!!%      USE input,                       ONLY: method
!!%      !      USE inner_hamiltonian
!!%      !      
!!%      IMPLICIT NONE
!!%      ! arguments
!!%
!!%      REAL(DPK),   INTENT(IN) ::  ce(n_1st:n_last,n_1st:n_last) ! b-splines  coeff :  
!!%      COMPLEX(DPK),INTENT(IN) :: psi(n_1st:n_last)              ! eigenstate coeff : 
!!%      REAL(DPK),   INTENT(in) :: r
!!%                                                                           ! psi = S_n ce_t *Pnl(x)
!!%      ! locals
!!%      COMPLEX(DPK)                            :: p
!!%      REAL(DPK)                               :: cb(bs%n)          ! b-spline coeff
!!%      REAL(DPK)                               :: bvalue       
!!%      INTEGER                                 :: deriv
!!%      INTEGER                                 :: ifile
!!%      INTEGER                                 :: i, j, jj, l
!!%      REAL(DPK)                               :: psi2
!!%      !
!!%      !      INTEGER                                  :: nof_points
!!%      !EXE!
!!%      !
!!%
!!%      !      nof_points = 1000
!!%
!!%      !
!!%      
!!%      ifile  = 1
!!%      DERIV = 0            ! only the function values B(x) are needed
!!%
!!%      CALL asciifile(ifile, it, "psi_x")
!!%
!!%      p = 0.0_dpk
!!%      sum_eigenstates: DO i = n_1st, n_last    ! i eigenstates index    
!!%            
!!%         cb(1) = 0.0D+00 
!!%         DO jj = 1, n_last            !jj B-splines index
!!%            cb(jj + 1) = ce(i, jj)
!!%         END DO
!!%         p = p +  psi(i)  *  bvalue( bs%t, cb, bs%n, bs%k, r, deriv)
!!%      END DO sum_eigenstates
!!%
!!%      WRITE(ifile, '(E15.5, <lmax-lmin+1>E25.15)')r, ABS(p)**2
!!%      
!!%      CLOSE(ifile) 
!!%      DEALLOCATE(p)
!!%
!!%      
!!%       END SUBROUTINE STORE_RADIAL_PSI_BS
!!%         
  END MODULE WF_1E
!!!##########################################################
!!!EOF
             


