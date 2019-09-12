!
!  6-9/12/2007/QUB
!
! several adds: (1) overlap integrals  <P_nl|P_n'l> (GL integration) 
!               (2) projection coefficients of a function f(x) on 
!                   the B-splines basis and/or the P_nl eigenstates
!               (3) radial dipole matrix elements <l|r|l+1> (GL integration)
!
! use as: (all integrations are performed with a GL integration quadrature)
!
! >> Rw1e la= 0         (saves all P_nl(x), n = 1,...,n_max in hwf1e-l.dat)
!
!
! >> Rw1e la= 0 na= 1    (output only the P_(10)(x) radial wavefunction)
!                       (necessitates the Rwf1e l= 0 to be run first  ) 
!
!
! >> Rw1e l= 0 a= orthornomality (calculates all the overlaps integrals <P_nl|P_ml> for the symmetry l)
!                                  it should be found as  delta_nm = <P_nl|P_ml>
!                        
!
! >> Rw1e l = 0 a= projection (calculates the projection coefficient for a given function f(x) 
!                               over the B-splines and over the P_nl(x) orthonormal basis 
!
!                              cb_f(i) = < f(x)|B_i(x) >     
!                              ce_f(n) = < f(x)|P_nl(x) >     
!
!
! >> Rw1e  l= 0 a= dipole (saves all radial dipole matrix elements <nl|r|ml+1> )
!
!


PROGRAM wf1e

  USE param,        only: nl
  USE PRECISION,    ONLY: dpk
  USE DATA,         ONLY: read_mx, read_mx_p, read_v_mx, write_mx, read_v
  USE wf_1e_adam,        ONLY: read_wf1e, write_wf1e, write_lopt_wf1e, write_all_wf1e_adam, save_wf1e_nl, write_fx, write_fx_bsp
  USE wf_1e_adam,        ONLY: evaluate_dipoles, evaluate_overlaps, evaluate_overlap_dipoles
  USE wf_1e_adam,        ONLY: evaluate_dipoles_diatomics, evaluate_dipoles_diatomics_2
  USE functions,    ONLY: fx_gauss_real, fx_gauss_imag
  USE utils,        ONLY: print_mx, banded_mat_v_mul, asciifile
  USE one_e_matrix !, ONLY: bsp_overlap
  USE my_types,     ONLY: state_index, basis, init_basis
  !
  IMPLICIT NONE
  !
  !
  TYPE(state_index)                      :: state_a
  TYPE(state_index)                      :: state_b
  TYPE(basis)                            :: bsplines
  !
  REAL(DPK), DIMENSION(:,:), POINTER     :: ce          ! < B_i  | P_nl(x) > 
  REAL(DPK), DIMENSION(:),   POINTER     :: en          !
  !  REAL(DPK), DIMENSION(:,:), POINTER     :: pl,dpl      ! p_l(r), p'_l(r)
  !  REAL(DPK), DIMENSION( :),  POINTER     :: r, dr       ! r_i, dr_i
  !
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: b           ! < B_i  | B_j  > 
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: cb_f        ! < B_i  | f(x) > 
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: ce_f        ! < P_nl | f(x) > 
  !
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: vr          ! 
  INTEGER                                :: icn, ndim, ne
  INTEGER                                :: i,j,l
  !
  INTEGER                                :: nascii, ifile
  !
  ! command line arguments
  INTEGER                                :: n_a, l_a, m_a, p_a
  INTEGER                                :: n_b, l_b, m_b, p_b
  CHARACTER(len=15)                      :: action 
  CHARACTER(len=6)                       :: gauge
  CHARACTER(len=6)                       :: interaction
  !
  INTERFACE
     SUBROUTINE get_command_line_arg(action,nof_l,na,la,nb,lb,gauge,interaction)
       IMPLICIT NONE
       INTEGER,           INTENT(inout) :: nof_l
       INTEGER,           INTENT(inout) :: na
       INTEGER,           INTENT(inout) :: nb
       INTEGER,           INTENT(inout) :: la
       INTEGER,           INTENT(inout) :: lb
       CHARACTER(len=15), INTENT(inout) :: action 
       CHARACTER(len=6),  INTENT(inout) :: gauge
       CHARACTER(len=6),  INTENT(inout) :: interaction
     END SUBROUTINE get_command_line_arg
  END INTERFACE
  !EXE!

  ! preliminaries

  nascii = 1
  ifile  = 2   
  

  ! default parameters (if not command line present)

  nl      = 1
  !
  n_a     = 1
  l_a     = 0

  !
  n_b     = 1
  l_b     = l_a + 1

  !
  action  ='save'
  gauge   ='l'
  interaction ='lopt'

  !
  m_a = 0
  m_b = 0
  p_a = 1
  p_b = 1
  !

  ! command line parameters

  CALL get_command_line_arg(action, nl, n_a,l_a, n_b,l_b, gauge, interaction)

  !init basis 

  CALL init_basis(bsplines)       
  

  ! main program statements


!
!!!     save      E_nl, P_nl(x) (in the physical grid) in binary format
!
 
  save_wf1e_binary_grid:IF(action=='save') THEN     ! Rwf1e l= 0  


     CALL read_v_mx( l_a, en, ce, "h1e-")           ! read en + coeff
     CALL asciifile(ifile,l_a,"wb-")                ! write in ascii
     !
     WRITE(ifile,*) SIZE(en)
     DO i= 1, SIZE(en)
        WRITE(ifile,'(2E25.14)') en(i),ce(i,SIZE(ce,dim=2))
     ENDDO
     CLOSE(ifile)
     !
     WRITE(*,'(a60, i10)')  ' save::   wb-l files saved in dat/ :'
     WRITE(*,'(a60, i10)')  ' save::                       l = ', l_a
     WRITE(*,'(a60, 2i10)') ' save::                     dim = ', SIZE(ce,dim=2), SIZE(en)
     WRITE(*,'(a60, 2i10)') ' save::                      nl = ', nl
     WRITE(*,'(a60, 2i10)') ' save::                  na, la = ', n_a, l_a


     !      CALL write_wf1e(n_a, l_a, en, ce)    ! store radial wfs p_nl(r),p'_nl(r)

     ! used from v2e program
     CALL write_lopt_wf1e(n_a, l_a, en, ce)   ! store radial wfs p_nl(r),p'_nl(r)

     !used from tdse_bs_rho program
     !CALL write_all_wf1e_adam(n_a, l_a, en, ce)   ! store radial wfs p_nl(r),p'_nl(r)
        
WRITE(*,'(a60, i10)')  ' I used adams code'
     DEALLOCATE( en, ce )
     
  ENDIF save_wf1e_binary_grid


!
!
!
!
!

  merge_target_energy_files:IF(action=='merge') THEN

     !determine size by reading the l_a file 
     !(not very elegant but quick)

     CALL asciifile(ifile,l_a,"en1e-")
     READ(ifile,*) icn
     CLOSE(ifile)

     ALLOCATE(vr(icn))
     icn  = 1
     OPEN(15, file='dat/en1e.dat')
     WRITE(15,'(6i5)') l_a+1, l_b+1
     DO l = l_a, l_b
        CALL read_v(l,vr,'en1e-')
        
        WRITE(15,'(6i5)') SIZE(vr), icn
        DO i= 1, SIZE(vr)
           WRITE(15,'(E25.14)') vr(i)*2.0_dpk
        ENDDO
     ENDDO
     CLOSE(15)

  ENDIF merge_target_energy_files

  !
  !     save P_nl(x) in physical grid in ascii     ! first run as (Rwf1e l= 0  save)
  !

  !  get_wf1e_nl:IF (na.GT.0) THEN 
  !     CALL READ_WF1E(la, pl, dpl, en, r, dr)
  !   CALL save_wf1e_nl(na, la, en(n), pl(n,:), dpl(n,:), r, SIZE(r))
  !   WRITE(*,*) '# wf1e::     save P(n,l) = ', na,la
  !  ENDIF get_wf1e_nl
  

  !
  ! calculate overlap integrals:  < P_nl | P_n'l >
  !

                                                     
  calculate_orthnormality_wf1e:IF(action=='orthonormality') THEN    !GL integration: 
     
     CALL READ_V_MX(l_a, en, ce, "h1e-")     ! read en + coeff     
     ndim = SIZE(ce(1,:))
     
     ALLOCATE( b(ndim,ndim) )
     CALL READ_MX(l_a, b, "bb-")             ! read (upper part) Bsp-overlap mx

     ALLOCATE( vr(ndim)     )                ! temporary storage vector

     CALL asciifile(nascii,l_a,'wf1e-orth-')

     perform_v_mat_v_mul: do i = 1, SIZE(en)                     ! <c_l | B | c_r >
        DO j = i, SIZE(en)
           vr = ce(j,:)
           CALL banded_mat_v_mul(b, vr, bsplines%k)                                ! c_r <--- b*c_r 
           WRITE(nascii,'(2I4,3X,1E30.15)') i, j, DOT_PRODUCT(ce(i,:),vr)  ! c_l*c_r 
        ENDDO

     ENDDO perform_v_mat_v_mul
     
     DEALLOCATE(en,ce,b,vr)
     CLOSE(nascii)
        
  ENDIF calculate_orthnormality_wf1e


  !
  !
  ! calculate dmx integrals:  < P_nl | T_g | P_n'l' >,  T_G =  l,v,a
  !
  !

    
  calculate_dipoles_pp:IF(action=='dipole') THEN       !GL integration! 


     IF((l_b-l_a).NE.1) THEN       !if not the default case

        all_dipoles_partial_waves: DO l = l_a, l_b-1
           
           state_a = state_index(n_a, l,  m_a, p_a, nl)
           state_b = state_index(n_b, l+1,m_b, p_b, nl)
           !
           IF(gauge=='lva') THEN
              gauge = 'l'
              CALL evaluate_dipoles(state_a, state_b, gauge, bsplines, interaction)
              gauge = 'v'
              CALL evaluate_dipoles(state_a, state_b, gauge, bsplines, interaction)
              gauge = 'a'
              CALL evaluate_dipoles(state_a, state_b, gauge, bsplines, interaction)
           ELSE 
              CALL evaluate_dipoles(state_a, state_b, gauge, bsplines, interaction)
           ENDIF
           
        ENDDO all_dipoles_partial_waves
        
     ELSE
        
           
        state_a = state_index(n_a,l_a,m_a,p_a, nl)
        state_b = state_index(n_b,l_b,m_b,p_b, nl)
           !
        IF(gauge=='lva') THEN
           gauge = 'l'
           CALL evaluate_dipoles(state_a,state_b,gauge,bsplines, interaction)
           gauge = 'v'
           CALL evaluate_dipoles(state_a,state_b,gauge,bsplines, interaction)
           gauge = 'a'
           CALL evaluate_dipoles(state_a,state_b,gauge,bsplines, interaction)
        ELSE 
           CALL evaluate_dipoles(state_a,state_b,gauge,bsplines, interaction)
        ENDIF
        

     ENDIF
     

  ENDIF calculate_dipoles_pp
  !
  !
  !
  !
  calculate_dipoles_diatomic:IF(action=='mdipole') THEN     

     
     state_a = state_index(n_a, l_a,     m_a, p_a, nl)
     state_b = state_index(n_b, l_a + 1, m_b, p_b, nl)
     !
     gauge = 'l'
     CALL evaluate_dipoles_diatomics(state_a, state_b, gauge, bsplines, interaction)


!     IF(gauge=='lva') THEN
!        gauge = 'l'
!        CALL evaluate_diatomic_dipoles(state_a, state_b, gauge, bsplines, interaction)
!        gauge = 'v'
!        CALL evaluate_diatomic_dipoles(state_a, state_b, gauge, bsplines, interaction)
!        gauge = 'a'
!        CALL evaluate_diatomic_dipoles(state_a, state_b, gauge, bsplines, interaction)
!     ELSE 
!        CALL evaluate_diatomic_dipoles(state_a, state_b, gauge, bsplines, interaction)
!     ENDIF

  ENDIF calculate_dipoles_diatomic




     !
     ! overlap:    <P_(na,la) | B_j >          (P_(na,la) fxd orbital P(R) = 0)
     !                          na = 1,..,nb-2, 
     !                          j  = 1,..,nb
     !
     !  P_(na,la):    fxd orbital:    P(R) = 0
     !



  calculate_overlaps_pb:IF(action=='overlap') THEN
     

     IF((l_b-l_a).NE.1) THEN       !if not the default case

        all_overlaps_partial_waves: DO l = l_a, l_b

           state_a = state_index(n_a, l, m_a, p_a, nl)
           CALL evaluate_overlaps(state_a,bsplines)

        ENDDO all_overlaps_partial_waves
     ELSE

        state_a = state_index(n_a, l_a, m_a,p_a, nl)
        CALL evaluate_overlaps(state_a, bsplines)
     ENDIF


     
  ENDIF calculate_overlaps_pb



  !
  !
  ! calculate dmx-overlap integrals:  < P_nl | T_g | B_j >,  T_G =  l,v,a
  !
  !           n = 1,...,nb-2
  !           j = 1,....nb
  !
  !           P_(na,la):    fxd orbital:    P(R) = 0
  !
  !          T_G = 1,lva
  !



  calculate_overlap_dipoles_pb:IF(action=='overlap_dipole') THEN       !GL integration! 


     IF((l_b-l_a).NE.1) THEN       !if not the default case

        all_overlap_dipoles_partial_waves: DO l = l_a, l_b-1

           state_a = state_index(n_a, l,  m_a, p_a, nl)
           state_b = state_index(n_b, l+1,m_b, p_b, nl)
     
           CALL evaluate_overlap_dipoles(state_a, state_b,gauge, bsplines)
           CALL evaluate_overlap_dipoles(state_b, state_a,gauge, bsplines)

           IF(gauge=='lva') THEN
              gauge = 'v'
              CALL evaluate_overlap_dipoles(state_a, state_b,gauge, bsplines)
              CALL evaluate_overlap_dipoles(state_b, state_a,gauge, bsplines)
              gauge = 'a'
              CALL evaluate_overlap_dipoles(state_a, state_b,gauge, bsplines)
              CALL evaluate_overlap_dipoles(state_b, state_a,gauge, bsplines)
           ENDIF
        ENDDO all_overlap_dipoles_partial_waves

     ELSE ! only the la --> lb (lb-->la) transitions

        state_a = state_index(n_a, l_a, m_a, p_a, nl)
        state_b = state_index(n_b, l_b, m_b, p_b, nl)
     
        CALL evaluate_overlap_dipoles(state_a, state_b,gauge, bsplines)
        CALL evaluate_overlap_dipoles(state_b, state_a,gauge, bsplines)
        
        IF(gauge=='lva') THEN
           gauge = 'v'
           CALL evaluate_overlap_dipoles(state_a, state_b,gauge, bsplines)
           CALL evaluate_overlap_dipoles(state_b, state_a,gauge, bsplines)
           gauge = 'a'
           CALL evaluate_overlap_dipoles(state_a, state_b,gauge, bsplines)
           CALL evaluate_overlap_dipoles(state_b, state_a,gauge, bsplines)
        ENDIF
           
     ENDIF

  ENDIF calculate_overlap_dipoles_pb

  !
  ! calculate projection coefficients cb_i and ce_i of a function f(x):
  !
  !               f(x) = S_i cb_i * B_i(x)        (1)
  !
  !               f(x) = S_i ce_i * P_nl(x)       (2) 
  !


  calculate_projection_wf1e:IF(action=='projection') THEN    !GL integration: < P_nl|P_n'l>


     !# CALL evaluate_projection(state_a,fx)


     ! calculate projection of f(x) over the non-orthogonal bsplines B_i(x) 

                                                     ! P_nl(x) = S_i ce_i B_i(x)
     CALL READ_V_MX( l_a, en, ce, "h1e-")              !   cb_n = < B_i | P_nl >
     
     ndim = SIZE(ce(1,:))
     ne   = SIZE(en)

     ALLOCATE(   cb_f(ndim) )                        ! projection over B-splines
     ALLOCATE( b(ndim,ndim) )
     
     CALL read_mx(l_a, b, "bb-")                     ! read (upper part) Bsp-overlap mx
     CALL bsp_overlap(fx_gauss_real, b, cb_f)        ! f(x) = S_i cb_f(i) * B_i(x)
     CALL write_fx_bsp(cb_f)                         ! store radial wfs p_nl(r),p'_nl(r)


     ! calculate projection of f(x) over orthogonal eigenstates P_nl(x) 
     
     ALLOCATE( ce_f(ne) )
     ALLOCATE( vr(ndim) )                            ! 

     CALL asciifile(nascii,l_a,'wf1e-proj-')


     vr = 0.0D+00
     perform_projection: DO i = 1, ne                 ! f(x) = S_n cb_n(n) P_n(x)
                                                      ! cb_n = < P_nl | f(x) >
        vr = cb_f
        CALL banded_mat_v_mul(b, vr, bsplines%k)              ! c_r <-- b * c_r 
        ce_f(i) = DOT_PRODUCT( ce(i,:), vr)           ! ce_f = <ce(i,:)|B|cb_f> 
        WRITE(nascii,'(1I4,3X,1E30.15)') i, ce_f(i)

     ENDDO perform_projection
     
     CALL write_fx(l_a, ce, ce_f)
     
     DEALLOCATE(en,ce,b,vr,cb_f,ce_f)
     CLOSE(nascii)
     
  ENDIF calculate_projection_wf1e

!CALL OVERLAP(l, pl(i,:), pl(j,:), r, dr, nx, integral )


!..............

END PROGRAM WF1E
!
!
!
SUBROUTINE get_command_line_arg(action,nof_l, na,la,nb,lb,gauge,interaction)
  !
  USE PRECISION, ONLY: DPK
  !
  IMPLICIT NONE
  !ARG!
  INTEGER,           INTENT(inout) :: nof_l
  INTEGER,           INTENT(inout) :: na
  INTEGER,           INTENT(inout) :: nb
  INTEGER,           INTENT(inout) :: la
  INTEGER,           INTENT(inout) :: lb
  CHARACTER(len=15), INTENT(inout) :: action 
  CHARACTER(len=6),  INTENT(inout) :: gauge
  CHARACTER(len=6),  INTENT(inout) :: interaction
  !LOC!
  CHARACTER*180 LINE, EXE, WHAT
  CHARACTER*40  CMD 
  INTEGER       NARG,IARG, NXT_ARG
  INTEGER       LENGTH, ISTATUS
  INTEGER       COMMAND_ARGUMENT_COUNT
  CHARACTER(len=10)               :: date,time,zone
  INTEGER, DIMENSION(dpk)         :: values 
  INTEGER nascii
  !EXE!
 
  nascii = 1
  
  CALL DATE_AND_TIME(date, time, zone, values)
  
  NARG = COMMAND_ARGUMENT_COUNT()    
  CALL GET_COMMAND( LINE, LENGTH, ISTATUS)     
  CALL GET_COMMAND_ARGUMENT(0,EXE,LENGTH,ISTATUS) 


  WRITE(*,*)'#'
  WRITE(*,'(a45,1X,a40)')'# w1e::            executable command    exe = ', line
  WRITE(*,'(a45,1X,i2)') '# w1e::                 nof arguments   narg = ', narg
  WRITE(*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'    
  WRITE(*,*)'#'

  IF(narg.LE.1) THEN
     WRITE(*,'(a45,1X,i3 )') '# wf1e::                                na = ', na
     WRITE(*,'(a45,1X,i3)')  '# wf1e::                                la = ', la
     WRITE(*,'(a45,1X,i3)')  '# wf1e::                                nb = ', nb
     WRITE(*,'(a45,1X,i3 )') '# wf1e::                                lb = ', lb
     WRITE(*,'(a45,1X,a40)') '# wf1e::                            action = ', action
     WRITE(*,'(a45,1X,a40)') '# wf1e::                             gauge = ', gauge
     WRITE(*,'(a45,1X,a40)') '# wf1e::                       interaction = ', interaction
     WRITE(*,'(a45,1X,i3 )') '# wf1e::                             nof_l = ', nof_l
  ENDIF

  nxt_arg = 0
  get_arguments: DO iarg = 1, narg

     CALL GET_COMMAND_ARGUMENT(IARG,CMD,LENGTH,ISTATUS)

     IF(MOD(iarg,2)==1) THEN    !iarg = 1, 3, 5, ...

        READ(cmd,*) what

        IF(what=='-o') THEN 
           WRITE(*,*) '# w1e::xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
           WRITE(*,*) '# w1e:: Available options:'
           WRITE(*,*) '# w1e:: 0. -o information '
           WRITE(*,*) '# w1e:: 1. na= eigenstate index (integer,   > 0),        [1]'
           WRITE(*,*) '# w1e:: 2. la= angular symmetry (integer,  >= 0),        [0]'
           WRITE(*,*) '# w1e:: 3. nb= eigenstate index (integer,   > 0),        [1]'
           WRITE(*,*) '# w1e:: 4. lb= angular symmetry (integer,  >= 0),     [la+1]'
           WRITE(*,*) '# w1e:: 5. a= action           (character,l=40),  [save]'
           WRITE(*,*) '# w1e::        (save)   : stores all eigenstates on a physical grid'
           WRITE(*,*) '# w1e::(orthornormality): calculates all <P_nl| 1 |P_ml> integrals'  
           WRITE(*,*) '# w1e::         (dipole): calculates all <P_nl| r |P_ml> integrals'  
           WRITE(*,*) '# w1e::        (mdipole): calculates all diatomic dipoles S_{l,lp} = k_(l,lp)<f_nl| r |f_nplp> '  
           WRITE(*,*) '# w1e::     (projection): calculates coeff f(x) = S_i c_n P_nl(x)'  
           WRITE(*,*) '# w1e::        (overlap): calculates all <P_nl|B_j> integrals'  
           WRITE(*,*) '# w1e:: (overlap_dipole): calculates all <P_nl|T_G|B_j> integrals, g=l,v,a'  
           WRITE(*,*) '# w1e:: 6. g= gauge           (character,l=6),                 [l],v,a,lva'
           WRITE(*,*) '# w1e:: 7. i= interaction     (character,l=6),   save data for [tdse],lopt'
           WRITE(*,*) '# w1e:: 8. nl= nof of coupled partial waves           (integer, >= 1), [1]'
           WRITE(*,*) '# w1e:: happy end.'  
           STOP
        ENDIF

        nxt_arg = iarg + 1

        WRITE(*,'(a45,1X,a40)') '# wf1e::                 what = ', what
        IF(     (what.NE.'na=').AND.     &
             &  (what.NE.'la=').AND.     &
             &  (what.NE.'nb=').AND.     &
             &  (what.NE.'lb=').AND.     &
             &  (what.NE.'nl=').AND.     &
             &  (what.NE.'a=').AND.     &
             &  (what.NE.'g=').AND.     &
             &  (what.NE.'i=')) THEN     
           WRITE(*,'(a60)') ' available options for :'
           WRITE(*,'(a60)') ' (na= principal q. number)'
           WRITE(*,'(a60)') ' (la= angular momenta)'
           WRITE(*,'(a60)') ' (nb= principal q. number)'
           WRITE(*,'(a60)') ' (lb= angular momenta)'
           WRITE(*,'(a60)') ' (nl= nof coupled partial waves)'
           WRITE(*,'(a100)') '(a= (save), orthornormality, projection, dipole, mdipole, overlap, overlap_dipole)'
           WRITE(*,'(a60)') ' (g= (l),v,a,lva)'
           WRITE(*,'(a60)') ' (i= (tdse),lopt)'
           STOP
        ENDIF
         
        CYCLE
     
     ENDIF

!...................     

     IF(iarg == nxt_arg ) THEN
        IF(what=='na=') THEN 
           READ(cmd,'(I3)')  na    
           WRITE(*,'(a45,1X,i3)') '# wf1e::                                na = ', na
        ELSE IF(what=='la=') THEN
           READ(cmd,*)  la 
           WRITE(*,'(a45,1X,i3 )') '# wf1e::                               la = ', la
        ELSE IF(what=='nb=') THEN
           READ(cmd,*)  nb 
           WRITE(*,'(a45,1X,i3 )') '# wf1e::                               nb = ', nb
        ELSE IF(what=='lb=') THEN
           READ(cmd,*)  lb 
           WRITE(*,'(a45,1X,i3 )') '# wf1e::                               lb = ', lb
        ELSE IF(what=='nl=') THEN
           READ(cmd,*)  nof_l 
           WRITE(*,'(a45,1X,i3 )') '# wf1e::                               nl = ', nof_l
        ELSE IF(what=='a=') THEN
           action = cmd
           IF(  (action.NE.'overlap'         ).AND.&
                (action.NE.'projection'      ).AND.&
                (action.NE.'save'            ).AND.&
                (action.NE.'overlap_dipole'  ).AND.&
                (action.NE.'orthonormality'  ).AND.&
                (action.NE.'mdipole').AND.&
                (action.NE.'dipole'     )) THEN
              WRITE(*,'(a60)') ' available options for a='
              WRITE(*,'(a60)') ' (save)'
              WRITE(*,'(a60)') ' (dipole)'
              WRITE(*,'(a60)') ' (mdipole)'
              WRITE(*,'(a60)') ' (projection)'
              WRITE(*,'(a60)') ' (overlap)'
              WRITE(*,'(a60)') ' (orthornormality)'
              WRITE(*,'(a60)') ' (overlap_dipole)'
           ENDIF
           WRITE(*,'(a45,1X,a40)') '# wf1e::                           action = ', action
        ELSE IF(what=='g=') THEN
          gauge = cmd
          WRITE(*,'(a45,1X,a40)') '# wf1e::                           gauge = ', gauge
        ELSE IF(what=='i=') THEN
          interaction = cmd
          WRITE(*,'(a45,1X,a40)') '# wf1e::                      interaction = ', interaction
        ELSE
           CYCLE
        ENDIF
        
     ENDIF
  ENDDO get_arguments
  
  !............... save history
  OPEN(nascii, file='log/wf1e_history.log', position='append')
  WRITE(nascii,'(a10,1x,a10,2x,a60)') date,time, line    
  CLOSE(nascii)


  RETURN
END SUBROUTINE get_command_line_arg
!eof
             
