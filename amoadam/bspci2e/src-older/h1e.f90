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


PROGRAM H1E

  USE param
  USE set_grid
  USE one_e_matrix
  USE DATA, ONLY: write_mx, write_v_mx, write_v
  USE utils, Only:print_mx
  USE functions
  USE bs_h1e

  !..................

  IMPLICIT NONE
  INTEGER  i,j
  ! I/O
  INTEGER NOUT                            
  !
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: h,b       ! H_ij, B_ij matrices
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: t         ! b-splines grid 
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: en        ! energy vector i = 1,2,..,ns
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: ce        ! coefficients for state en(i)
  !
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: bmx          ! <B_i|x^q|b_j> or <B_i|d/dx|b_j> or <db_i|db_j> 
  INTEGER                                :: deriv, q
  INTEGER,   ALLOCATABLE, DIMENSION(:)   :: nvalence
  !
  !
  !ARG!
  INTEGER                                :: n_a, l_a
  CHARACTER(len=15)                      :: action 
  CHARACTER(len=6)                       :: gauge
!!
  INTERFACE
     SUBROUTINE get_command_line_arg(action,na,la,gauge)
       IMPLICIT NONE
       INTEGER,           INTENT(inout) :: na
       INTEGER,           INTENT(inout) :: la
       CHARACTER(len=15), INTENT(inout) :: action 
       CHARACTER(len=6),  INTENT(inout) :: gauge
     END SUBROUTINE get_command_line_arg
  END INTERFACE
  !EXE!


  nout  = 16  
  l_a = 0
  n_a = 0
  action='solve'
  gauge='no'

  call get_command_line_arg(action,n_a,l_a, gauge)
  CALL input

  OPEN(nout, file='out/h1e.log')

  !,,,,,


  ALLOCATE( t(1:nb + kb)     )
  call mkgrid(t)                     ! make b-splines grid


!  CALL mkgrid_mixed(t)

  
  ! evaluate and store <B_i|O(r)|B_j>     O(r) == d/dr or O(r) = r^q 


  print_B_splines_matrices:IF(action=='test')  THEN
     
     ALLOCATE( bmx(ndim, ndim)  )
 
     WRITE(*,'(a60)') '<B_i|1|B_j> evaluation.'        
     bmx = 0.0_dpk
     deriv = 0
     q = 0
     CALL bsp_integral(nb, bmx, t, deriv, q)        !<B_i|1|B_j>
     CALL print_mx(ndim,bmx,'b_1_b','f')
     !
     WRITE(*,'(a60)') '<B_i|r|B_j> evaluation.'        
     bmx = 0.0_dpk  
     deriv = 0
     q = 1
     CALL bsp_integral(nb, bmx, t, deriv, q)        !<B_i|r|B_j>
     CALL print_mx(ndim,bmx,'b_r_b','f')
     !     
     WRITE(*,'(a60)') '<B_i|d/dr|B_j> evaluation.'        
     bmx   = 0.0_dpk
     deriv = 1
     q = 0
     CALL bsp_integral(nb, bmx, t, deriv, q)      !<B_i|d/dr|B_j>
     CALL print_mx(ndim,bmx,'b_1_db','f')
     !     
     WRITE(*,'(a60)') '<dB_i/dr|dB_j/dr> evaluation.'        
     bmx   = 0.0_dpk
     deriv = 2
     q = 0
     CALL bsp_integral(nb, bmx, t, deriv, q)      !<B_i|d/dr|B_j>
     CALL print_mx(ndim,bmx,'db_1_db','f')
     !
     WRITE(*,'(a60)') '<B_i|1/r|B_j>  evaluation.'        
     bmx   = 0.0_dpk
     deriv = 0
     q = -1
     CALL bsp_integral(nb, bmx, t, deriv, q)        !<B_i|1/r|B_j>
     CALL print_mx(ndim,bmx,'b_1_r_b','f')
     WRITE(*,'(a60)') '<B_i|1/r|B_j>  evaluation.'
     !
     bmx   = 0.0_dpk
     deriv = 0
     q = -2
     CALL bsp_integral(nb, bmx, t, deriv, q)        !<B_i|1/r^2|B_j>
     CALL print_mx(ndim,bmx,'b_1_r2_b','f')

     DEALLOCATE(bmx)
     STOP
  ENDIF print_B_splines_matrices


  evaluate_bsplines_matrices:IF( l_a < 0 ) THEN 
     !
     ! the overlaps are assumed independent on the angular symmetry
     ! implied that for all l's the same B-splines basis is used
     ! 
     ! (nb,kb,tb) == same for all l's
     !

     ALLOCATE( bmx(ndim, ndim)  ) 
     WRITE(*,*) '# allocation of bmx done ndim = ', ndim

     IF(gauge=='no') THEN 

        bmx = 0.0_dpk
        
        deriv = 0
        q = 0
        CALL bsp_integral(nb, bmx, t, deriv, q)        !<B_i|1|B_j>
        CALL write_mx(-1, bmx, "bb-1")
        
     ENDIF
     
     !,,,,,,, length matrix
     
     IF(gauge=='l'.OR.gauge=='lva') THEN 
        
        bmx   = 0.0_dpk
        deriv = 0
        q = 1
        !    

        CALL bsp_integral(nb, bmx, t, deriv, q)        !<B_i|r|B_j>
        CALL write_mx(-1, bmx, "bb-r")
        
        WRITE(*,'(a60)') '<B_i|r|B_j> evaluated.'
     ENDIF
     
     !,,,,,,, velocity matrices

     IF(gauge=='v'.OR.gauge=='lva') THEN

        bmx = 0.0_dpk
        
        deriv = 1
        q = 0
        CALL bsp_integral(nb, bmx, t, deriv, q)      !<B_i|d/dr|B_j>
        CALL write_mx(-1, bmx, "bb-dr")               
        
        bmx = 0.0_dpk
        deriv = 0
        q = -1
        CALL bsp_integral(nb, bmx, t, deriv, q)        !<B_i|1/r|B_j>
        CALL write_mx(-1, bmx, "bb-1r")
        

        WRITE(*,'(a60)') '<B_i|d/dr|B_j> evaluated.'
        WRITE(*,'(a60)') '<B_i|1/r|B_j>  evaluated.'
     ENDIF
  
     !,,,,,,, acceleration matrix    
     
     IF(gauge=='a'.OR.gauge=='lva') THEN
        
        bmx = 0.0_dpk
        deriv = 0
        q = -2
        CALL bsp_integral(nb, bmx, t, deriv, q)         !<B_i|1/r^2|B_j>
        CALL write_mx(-1, bmx, "bb-1r2")
        WRITE(*,'(a60)') '<B_i|1/r^2|B_j>  evaluated.'
     ENDIF
     

     DEALLOCATE(bmx)
     
  ENDIF evaluate_bsplines_matrices


  !solve
  !       ( h_l - e(l) * b ) * c_l =  0
  !
  !
  !,,,,, h_l = -1/2 * d^2/dr^2 + l(l+1)/2r^2 - Z/r   
  !
  !

  make_1e_se_on_bsplines: IF(l_a.GE.0) THEN

     !
     ! evaluate and store <B_i|h0|B_j> , <B_i|B_j>
     !

     ALLOCATE( nvalence( 1 )  )
     ALLOCATE( h(ndim, ndim)  )         ! hamiltonian matrix
     ALLOCATE( b(ndim, ndim)  )         ! b-splines overlap matrix
     !
     nvalence    = 1                    ! evaluate on B-splines hydrogen hamiltonian  
     CALL setmat(nb, l_a, h, b, t)      ! h0 = -d^2/dr^2 + la(la+1)/2r^2 - znuc/r 


     !     CALL print_mx(ndim,b,'b','f')
     !     CALL write_mx(l_a, h, "hb-")
     !     CALL write_mx(l_a, b, "bb-")
     !     CALL print_mx(ndim,h,'h','f')





     get_the_eigensolutions:IF(action=='solve')  THEN

        ALLOCATE( en(ndim) )
        ALLOCATE( ce(ndim, ndim) )     

        CALL solve_banded_diag_nag77(l_a, h, b, en, ce) 
        CALL write_v(l_a,         en, "en1e-")
        CALL write_v(l_a, ce(:,ndim), "wfb-" )
        CALL write_v_mx(l_a, en, ce,  "h1e-" )    
       
        DEALLOCATE( en,ce )
        DEALLOCATE( nvalence )
     ENDIF get_the_eigensolutions
        
     DEALLOCATE( b,h)

     !     
  ENDIF make_1e_se_on_bsplines

  DEALLOCATE( t )
  CLOSE(nout)



  !
  ! code for the free-boundaries solution of SE (linear equation solver)
  !


!!%  !,,,,,,,,,             solve 1-e SE for partial wave l
!!%
!!%  ALLOCATE( C0(NB-1, 1)    )
!!%  C0          = 0.0D+00
!!%  C0(NB-1, 1) = 1.0D+00
!!%     
!!%  IF(method.EQ.'l') THEN  ! linear equations method    ( P(R) free )
!!%      
!!%     ALLOCATE( en(ndim) )
!!%     ALLOCATE( ce(ndim, nb-1) )
!!%     
!!%     CALL solve_banded_diag_nag77(l, h, b, en, ce) 
!!%     CALL write_v(l,         en, "en1e-")
!!%     CALL write_v(l, ce(:,nb-1), "wfb-" )
!!%
!!%  ELSE IF(method == 'd') THEN                        ! diagonalization method ( P(R) fxd  (==0) )
!!%
!!%     ns = nbs + ncs
!!%
!!%     ALLOCATE( EN(NS) )
!!%     ALLOCATE( CE(NS, NB-1) )
!!%
!!%     CALL ENERGY_SPECTRUM(L, EN, 'OUT')              ! calculate en(i) and store
!!%     CALL SOLVE_BANDED_LIN_NAG77(L, H, B, C0, EN, CE)
!!%
!!%
!!%
!!%  ELSE IF(method == 'dl') THEN ! diagonalization (bound states) + linear (continuum states) method
!!%     
!!%     ALLOCATE( EN_B(NB-2)       )
!!%     ALLOCATE( CE_B(NB-2, NB-1) )
!!%
!!%     CALL SOLVE_BANDED_DIAG_NAG77(L, H, B, EN_B, CE_B)   !solve for en<0 (bound spectrum)
!!%
!!%     ! get nof of bound states by fxd conditions
!!%
!!%     nbs = 0
!!%     DO i = 1, SIZE(en_b)
!!%        IF(en_b(i) < 0.0D+00) nbs = nbs + 1
!!%     ENDDO
!!%
!!%     WRITE(*,*) ' diag nbs = ', nbs
!!%
!!%     IF(spectrum=='fxd') THEN 
!!%
!!%        ! note that here we need only eigenvalues of the fxd problem. To be accomodate it
!!%
!!%        ns = SIZE(en_b)
!!%
!!%        WRITE(*,*) '         ns = ', ns
!!%
!!%        ALLOCATE( EN(NS) )
!!%        ALLOCATE( CE(NS, NB-1) )
!!%
!!%        en = en_b
!!%
!!%        CALL SOLVE_BANDED_LIN_NAG77(L, H, B, C0, EN, CE)  !solve for en>0  (continuum spectrum)
!!%
!!%        ce(1:nbs,:)    = ce_b(1:nbs,:)
!!%
!!%        DEALLOCATE( EN_B )    !en_b
!!%        DEALLOCATE( CE_B )    !ce_b
!!%
!!%     ELSE IF (spectrum=='mxd') THEN
!!%
!!%        ALLOCATE( EN_C(NCS) )
!!%        ALLOCATE( CE_C(NCS, NB-1) )
!!%
!!%        CALL ENERGY_SPECTRUM(L, EN_C, 'INOUT')               ! calculate en_c(i) (do not store)
!!%        CALL SOLVE_BANDED_LIN_NAG77(L, H, B, C0, EN_C, CE_C) !solve for en>0  (continuum spectrum)
!!%        ns = nbs + ncs
!!%        
!!%        WRITE(*,*) "# h1e::                       nbs = ", nbs
!!%        WRITE(*,*) "# h1e::                       ncs = ", ncs
!!%
!!%        ALLOCATE( EN(NS) )
!!%        ALLOCATE( CE(NS, NB-1) )
!!%        
!!%! energy spectrum (b + c)
!!%
!!%        en(     1 : nbs ) = en_b( 1 : nbs)
!!%        en( nbs+1 :  ns ) = en_c( 1 : ncs)
!!%
!!%! coefficients    (b + c)
!!%
!!%        ce(    1: nbs,  : ) = ce_b( 1 : nbs, : )
!!%        ce(nbs+1:  ns , : ) = ce_c( 1 : ncs, : )
!!%
!!%        DEALLOCATE( EN_B )    !en_b
!!%        DEALLOCATE( CE_B )    !ce_b
!!%        DEALLOCATE( EN_C )    !en_c
!!%        DEALLOCATE( CE_C )    !ce_c
!!%     ENDIF
!!%     CALL ENERGY_SPECTRUM(L, EN, 'IN')  ! get en(i) and store
!!%  ENDIF
!!%
!!% 
!!%  CALL WRITE_V_MX(L, EN, CE, "h1e-")
!!%  DEALLOCATE( C0 )

!  CLOSE(nout)
!  DEALLOCATE( nvalence)
!  DEALLOCATE( t )
!  DEALLOCATE( h )
!  DEALLOCATE( b )
!  DEALLOCATE( en)
!  DEALLOCATE( ce)
     

  WRITE(*,*) '@ main program h1e end.'

END PROGRAM h1e
!
!
!
SUBROUTINE get_command_line_arg(action,na,la, gauge)
  !
  USE PRECISION, ONLY: DPK
  !
  IMPLICIT NONE
  !ARG!
  INTEGER,           INTENT(inout) :: na
  INTEGER,           INTENT(inout) :: la
  CHARACTER(len=15), INTENT(inout) :: action 
  CHARACTER(len=6),  INTENT(inout) :: gauge
  !LOC!
  CHARACTER*180 LINE, EXE, WHAT
  CHARACTER*40  CMD 
  INTEGER       NARG,IARG, NXT_ARG
  INTEGER       LENGTH, ISTATUS
  INTEGER       COMMAND_ARGUMENT_COUNT
  CHARACTER(len=10)       :: date,time,zone
  INTEGER, DIMENSION(dpk) :: values 
  INTEGER nascii
  !EXE!
 
  nascii = 1
  
  CALL DATE_AND_TIME(date, time, zone, values)
  
  NARG = COMMAND_ARGUMENT_COUNT()    
  CALL GET_COMMAND( LINE, LENGTH, ISTATUS)     
  CALL GET_COMMAND_ARGUMENT(0,EXE,LENGTH,ISTATUS) 


  WRITE(*,*)'#'
  WRITE(*,'(a45,1X,a40)')'# w1e::            executable command exe = ', line
  WRITE(*,'(a45,1X,i2)') '# wf1e::        nof arguments        narg = ', narg  
  IF(narg.LE.1) THEN
     WRITE(*,'(a45,1X,i3 )') '# wf1e::                                na = ', na
     WRITE(*,'(a45,1X,i3)')  '# wf1e::                                la = ', la
     WRITE(*,'(a45,1X,a40)') '# wf1e::                            action = ', action
     WRITE(*,'(a45,1X,a40)') '# wf1e::                             gauge = ', gauge
  ENDIF

  nxt_arg = 1
  get_arguments: DO iarg = 1, narg

     CALL GET_COMMAND_ARGUMENT(IARG,CMD,LENGTH,ISTATUS)

     IF(MOD(iarg,2)==1) THEN    !iarg = 1, 3, 5, ...

        READ(cmd,*) what

        IF(what=='-o') THEN 
           WRITE(*,*) '# wf1e::'
           WRITE(*,*) '# wf1e:: Available options:'
           WRITE(*,*) '# wf1e:: 0. -o information '
           WRITE(*,*) '# wf1e:: 1. n= eigenstate index (integer,   > 0),        [1]'
           WRITE(*,*) '# wf1e:: 2. l= angular symmetry (integer,  >= 0),        [0]'
           WRITE(*,*) '# wf1e:: 3. a= action           (character,l=40),  [save]'
     
      WRITE(*,*) '# wf1e::        (save)   : stores all eigenstates'
           WRITE(*,*) '# wf1e:: 6. g= gauge           (character,l=6),  [l],v,a,lva'
           WRITE(*,*) '# wf1e:: happy end.'  
           STOP
        ENDIF

        nxt_arg = iarg + 1

        WRITE(*,'(a45,1X,a40)') '# wf1e::                             what = ', what
        IF(     (what.NE.'na='        ).AND.     &
             &  (what.NE.'la='        ).AND.     &
             &  (what.NE.'a='        ).AND.     &
             &  (what.NE.'g='        ) ) THEN   
           WRITE(*,'(a60)') 'available options for :'
           WRITE(*,'(a60)') ' (n= principal q. number)'
           WRITE(*,'(a60)') ' (l= angular momenta)'
           WRITE(*,'(a60)') ' (a= (save) )'
           WRITE(*,'(a60)') ' (g= (l),v,a,lva)'
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
        ELSE IF(what=='a=') THEN
           action = trim(cmd)
           IF(  (action.NE.'save'.and.action.ne.'solve'.and.action.ne.'test') ) THEN
              WRITE(*,'(a60)') ' available options for a='
              WRITE(*,'(a60)') ' (save)'
           ENDIF
           WRITE(*,'(a45,1X,a40)') '# wf1e::                           action = ', action
        ELSE IF(what=='g=') THEN
          gauge = trim(cmd)
          WRITE(*,'(a45,1X,a40)') '# wf1e::                             gauge = ', gauge
        ELSE
           CYCLE
        ENDIF
        
     ENDIF
  ENDDO get_arguments
  
  !............... save history
  OPEN(nascii, file='log/h1e_history.log', position='append')
  WRITE(nascii,'(a10,1x,a10,2x,a60)') date,time, line    
  CLOSE(nascii)


  RETURN
END SUBROUTINE get_command_line_arg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! EOF


