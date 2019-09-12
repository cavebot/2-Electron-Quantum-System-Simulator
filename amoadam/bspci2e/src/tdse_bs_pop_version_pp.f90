!
! [1] PRA 68,013409, Laulan and Bachau,
! 'Correlation effects in two-photon single and Double ionization of helium'
!

!
! calculates ionization probabilities for the 2-e Hamiltonian problem
!
!
!  id/dt psi(r1,r2;t) = [ H(r1,r2) + V(r1,r2;t)] psi(r1,r2;t) 
!
!    H = h(r1) + h(r2) + 1/|r_1-r_2| = H_0 + V12
!
!    if |i> defined as the CI 2e states (E;LS)  (ignore here M_L,M_S)
!
!    H|i> = E |i>,  S_i |i><i| = 1 and <i|i'> = delta_ii' 
!
!    and |j> the zero-order 2e states (uncorrelated)    (n1l1;n2l2;E_0=e1+e2;LS)
!
!          H_0|j> = E_0 |j>,  S_i |j><j| = 1 and <j|j'> = delta_jj', E_0 = e_1 + e_2
!
!     with e1,e2 the eigenenergies of the 1-e SE:
!
!         h(r1)|e1> = e1 |e1> and h(r2) |e2> = e2|e2>
!
! then this program calculates two different ionization probabilities are calculated
!
! - 
!
!   (1) psi(t) = S_i c_i(t) |i>,   TDSE calculation provides c_i(t) due to the laser interaction
!
!  with    i == (nL) 
!
!  |i> = S_{j} v_ij |j>,           CI calculation provides elements v_ij of the V12 interaction
!
!    j  = n1l1;n2l2 and Phi_j the zero-order  
! v_ij  = the  CI-mixing coefficients
! 
!  |j > = the zero-order (uncorrelated) 2e-hamiltonian
!
! we have two expansions for the TD wf:
!
!  psi(t) = S_i c_i(t) |i> = S_j v_j(t) |j>,   v_j = S_i c_i(t) v_ij
!
!
! then
!       P_j(t) = |<j|psi(t)>| = ...= |S_i c_i(t)v_ij|^2
!
!  probability that the system will be found at state |j>
!
!  with energy E= e_1 + e_2, 
!
! and
!
!      P_i(t)  = |<i|psi(t)>| = ...= |c_i(t)|^2
!
! probability that the system will be found in state |i>
!
! with energy E_n and angular momentum L.
!
!

MODULE atom_2e
  !
  USE PRECISION, ONLY: dpk
  !
  IMPLICIT NONE


  TYPE  symmetry_LS2e
     ! 
     INTEGER                                :: l
     INTEGER                                :: s     
     !TD
     INTEGER                                :: net      !  brief nof of states included in the TDSE calculations
     ! note:  not always equal to 'nev' as we can restrict TDSE calculations below an E_CUT threshold     
     REAL(dpk),    ALLOCATABLE, DIMENSION(:):: e2e      !  energies in symmetry L
     COMPLEX(dpk), ALLOCATABLE, DIMENSION(:):: ct       !  time-dependent coefficients
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pt       !  brief time-dependent probability pt = ctr^2 + crti^2

     !CI
     INTEGER                                :: nev      !  nof of states included in the CI calculations
     INTEGER                                :: ncf      !  nof of configurations  included in L2e
     INTEGER                                :: nch      !  nof of channels (nll') included in L2e
     !
     INTEGER, ALLOCATABLE, DIMENSION(:)     :: n1       !  cfg electron '1'  radial q.n
     INTEGER, ALLOCATABLE, DIMENSION(:)     :: l1       !      electron '1' angular q.n
     INTEGER, ALLOCATABLE, DIMENSION(:)     :: n2       !  cfg electron '2'  radial q.n
     INTEGER, ALLOCATABLE, DIMENSION(:)     :: l2       !      electron '2' angular q.n
     INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: l12      !  existing(1) or non-existing (0) combination of l1,l2
     !
     INTEGER                                :: chl_n   !  channels icl = 1,...  = channel = (l1,l2)
     INTEGER, ALLOCATABLE, DIMENSION(:)     :: chl_1    !  l1 of ch(ich) 
     INTEGER, ALLOCATABLE, DIMENSION(:)     :: chl_2    !  l2 of ch(ich)
     !
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: ev       !  Unperturbed energies  ev = e1 + e2
     REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: cv       !  CI  coefficients
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt_0    !  brief time-dependent probability pt = ctr^2 + crti^2
     
!     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt      !>  \f$ P(ic) =  \sum_{ie} | cv(ie,ic) * ct(ie) |^2 \f$     
!     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt_b    !  
!     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt_s    !   
!     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt_d    !

     !
     REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: pop_c2e      ! population in cfg of  L

     !
     !                                    E > 0  &  e1>0,e2>0
     !                                    
     REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:) :: popd_ee_Lll  ! DI (icl,k1,k2,) for icl -> (l1,l2) channels
     REAL(dpk), ALLOCATABLE, DIMENSION(:,:)   :: popd_e1_Lll  ! DI (icl,k1)     for icl -> (l1,l2) channels 
     REAL(dpk), ALLOCATABLE, DIMENSION(:,:)   :: popd_e2_Lll  ! DI (icl,k2)     for icl -> (l1,l2) channels
     REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: popd_Lll     ! DI (icl)        for icl -> (l1,l2) channels
     REAL(dpk)                                :: popd_L       ! DI,
     !
     !                                    E+ < E < 0
     !
     REAL(dpk), ALLOCATABLE, DIMENSION(:,:)   :: pops_e_Ll    ! SI (icl), E<0   for icl -> (l1,l2) channels
     !
     !                                    E < E+ < 0
     !
     REAL(dpk), ALLOCATABLE, DIMENSION(:,:)   :: popb_e_Ll    ! SI (icl), E<A+  for icl -> (l1,l2) channels 
     !
  END TYPE symmetry_LS2e



  
CONTAINS

  !     
  ! 1-e configurations included     ! ev(n),      n = 1, 2,  ..., N_L
  ! cv(ie,ic)   i = 1, 2, ...., N_L
  !
  !

  SUBROUTINE read_w2e_ci(this,li)   ! Coulombic interaction
    !
    IMPLICIT NONE
    !
    TYPE(symmetry_ls2e), INTENT(inout) :: this
    INTEGER                            :: li
    !CFG
    INTEGER                            :: l, ls
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nhf, lhf, ll          !ncs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nllmin, nllmax, nd    !ncs nd_2 = nl2_max - nl2_min + 1 (same n1,l1,l2) 
    INTEGER                            :: l2e, n2e 
    INTEGER                            :: nhx, ncs
    !
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ch_tmp, ch_l1_tmp, ch_l2_tmp  
    
    !
    INTEGER                            :: ic      !index for total nof (k1,l1; k2,l2)  included
    INTEGER                            :: icl,ncl !index and total nof (l1,l2)         included
    INTEGER                            :: ie      !index for the k1 or k2              included
    !    
    INTEGER                            :: l1_max, l2_max, l1e_max
    INTEGER                            :: it1, it2, it3, it4, k1,  k2
    !I/O    
    INTEGER                            :: nwf2e

    ! exe
    !

    nwf2e = 9
    
    !    loop_L_symmetries: DO li = 0, lmax
       
  
       WRITE(*,'(a1,a60,G15.3,a2)')'&',' open file ' 
       CALL hfile(nwf2e,"dat","w2eb","bin",li)                       !L symmetry

        
       READ(unit=nwf2e)    l, ls
  
       IF(li.NE.l) THEN
          WRITE(*,*) '& incosistent files'
          WRITE(*,*) '&               li = ', li
          WRITE(*,*) '&                l = ', l
          STOP
       ENDIF
     
        
       !  read configurations data ( cfg-L.inp)
       
       READ(unit=nwf2e) ncs    ! nof channels (n1,l1,l2;n2)
        
       ALLOCATE(    nhf(ncs) ) 
       ALLOCATE(    lhf(ncs) )
       ALLOCATE(     ll(ncs) )
       ALLOCATE( nllmin(ncs) ) 
       ALLOCATE( nllmax(ncs) )
       ALLOCATE(     nd(ncs) )
        
     !f90
       READ(nwf2e)  nhf
       READ(nwf2e)  lhf
       READ(nwf2e)  ll
       READ(nwf2e)  nllmin
       READ(nwf2e)  nllmax
       READ(nwf2e)  nd
       READ(nwf2e)  nhx                ! nof configurations (n1,l1;n2,l2)
       READ(nwf2e)  n2e

       
       
       ! assign wf2e(l) now:
       !
        this%l   = l
        this%s   = ls
        this%nch = ncs   !channels
        this%ncf = nhx   !cfg ..
        this%nev = n2e
       


       ! read energies/coefficients from CI calculation


       
        ALLOCATE( this%ev(1:n2e)  )          ! ev = ev1 + ev2 (non-interacting e-)
        ALLOCATE( this%cv(1:n2e,1:nhx))     !'n' energy eigenstate index, 'c' --> configuration index

        read_energies_ci_2e:DO ie = 1, n2e
           READ(nwf2e) this%ev(ie)                      
           READ(nwf2e) this%cv(ie, :)    ! v(ie,jc)
        ENDDO read_energies_ci_2e
        
        ! re_index cfg(s) according to ie,ic   (n1,l1,n2,l2)


        ALLOCATE( this%n1(1:nhx)   )    !> orbital '1'
        ALLOCATE( this%l1(1:nhx)   )   
        
        ALLOCATE( this%n2(1:nhx)   )    !> orbital '2'
        ALLOCATE( this%l2(1:nhx)   )

        
        ! allocate with the (theoretically maximum possible) # of channels = ncs
        ALLOCATE(ch_l1_tmp(ncs),ch_l2_tmp(ncs)) 
        
                
        ic  = 0
        ncl = 1
        register_cfgs:DO k1 = 1, this%nch

           
           it1 = nhf(k1)
           it2 = lhf(k1)
           
           it3 =  ll(k1)
           it4 = nllmin(k1)
 
           
           DO k2 = 1, nd(k1)
              
              ic = ic + 1
              
             this%n1(ic) = it1
             this%l1(ic) = it2
             this%l2(ic) = it3
             this%n2(ic) = it4 + k2  - 1
              
          ENDDO

          !
          ! now identify the (l1,l2) partial channels for the L-symmetry
          !
          register_l_channels:IF( k1.EQ.1) THEN
             
             ch_l1_tmp(1) = it2
             ch_l2_tmp(1) = it3
          ELSE
             
             DO icl = 1, ncl          !check if the same channel has registered...
                IF( (it2.EQ.ch_l1_tmp(icl) ).AND.( it3.EQ.ch_l2_tmp(icl) ) ) CYCLE register_cfgs
                IF( (it2.EQ.ch_l2_tmp(icl) ).AND.( it3.EQ.ch_l1_tmp(icl) ) ) CYCLE register_cfgs
             ENDDO
             !             PRINT*, " pass "
             
             ncl = ncl + 1
              
             ch_l1_tmp(icl) = it2
             ch_l2_tmp(icl) = it3
             
          END IF register_l_channels
          !
       ENDDO register_cfgs


       
        ! total number of configurations included ic

       PRINT*, " total # of configurations included  ic = ", ic
       PRINT*, " total # of channels       included nch = ", this%nch
       PRINT*, " total # of (l1,l2)        included ncl = ", ncl
       PRINT*, " L = ", l

        
        ! ready to assign the channels to the module structure

       ALLOCATE(this%chl_1(1:ncl))
       ALLOCATE(this%chl_2(1:ncl))

       this%chl_n  = ncl
       this%chl_1(1:ncl) = ch_l1_tmp(1:ncl)
       this%chl_2(1:ncl) = ch_l2_tmp(1:ncl)
        

       !
       PRINT*, "(l1,l2) channels in symmetry L"
       PRINT*, this%chl_n
       PRINT*, this%chl_1
       PRINT*, this%chl_2
       !
        
       l1_max  = MAXVAL(this%l1)
       l2_max  = MAXVAL(this%l2)
       l1e_max = MAX(l1_max,l2_max)
       
       !obsolete now?
       ALLOCATE( this%l12(0:l1e_max,0:l1e_max)  )    !> channels present
       this%l12 = 0
       DO ic = 1, SIZE(this%l1)           
          this%l12(this%l1(ic),this%l2(ic)) = 1
          this%l12(this%l2(ic),this%l1(ic)) = 1
       ENDDO
       

       DO k1 = 0, l1e_max
          DO k2 = 0, l1e_max
             PRINT*, "l1_max, l2_max, k1,k2,l12(k1,k2) = ", k1, k2, this%l12(k1,k2), l1e_max
          ENDDO
       END DO
       
        
       !if l12(l1,l2) = 1 then the combination (l1,l2) is present!)
        
        this%l   = l
        this%s   = ls
        this%nch = ncs   !channels
        this%ncf = nhx   !cfg ..
        this%nev = n2e
        
        !
        !
        
        DEALLOCATE(nhf,lhf,ll,nllmin,nllmax,nd)
        DEALLOCATE(ch_l1_tmp,ch_l2_tmp)        
        
        !    END DO loop_L_symmetries
  
     RETURN
     
   END SUBROUTINE read_w2e_ci
   !
   !
   !
   !
   !
   !
   SUBROUTINE  read_w2e_ct(this, lmax)
     !
     USE netcdf
     USE ncFile_new
     USE units,          ONLY: zi
     
     !
     IMPLICIT NONE     
     !
     TYPE(symmetry_ls2e), ALLOCATABLE, DIMENSION(:), INTENT(inout) :: this                                         
     INTEGER,                                        INTENT(inout) :: lmax
     !CFG
     REAL(dpk)                              :: yr,yi
     !netCDF
     INTEGER                                :: dim_1, dim_2
     CHARACTER( LEN = 100 )                 :: ncFilename
     INTEGER,   ALLOCATABLE, DIMENSION(:)   :: iv_nc
     REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: rm_nc
     !netCDF
     INTEGER,   ALLOCATABLE, DIMENSION(:)   :: lstart, lend
     INTEGER                                :: l, ls
     INTEGER                                :: ie
     !I/O        

       
     ncFilename="tdat/hev.nc"!//atomic_name//TRIM(gauge)//ncFile_end

     !CALL create_ncFile(ncFilename)          ! create file     
     !CALL read_rs_ncFile(ncFilename,    i0,  "intensity")           !read 'intensity'
     !CALL read_rs_ncFile(ncFilename, omega,  "frequency")           !read 'frequency'
     !CALL read_rs_ncFile(ncFilename,   tau,  "duration")            !read 'duration'

     CALL read_is_ncFile(ncFilename,  lmax, "maximum_partial_wave") 
     ALLOCATE(this(0:lmax) )    
     CALL read_iv_dim_ncFile(ncFilename, iv_nc, "energy_index","l_dim")               !
     CALL read_rm_dim_ncFile(ncFilename, rm_nc, "energy_matrix","e_dim_1","e_dim_2")  !read energy matrix

  
     dim_1 = SIZE(rm_nc,dim=1)
     dim_2 = SIZE(rm_nc,dim=2)

     PRINT*, "energy_matrix:  dim_1 =  ", dim_1
     PRINT*, "energy_matrix:  dim_2 =  ", dim_2
     
 
     assign_energies_of_block_L:DO l = 0, lmax

        this(l)%net = iv_nc(l+1)
        
        ALLOCATE( this(l)%e2e( 1:this(l)%net ) )
        
        this(l)%e2e(:) = rm_nc(:,l+1)
        
     ENDDO assign_energies_of_block_L
     
     DEALLOCATE(rm_nc)                               !de-allocate to re-use
     
     
     !
     !sof: read coefficients for each 'L':  w2e(l)%ct(ie) 
     !
     
     
     ! first calculate indices for the block symmetries 'L' 
     
     ALLOCATE( lstart( 0:lmax) )
     ALLOCATE(   lend( 0:lmax) )          
     lstart(0) = 1
     DO l = 1, lmax
        lstart(l) = lstart(l-1) + this(l)%net 
     END DO     
     lend(0:lmax) = lstart(0:lmax) + this(0:lmax)%net - 1
     
     
     
     ! read real part of the coefficients array
     
     
     CALL read_rm_dim_ncFile(ncFilename, rm_nc, "re","re_dim_1","re_dim_2") 
     
     dim_1 = SIZE(rm_nc,dim=1)
     dim_2 = SIZE(rm_nc,dim=2)
     
     PRINT*, "dim_1 =  ", dim_1
     PRINT*, "dim_2 =  ", dim_2
     
     DO l = 0, lmax        
        ALLOCATE( this(l)%ct( 1:this(l)%net ) )
        ALLOCATE( this(l)%pt( 1:this(l)%net ) )
     ENDDO
        
     DO l = 0, lmax        
        re_part_of_ct:DO ie = 1, this(l)%net
           
           yr = rm_nc( dim_1, ie + lstart(l) - 1)  ! assign the real part of ct(ie)
           
           this(l)%ct(ie) = yr
           this(l)%pt(ie) = yr*yr
           
        ENDDO re_part_of_ct
        
        
     ENDDO
     
     DEALLOCATE(rm_nc)   !de-allocate to re-use
     

     !
     CALL read_rm_dim_ncFile(ncFilename, rm_nc, "im","im_dim_1","im_dim_2")  !> read imag part
     
     DO l = 0, lmax
        im_part_of_ct:DO ie = 1, this(l)%net
           
           yi =  rm_nc( dim_1, ie + lstart(l) - 1)
           
           this(l)%ct(ie) = this(l)%ct(ie) + zi * yi  ! add now the imag part of ct(ie)
           this(l)%pt(ie) = this(l)%pt(ie) + yi*yi
           
        ENDDO im_part_of_ct
     ENDDO
     
     
     !     PRINT*, i0,omega,tau
     
     
     
     DEALLOCATE(iv_nc)
     DEALLOCATE(rm_nc)
     
     DO l = 0 , lmax
        WRITE(*,'(a1,a60)') '&','................................'
        WRITE(*,'(a1,a60, I6)') '&','    angular symmetry block     l2e = ', l
        WRITE(*,'(a1,a60, I6)') '&',' nof 2e states at l2e         ne_l = ', this(l)%net 
        WRITE(*,'(a1,a60, I6)') '&','                            lstart = ', lstart(l)
        WRITE(*,'(a1,a60, I6)') '&','                              lend = ', lend(l) 
     ENDDO
     
         
     !check
     !  DO  l = 0, lmax
     !PRINT*,"---------------- l = ", l
     !DO ie = 1, 10
     !   !        PRINT*,ie, w2e(l)%e2e(ie), cct(ie+lstart(l)-1), cct(ie + lstart(l) - 1 + ntot)
     !   PRINT*,ie, this(l)%e2e(ie), REAL(this(l)%ct(ie)), AIMAG(this(l)%ct(ie))
     !ENDDO
     !ENDDO
  

  RETURN
  
END SUBROUTINE read_w2e_ct
     
   
  
END MODULE atom_2e
!
!
!
PROGRAM tdse_pop_sdi_pp_version !only the l1=1,l2=1 channels are fully analyzed
  !
  USE PRECISION
  USE parameter_tdse_fxd
  USE bs_frb_2e,      ONLY: read_target_energies
  USE io
  USE atom_2e!, ONLY: read_w2e_ci, read_w2e_ct
  !  USE atom_1e!, ONLY: read_w1e  ! to be implemented instead of read_target_energies

  IMPLICIT NONE
  !
  !
  INTEGER                                     :: i, j
  INTEGER                                     :: k, k1, k2
  REAL(dpk)                                   :: itmp
  INTEGER                                     :: li, l
  INTEGER                                     :: ie, ie_min, ie_max
  INTEGER                                     :: ic, ic_min, ic_max
  INTEGER                                     :: icl
  ! coefficients
  !field! 

  !
  !grid!

  !field-free atomic configurations
  TYPE(symmetry_ls2e), ALLOCATABLE, DIMENSION(:) :: w2e

  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: nhf, lhf, ll          !ncs
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: nd    !ncs nd_2 = nl2_max - nl2_min + 1 (same n1,l1,l2) 
  INTEGER                                     :: l2e, n2e, ls
  INTEGER                                     :: l1_max, l2_max,l1e_max 

  !
  REAL(dpk),   DIMENSION(:,:), POINTER        :: en1e                !nl,ns
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: e1e                 !nl,ns
  !  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: dk1e                !nl,ns
  !  REAL(dpk)                                   :: dk_min, factor
  !  INTEGER                                     :: ne_grid
     
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: ek1                 !nl,ns
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: rho_e1e             !nl,ns
  INTEGER                                     :: ne_1e, nl_1e
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:)    :: popl_e1e2          ! pop_e1e2(L2e, e1,e2) 
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: pop_e1e2           ! pop_e1e2(e1,e2)
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:,:,:):: pop_d_e12_channels_0 !(L2e, le1,le2,e1,e2)
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:,:)  :: pop_d_e1_channels_0  !(L2e, le1,le2,e1)
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:)    :: pop_d_channels_0     !(L2e, le1,le2,e1) 
  REAL(dpk)                                   :: sum_this, sum_this_1, sum_this_2
  REAL(dpk)                                   :: sum_d_ss, sum_d_pp, sum_d_dd, sum_d_ff
  REAL(dpk)                                   :: sum_d_sp, sum_d_pd, sum_d_df
  REAL(dpk)                                   :: sum_d_sd, sum_d_pf, sum_d_sf
  REAL(dpk)                                   :: pe
  COMPLEX(dpk)                                :: zsum
  REAL(dpk)                                   :: pop_b_0, pop_s_0, pop_d_0
  REAL(dpk)                                   :: e1,e2               !e1 = e1e(l1,n1), e2 = e1e(l2,n2)
  INTEGER                                     :: le1,le2
  !
  INTEGER                                     :: ne_1e_cut       !ne_1e_cut = ne_1e - nof_states_to_exclude
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: popl_e1,popl_e2 !
  REAL(dpk)                                   :: norm_0, norm_v  ! all population
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: popl_0, popl_v  !population in L-channel
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_di_0
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: popl_b_v, popl_b_0        !
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: popl_s_v, popl_s_0        ! E+  < E < E++
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: popl_d_v, popl_d_0        ! E++ < E
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_d_ss, pop_d_pp, pop_d_dd, pop_d_ff
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_d_sp, pop_d_pd, pop_d_df
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_d_sd, pop_d_pf
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_d_sf             !L=3
  !var
  REAL(dpk)                                   :: zz              !
  INTEGER                                     :: n_xyz           !         
  !
  !
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: ndi_l           ! pop with  E(ndi_l,l) = E++
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: nsi_l           ! pop with  E(ndi_l,l) = E+
  INTEGER                                     :: n_l2e           ! nof CI states of symmetry l2e
  !I/O

  CHARACTER(LEN=25)                           :: pop_file         !i
  CHARACTER(LEN=25)                           :: ion_file         !o
  CHARACTER(LEN=25)                           :: inp_file         !i
  CHARACTER(len = 60 )                        :: data_file
  CHARACTER(len = 6  )                        :: sl,sl1,sl2

  !
  !  CHARACTER(LEN=25)                           :: filetype         !f90,pes,cdf
  !
  
  !f90
  !INTEGER                                     :: ncoe  
  !LOGICAL                                     :: coe_file_exists
  !INTEGER                                     :: it_index, nt_out
  !INTEGER                                     :: nof_1e_states, nof_r_points
  ! read input partial pes file

   !
  !  CALL getarg(1, argv)             ! nof partial waves
  !  READ(argv,*)   lmax              !
  
  ! define i/o id's and file names


  !  ncoe  = 12
  
  !
  ! now read
  !
  
  pop_file = "tdat/pop_d_"
  ion_file = 'tdat/pop_di.dat'
  !  coe_file = 'tdat/coe.dat'

  !
  
  CALL input_tdse_fxd              ! tinp/tdse_bs_fxd_2e.inp
  CALL output_tdse_fxd             ! tout/tdse_bs_fxd_2e.out
  CALL read_target_energies(en1e)  ! read 1e energies (in Ryd)
  CALL read_w2e_ct(w2e,lmax)       ! read TDSE wavefunction $\Psi(t) = \sum_i c_{i}(t) \Phi_{i}$
  DO l = 0, lmax
     CALL read_w2e_ci(w2e(l),l)    !> read coefs $\Phi_j = \sum_i v_{ij} \Phi^{(0)}_{i}$
  ENDDO

  
  
  !
  WRITE(*,'(a1)')'&'
  WRITE(*,'(a1,a60,G15.3)')'&',' input w1e, w2e read finished ' 
  WRITE(*,'(a1)')'&'
  WRITE(*,'(a1)')'&'

  !

  
  ! some checks
  !  filetype = "cdf"
 
         
     !check

  !      PRINT*,ie, w2e(l)%e2e(ie), REAL(w2e(l)%ct(ie)), AIMAG(w2e(l)%ct(ie))
  !  DO  l = 0, lmax
  !   PRINT*,"---------------- l = ", l
  !   DO ie = 1, 10
  !      PRINT*,ie, w2e(l)%e2e(ie), REAL(w2e(l)%ct(ie)), AIMAG(w2e(l)%ct(ie))
  !   ENDDO
  !ENDDO
  
     
        
        
     DO l = 0, lmax
        WRITE(*,'(a1,a60)') '&','----------------------------------'
        WRITE(*,'(a1,a60,I6)') '&','        ang.  symmetry      l2e = ', w2e(l)%l
        WRITE(*,'(a1,a60,I6)') '&','        spin  symmetry      s2e = ', w2e(l)%s
        WRITE(*,'(a1,a60,I6)') '&',' nof zero-2e-states         n2e = ', w2e(l)%nev
        WRITE(*,'(a1,a60,I6)') '&',' nof config. series         ncs = ', w2e(l)%nch
        WRITE(*,'(a1,a60)'   ) '&',' read configuration finished.'
        !!%        
     ENDDO

     !
     !  data processing now


     !
     ! work out the 1-e energies
     !          
     ! l1e_max : is the max angular channel of the 1e-states included in the cfg
     ! nl_1e-1 : the nof angular channels calculated for the 1e-states
     ! In a full basis l1e_max =lmax_1e but in a truncated basis might differ. 
     !
     !

     !
     !NOTE::    en1e pointer index from 1:nl_1e NOT 0:nl_1e
     !
     
     en1e = en1e * 0.5_dpk              ! convert energies in a.u.
     
     
     nl_1e = SIZE(en1e,dim=1)  ! e1e(l1_max, n1_max) 
     ne_1e = SIZE(en1e,dim=2)  ! e1e(l1_max, n1_max)

     
     
     ! for conformity reasons we rewrite
     ! en1e(1:nl_1e,1:ne_1e) --> e1e(0:nl_1e-1,1:ne_1e) -->



     !First find the largest angular momentum included in the cfg channels.

     l1e_max = 0
     DO l = 0, lmax   
        IF(SIZE(w2e(l)%l12,dim=1).GT.l1e_max) l1e_max = SIZE(w2e(l)%l12,dim=1)
     ENDDO
     l1e_max = l1e_max - 1   ! since l1e_max inside the loop returns the size of l12
     
     PRINT*,"&             l1e_max = ", l1e_max
     
     IF((nl_1e-1).LT.l1e_max) THEN
        PRINT*,"something wrong here:"
        PRINT*," max l1e calculated              ne_1e = ", nl_1e-1
        PRINT*," max l1e included in cfg files l1e_max = ", l1e_max
        STOP
     ENDIF

     
     ALLOCATE( e1e(0:l1e_max,1:ne_1e))

     
     DO l = 0, l1e_max
        DO ie = 1, ne_1e
           e1e(l,ie) = en1e(l+1,ie)
        ENDDO
     ENDDO
     DEALLOCATE(en1e)
     
!!%     

  ALLOCATE(rho_e1e(0:l1e_max, 1:ne_1e) )  ! e1e(l1_max, n1_max)

  energy_density_1e:DO l = 0, l1e_max
     rho_e1e(l,1) =  1.0_dpk / ABS(e1e(l,2)-e1e(l,1))
     DO ie = 2, ne_1e-1
        rho_e1e(l,ie) = 2.0_dpk / ABS(e1e(l,ie-1)-e1e(l,ie+1))
     ENDDO
     rho_e1e(l,ne_1e) = 1.0_dpk / ABS(e1e(l,ne_1e-1)-e1e(l,ne_1e))
     
  ENDDO energy_density_1e
  
     
  !
  !
  ! 2-electron work out now
  !
  !


     
     DO l = 0, lmax
        w2e(l)%e2e = 0.5_dpk * w2e(l)%e2e    !transform to a.u. ?  
     ENDDO

     
     ALLOCATE( ndi_l(0:lmax) )
     ALLOCATE( nsi_l(0:lmax) )


     ndi_l =  0
     nsi_l  =  0
     find_si_di_thresholds:DO l = 0, lmax
        
        n_l2e = w2e(l)%net

        si:DO  ie = 1, n_l2e     ! find SI threshold energy  for symmetry L        
           nsi_l(l) = ie
           !        IF(  w2e(l)%e2e(ie) > 0.0_dpk ) EXIT
           IF(  w2e(l)%e2e(ie) > en_ion_1 ) EXIT        
        ENDDO si
        
     
     di:DO  ie = 1, n_l2e     ! find DI threshold energy  for symmetry L       
        ndi_l(l) = ie       

        ! since 2-electron systemd by definition en_ion_2 = 0.0_dpk
        
        IF( ( w2e(l)%e2e(ie) ) > 0.0_dpk ) EXIT            
     ENDDO di
     WRITE(*,*) "tdse_sdi_pop:: si threshold              en_ion_1 = ", en_ion_1
     WRITE(*,*) "tdse_sdi_pop:: si threshold  E(", nsi_l(l), l, ") = ", w2e(l)%e2e(nsi_l(l))
     WRITE(*,*) "tdse_sdi_pop:: di threshold  E(", ndi_l(l), l, ") = ", w2e(l)%e2e(ndi_l(l))
     WRITE(*,*) "tdse_sdi_pop:: max energy    E(", n_l2e, l,    ") = ", w2e(l)%e2e(n_l2e)

     !     IF(nsi_l(l).EQ.1) THEN
     !        WRITE(*,*) "tdse_sdi_pop:: ie for SI threshold can't be ", nsi_l(l)
     !        STOP
     !     ENDIF
     
  ENDDO find_si_di_thresholds

  


  ! population based on  Eq.(11) of Ref [1]   (Laulan and Bachau)
  

  ALLOCATE(popl_0(0:lmax) )
  
  popl_0 = 0.0_dpk


  partial_probability_Lv:DO  l = 0, lmax

     PRINT*, '#          l2e = ',l
     PRINT*, '#          ncf = ',w2e(l)%ncf


     
     !     ie_min = nsi_l(l)
     ie_min = ndi_l(l) 
     WRITE(*,*) "tdse_sdi_pop::       ie_min(",l,") = ", ie_min
     

    ! ie_min = 1    ;
    ! IF(l.EQ.0) ie_min = 2   ! do not include population of GROUND state (L=0)
                             ! exclude from projection (see Eq.(11) of Ref [1])


     ALLOCATE( w2e(l)%pvt_0(1:w2e(l)%ncf ))
     
     
     ie_max = w2e(l)%net             
     population_of_cfg_a:  DO ic = 1, w2e(l)%ncf  ! cfg index  
        
        zsum = 0.0_dpk
           DO ie = ie_min, ie_max                              ! j = (n1l1;n2l2;LS)
              zsum = zsum + w2e(l)%cv(ie,ic) * w2e(l)%ct(ie)   !v_j(t) = S_i c_i(t)v_ij
           ENDDO
        
        w2e(l)%pvt_0(ic) = ABS(zsum)**2                  !p_j(t) = |v_j(t)|^2  
        
     ENDDO population_of_cfg_a
        
     popl_0(l)  = SUM( w2e(l)%pvt_0 )   ! sum over all cfgs (n1l1;n2l2) ionization in L channel

  ENDDO partial_probability_Lv
  
  !
  PRINT*, "# proceed with partial wave probabilities now"
  !
  
  !REMARKS:
  !
   ! probabilities in terms of the zero 2e-functions
  !
  ! \[
  !   \psi({\bf r}_1,{\bf r}_2;t) = \sum_{L;nln'l'} v_{L;nln'l'}(t) \Phi_{L;nln'l'}({\bf r}_1,{\bf r}_2)
  ! \]
  !with
  ! \[
  !  v_{L;nln'l'}(t) = \sum_{N} C_{NL}(t) 
  ! \]  
  !


  !
  !
  !  probabilities in terms of  (E = e1 + e2) and (n1,l1, n2,l2)
  !
  !


  ALLOCATE(popl_b_0(0:lmax) )     !bound
  ALLOCATE(popl_s_0(0:lmax) )     !single ionization
  ALLOCATE(popl_d_0(0:lmax) )     !double ionization
  !                              ! (E> E++ and e1>0, e2>0)
  ALLOCATE(pop_d_ss(0:lmax))     ! populations for partial waves 
  ALLOCATE(pop_d_pp(0:lmax))
  ALLOCATE(pop_d_dd(0:lmax))
  ALLOCATE(pop_d_ff(0:lmax))  
  ALLOCATE(pop_d_sp(0:lmax))
  ALLOCATE(pop_d_pd(0:lmax))
  ALLOCATE(pop_d_df(0:lmax))     
  ALLOCATE(pop_d_sd(0:lmax))
  ALLOCATE(pop_d_pf(0:lmax))
  ALLOCATE(pop_d_sf(0:lmax))
  
  ALLOCATE(pop_di_0(0:lmax))
  !
  popl_b_0 = 0.0_dpk
  popl_s_0 = 0.0_dpk
  popl_d_0 = 0.0_dpk
  ! 
  !  
  pop_d_ss = 0.0_dpk
  pop_d_pp = 0.0_dpk
  pop_d_dd = 0.0_dpk
  pop_d_ff = 0.0_dpk
  pop_d_sp = 0.0_dpk
  pop_d_pd = 0.0_dpk
  pop_d_df = 0.0_dpk
  pop_d_sd = 0.0_dpk
  pop_d_pf = 0.0_dpk
  pop_d_sf = 0.0_dpk


!  pop_e2d_l12_channels(l:icl:k1,k2)


  ALLOCATE(    pop_d_e12_channels_0( 1:ne_1e, 1:ne_1e,0:l1e_max, 0:l1e_max,0:lmax) )
  ALLOCATE(     pop_d_e1_channels_0( 1:ne_1e, 0:l1e_max, 0:l1e_max, 0:lmax) )
  ALLOCATE(        pop_d_channels_0( 0:l1e_max, 0:l1e_max, 0:lmax) )
  !
  ALLOCATE(               popl_e1e2( 0:lmax,  1:ne_1e, 1:ne_1e) )
  ALLOCATE(                pop_e1e2( 1:ne_1e, 1:ne_1e) )         
  ALLOCATE(                 popl_e1( 0:lmax,  1:ne_1e) )
  ALLOCATE(                 popl_e2( 0:lmax,  1:ne_1e) )
  
  
  popl_e1e2         = 0.0_dpk
  pop_e1e2          = 0.0_dpk
  popl_e1           = 0.0_dpk
  popl_e2           = 0.0_dpk
  pop_d_e12_channels_0 = 0.0_dpk
  !
  PRINT*, "zero-order projections:"
  !


  !  pop_e1d1_l12_channels = 0.0_dpk ;
  !  pop_e1d2_l12_channels = 0.0_dpk ;
  
  angular_symmetry_2e_zero_order:DO l = 0, lmax


     
     PRINT*, "l2e = ", l


     !di E > E++
     pop_b_0 = 0.0_dpk
     pop_s_0 = 0.0_dpk
     pop_d_0 = 0.0_dpk     
     
     sum_d_ss = 0.0_dpk
     sum_d_pp = 0.0_dpk
     sum_d_dd = 0.0_dpk
     sum_d_ff = 0.0_dpk
     sum_d_sp = 0.0_dpk
     sum_d_pd = 0.0_dpk
     sum_d_df = 0.0_dpk
     sum_d_sd = 0.0_dpk
     sum_d_pf = 0.0_dpk
     sum_d_sf = 0.0_dpk
     !
     
     sum_all_configurations: DO ic = 1, w2e(l)%ncf

        !
        ! based on (L,ic) --->  e1,l1 and e2,l2
        !

        k1  = w2e(l)%n1(ic)
        le1 = w2e(l)%l1(ic)
        
        k2  = w2e(l)%n2(ic)
        le2 = w2e(l)%l2(ic)

        
        e1  = e1e( w2e(l)%l1(ic), w2e(l)%n1(ic) )
        e2  = e1e( w2e(l)%l2(ic), w2e(l)%n2(ic) )

        
        !

        ! according Eq. (11) of ref [1]
        ! 
        ! population of (L,ic) from all states excluding ONLY the GROUND state
        !
        ! pop(L,e1,e2) = |S_i c_i(t) * v_ij|^2

        
        pe   = w2e(l)%pvt_0(ic)     !population p(L:e_1,e2) = |S_i c_i(t)v_ij|^2  


        
        pick_channels:IF( (e1>0).AND.(e2>0) )  THEN     

           
           pop_d_0 = pop_d_0  + pe        !  (e1,e2 > 0)           
           

                      
           pop_d_e12_channels_0(k1,k2,le1,le2,l) = rho_e1e(le1,k1) * rho_e1e(le2,k2) * pe
           pop_d_e12_channels_0(k2,k1,le2,le1,l) = pop_d_e12_channels_0(k1,k2,le1,le2,l) 


           
           IF((le1.EQ.0).AND.(le2.EQ.0)) sum_d_ss = sum_d_ss + pe
           
           IF((le1.EQ.1).AND.(le2.EQ.1)) THEN
              
              sum_d_pp = sum_d_pp + pe

              !channel 1
              popl_e1e2(l,k1,k2) = rho_e1e(le1,k1) * rho_e1e(le2,k2) * pe
              popl_e1e2(l,k2,k1) = popl_e1e2(l,k1,k2)

              
           ENDIF
           
           IF((le1.EQ.2).AND.(le2.EQ.2)) sum_d_dd = sum_d_dd + pe
           IF((le1.EQ.3).AND.(le2.EQ.3)) sum_d_ff = sum_d_ff + pe

           IF((le1.EQ.0).AND.(le2.EQ.1)) sum_d_sp = sum_d_sp + pe
           IF((le1.EQ.1).AND.(le2.EQ.0)) sum_d_sp = sum_d_sp + pe
           
           IF((le1.EQ.1).AND.(le2.EQ.2)) sum_d_pd = sum_d_pd + pe
           IF((le1.EQ.2).AND.(le2.EQ.1)) sum_d_pd = sum_d_pd + pe

           IF((le1.EQ.2).AND.(le2.EQ.3)) sum_d_df = sum_d_df + pe
           IF((le1.EQ.3).AND.(le2.EQ.2)) sum_d_df = sum_d_df + pe

           IF((le1.EQ.0).AND.(le2.EQ.2)) sum_d_sd = sum_d_sd + pe
           IF((le1.EQ.2).AND.(le2.EQ.0)) sum_d_sd = sum_d_sd + pe

           IF((le1.EQ.1).AND.(le2.EQ.3)) sum_d_pf = sum_d_pf + pe
           IF((le1.EQ.3).AND.(le2.EQ.1)) sum_d_pf = sum_d_pf + pe
           
           IF((le1.EQ.0).AND.(le2.EQ.3)) sum_d_sf = sum_d_sf + pe
           IF((le1.EQ.3).AND.(le2.EQ.0)) sum_d_sf = sum_d_sf + pe

           
        ELSE IF ( (e1<0).and.(e2<0) ) THEN

           pop_b_0 = pop_b_0 + pe        ! (e1,e2 < 0)              
        ELSE
           
           pop_s_0 = pop_s_0 + pe       !  (e1<0 and e2>0 or e2 < 0 and e1 > 0)           
        ENDIF pick_channels

        
     ENDDO sum_all_configurations

     !      
     popl_b_0(l)   = pop_b_0  ! ionization in L-channel with e1<0 and e2<0     
     popl_s_0(l)   = pop_s_0  ! ionization in L-channel with e1<0 and e2>0 or vice versa
     popl_d_0(l)   = pop_d_0  ! ionization in L-channel with e1>0 and e2>0
     !
     pop_d_ss(l) = sum_d_ss
     pop_d_pp(l) = sum_d_pp
     pop_d_dd(l) = sum_d_dd
     pop_d_ff(l) = sum_d_ff
     pop_d_sp(l) = sum_d_sp
     pop_d_pd(l) = sum_d_pd
     pop_d_df(l) = sum_d_df
     pop_d_sd(l) = sum_d_sd
     pop_d_pf(l) = sum_d_pf
     pop_d_sf(l) = sum_d_sf
     
  ENDDO angular_symmetry_2e_zero_order


  ne_1e_cut = ne_1e - 11



  !
  ! ONLY THE L1=L2=1 partial waves are treated below.
  !
  
  le2 = 1
  le1 = 1


  !
  !
  ! electron - 1
  !
  !
  popl_e1 = 0.0_dpk
  probability_e1:DO l = 0, 2
     
     electron_1:DO k1 = 1, ne_1e_cut
        
        sum_this = 0.0_dpk
        sum_over_electron_2:DO k2 = 1, ne_1e_cut
           
           IF(rho_e1e(1,k2).LT.1.0e-300_dpk) STOP
           
           sum_this = sum_this + popl_e1e2(l,k1,k2)/rho_e1e(le2,k2)
           
        ENDDO sum_over_electron_2
        
        popl_e1(l,k1) = sum_this
        
     ENDDO electron_1     
  ENDDO probability_e1
  !
  !
  ! electron - 2
  !
  !  
  popl_e2 = 0.0_dpk
  probability_e2:DO l = 0, 2
     
     electron_2:DO k2 = 1, ne_1e_cut
        
        sum_this = 0.0_dpk
        sum_over_electron_1:DO k1 = 1, ne_1e_cut
           IF(rho_e1e(1,k1).LT.1.0e-300_dpk) STOP
           sum_this = sum_this + popl_e1e2(l,k1,k2)/rho_e1e(le2,k1)
           
        ENDDO sum_over_electron_1
        
        popl_e2(l,k1) = sum_this
        
     ENDDO electron_2
  ENDDO probability_e2
  




  
  !
  !  pop_d_e12_channels_0 = dP(e1,e2)/de1de2            (for fixed l,le1,le2)
  !
  !
  !-  Marginal probability for e1  
  !  pop_d_e1_channels_0 = I_de2 [ dP(e1,e2)/de1de2]    (for fixed l,le1,le2)
  !
  !
  !-  Marginal probability for e2  
  !  pop_d_e2_channels_0 = I_de1 [ dP(e1,e2)/de1de2]    (for fixed l,le1,le2)
  !
  !
  !
  !  pop_d_channels_0 = I_de_1_de2 [ dP(e1,e2)/de1de2]  (for fixed l,le1,le2)
  !
  !
  !  S_{l,le1,le2} = pop_double_ionization (e1>0, e2>)
  !


  pop_d_e1_channels_0 = 0.0_dpk
  pop_d_channels_0    = 0.0_dpk
  
  probability_e1_new:DO l = 0, lmax
     loop_l1:DO le1 = 0, l1e_max
       loop_l2:DO le2 = 0, l1e_max

          
          sum_this_2 = 0.0_dpk
          loop_e1:DO k1 = 1, ne_1e_cut
              
             sum_this_1 = 0.0_dpk
             loop_e2:DO k2 = 1, ne_1e_cut
                IF(rho_e1e(le1,k2).EQ.0.0_dpk) STOP
                sum_this_1 = sum_this_1 + pop_d_e12_channels_0(k1,k2,le1,le2,l)/rho_e1e(le2,k2)
                 
             ENDDO loop_e2
             sum_this_2 = sum_this_2 + sum_this_1/rho_e1e(le1,k1)
             
             pop_d_e1_channels_0(k1,le1,le2,l) = sum_this_1
             
          ENDDO loop_e1
          
          pop_d_channels_0(le1,le2,l) = sum_this_2
          
       ENDDO loop_l2
    ENDDO loop_l1
    
 ENDDO probability_e1_new


 !
 ! make energy grid common to all (l1,l2)
 !

 
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   DI: uncorrelated basis projection'

  !
  !  P = P_{Ll1l2}(ek1) = dP/dek_1 
  !
  !
  

     DO le1 = 0, l1e_max
        DO le2 = 0, l1e_max
           
           DO l = 0, lmax
              
              existing_channels:IF(w2e(l)%l12(le1,le2).EQ.1) THEN !print only existing channels
    
                 WRITE(sl,'(I6)')  l
                 WRITE(sl1,'(I6)') le1
                 WRITE(sl2,'(I6)') le2
              
                 data_file = &
                      &  trim(pop_file)//TRIM(ADJUSTL(sl))//TRIM(ADJUSTL(sl1))//TRIM(ADJUSTL(sl2)) &
                      & //".dat"
                 
                 OPEN(nout, file=data_file)
                 WRITE(nout,'(A1,A30,1X,3I3)') '&','pop(l1,l2,l), SUM(pop(l1,l2,:)), SUM(pop)',l,le1,le2
                 WRITE(nout,'(A1,A30)') '&',' e1, dP/de1(l1,l2,l)' 
                 WRITE(nout,'(A1,2X, 3E20.6)') '&',       &
                      & pop_d_channels_0(le1,le2,l),      &
                      & SUM(pop_d_channels_0(le1,le2,:)), &
                      & SUM(pop_d_channels_0)
                 
                 DO k1 = 1, ne_1e_cut
                    IF(e1e(le1,k1).GT.0.0_dpk) THEN
                       WRITE(nout,'(2E20.6)') e1e(le1,k1), pop_d_e1_channels_0(k1,le1,le2,l)
                    ENDIF
                 END DO
                 
                 CLOSE(nout)
                 
              ENDIF existing_channels
              
           ENDDO
           
        ENDDO
     ENDDO
     
  
     ! dump out

     ! electron-1                P(ek_1) = dP/dek_1
     
     data_file = TRIM(pop_file)//"e1.dat"
     OPEN(nout,file=data_file)

     DO k1 = 1, ne_1e_cut
        !     WRITE(10,'(5E15.6)') e1e(1,k1), popl_e1(0:lmax,k1)!, SUM(popl_e1,dim=1)
        IF(e1e(1,k1)>0.0_dpk) THEN 
           WRITE(nout,'(15E15.6)') e1e(1,k1), popl_e1(0:lmax,k1), SUM(popl_e1(:,k1))
        ENDIF
     ENDDO
     CLOSE(nout)

     ! electron-2              P(ek_2) = dP/dek_2
     
     data_file = TRIM(pop_file)//"e2.dat"
     OPEN(nout,file=data_file)

     DO k2 = 1, ne_1e_cut
        !     WRITE(10,'(5E15.6)') e1e(1,k1), popl_e1(0:lmax,k1)!, SUM(popl_e1,dim=1)
        IF(e1e(1,k2)>0.0_dpk) THEN 
           WRITE(nout,'(15E15.6)') e1e(1,k2), popl_e2(0:lmax,k2), SUM(popl_e2(:,k1))
        ENDIF
     ENDDO
     CLOSE(nout)


     
   !!!!!!!!!!   find P(e1,e2) = d^2/dek1*dek2(ek1,ek2)  for all partial waves l1,l2

     
     DO k1 = 1, ne_1e_cut
        DO k2 = 1, ne_1e_cut
           pop_e1e2(k1,k2) = SUM( popl_e1e2(0:lmax,k1,k2) )
        ENDDO
     ENDDO

     
     !! save rho(e1,e2) as table, no (e1k,e2k) points are saved
     
     data_file = TRIM(pop_file)//"e1e2-1.dat"
     OPEN(nout,file=data_file)
     

     DO k1 = ne_1e_cut,1,-1  
        WRITE(nout,*) pop_e1e2(k1,1:ne_1e_cut)
     ENDDO
     CLOSE(nout)

     !
     !
     ! now save same data in  *.vtk format
     !
     !

     
     data_file = TRIM(pop_file)//"e1e2-2.vtk"
     OPEN(nout,file=data_file)
     
     ! first count number of (x,y,z) points
     n_xyz = 0 
     nof_e1e2_points:DO k1 = 1, ne_1e_cut
        DO k2 = 1, ne_1e_cut
           IF( (e1e(1,k1)>0).AND.(e1e(1,k2)>0) )  THEN
              n_xyz = n_xyz + 1
           ENDIF
        ENDDO
     ENDDO nof_e1e2_points
     PRINT*, 'n_xyz values are to be stored in', data_file

     ! since it is a surface plot we simply set z = 0.0


     ! grid
     zz = 0.0_dpk
     WRITE(nout,*) "#vtk DataFile Version 2.0"
     WRITE(nout,*) "rho(e1,e2) data"
     WRITE(nout,*) "ASCII"
     WRITE(nout,*) "DATASET POLYDATA"
     WRITE(nout,*) "POINTS",n_xyz, "DOUBLE"        !nof (x,y,z) grid points 
     grid_e1e2:DO k1 = 1, ne_1e_cut
        DO k2 = 1, ne_1e_cut
           IF( (e1e(1,k1)>0).AND.(e1e(1,k2)>0) )  THEN
              WRITE(nout,*) e1e(1,k1),e1e(1,k2), zz   !(x,y,z) points
           ENDIF
        ENDDO
     ENDDO grid_e1e2     
     !
     WRITE(nout,*)
     WRITE(nout,*) "SCALARS Z-dimension double"
     WRITE(nout,*) "POINT_DATA",n_xyz 
     WRITE(nout,*) "LOOKUP TABLE default"
     propability_e1e2:DO k1 = 1, ne_1e_cut      !function values on the grid
        DO k2 = 1, ne_1e_cut
           IF( (e1e(1,k1)>0).AND.(e1e(1,k2)>0) )  THEN
              WRITE(nout,*) pop_e1e2(k1,k2)          !f(x,y,z) value
           ENDIF
        ENDDO
     ENDDO propability_e1e2
     CLOSE(nout)
     
     !
     ! save for ROOT
     !
     
     data_file = TRIM(pop_file)//"e1e2-root.dat"
     OPEN(nout,file=data_file)

     save_e1e2_rho:DO k1 = 1, ne_1e_cut      !function values on the grid
        DO k2 = 1, ne_1e_cut
           IF( (e1e(1,k1)>0).AND.(e1e(1,k2)>0) )  THEN
              WRITE(nout,*) e1e(1,k1),e1e(1,k2), pop_e1e2(k1,k2)!*1.0e+18 !e1,e2, rho(e1,e2)
           ENDIF
        ENDDO
     ENDDO save_e1e2_rho
     CLOSE(nout)


     !  
  
  norm_0  = SUM(popl_0) ! sum over all L-channels  (normalization to 1-p_g)

 pop_di_0 = pop_d_ss + pop_d_pp + pop_d_dd + pop_d_ff + pop_d_sp + pop_d_pd  + pop_d_df + pop_d_sd + pop_d_pf + pop_d_sf
  
  !
  !> probabilities in terms of the CI 2e-functions  psi(t) = S_NL C_NL(t) Phi_NL(r_1,r_2)
  !
  
  ALLOCATE(   popl_v(0:lmax) )  
  ALLOCATE(   popl_b_v(0:lmax) ) !bound
  ALLOCATE(   popl_s_v(0:lmax) ) !single ionization
  ALLOCATE(   popl_d_v(0:lmax) ) !double ionization

  
  popl_b_v = 0.0_dpk
  popl_s_v = 0.0_dpk
  popl_d_v = 0.0_dpk     
  populations_relative_to_ci_energies:DO l = 0, lmax
     !
     popl_b_v(l) = SUM(w2e(l)%pt(        1: nsi_l(l)-1) )   ! E_0 < E < E+
     popl_s_v(l) = SUM(w2e(l)%pt( nsi_l(l): ndi_l(l)-1) )   ! E+  < E < E++
     popl_d_v(l) = SUM(w2e(l)%pt( ndi_l(l): w2e(l)%net) )   ! E++ < E
     !   
  ENDDO populations_relative_to_ci_energies

  popl_v  = popl_b_v + popl_s_v + popl_d_v     ! all population in L-channgel bound + single + double
  norm_v = SUM(popl_v)                   ! sum over L-channels  (normalization to 1)

  ! printout results
    
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   bound population'
  PRINT*, 'x'
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                        E  < E+     pop_b_0 = ',   popl_b_0,  SUM(popl_b_0)
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                        E  < E+     pop_b_v = ',   popl_b_v,  SUM(popl_b_v)
  PRINT*, 'x'
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   single ionization'
  PRINT*, 'x'
  WRITE(*,'(A60,2X, 1P10E12.4)')  ' (e1>0 or e2>0) and  E+  < E < E++    pop_s_0 = ', popl_s_0,  SUM(popl_s_0)
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                     E+  < E < E++    pop_s_v = ', popl_s_v,  SUM(popl_s_v)
  PRINT*, 'x'
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   double ionization'
  PRINT*, 'x'
  WRITE(*,'(A60,2X, 1P10E12.4)')  ' (e1,e2 > 0) and           E++ < E    pop_d_0 = ', popl_d_0,  SUM(popl_d_0)
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                           E++ < E    pop_d_v = ', popl_d_v,  SUM(popl_d_v)
  PRINT*, 'x'  
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   total population'
  PRINT*, 'x'
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                                       norm_0 = ', norm_0, SUM(popl_0)  !`
  WRITE(*,'(A60,2X, 1P10E12.4)')  '                                       norm_v = ', norm_v, SUM(popl_v)  !
  PRINT*, 'x'
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx '
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   double ionization'
  PRINT*, 'x'
  WRITE(*,'(A30,2X, 1P10E12.4)')  ' (e1,e2 > 0) and           E++ < E    pop_dd = ', popl_d_0,   SUM(popl_d_0)
  PRINT*, 'x          L = 0,1,...,Lmax,       sum over L'
  PRINT*, 'x'
  WRITE(*,'(A15,2X, 1P10E12.4)')  '   pop_d_ss  = ', pop_d_ss,   SUM(pop_d_ss)
  WRITE(*,'(A15,2X, 1P10E12.4)')  '   pop_d_pp  = ', pop_d_pp,   SUM(pop_d_pp)
  WRITE(*,'(A15,2X, 1P10E12.4)')  '   pop_d_dd  = ', pop_d_dd,   SUM(pop_d_dd)
  WRITE(*,'(A15,2X, 1P10E12.4)')  '   pop_d_ff  = ', pop_d_ff,   SUM(pop_d_ff)
  WRITE(*,'(A15,2X, 1P10E12.4)')  '   pop_d_sp  = ', pop_d_sp,   SUM(pop_d_sp)
  WRITE(*,'(A15,2X, 1P10E12.4)')  '   pop_d_pd  = ', pop_d_pd,   SUM(pop_d_pd)
  WRITE(*,'(A15,2X, 1P10E12.4)')  '   pop_d_df  = ', pop_d_df,   SUM(pop_d_df)
  WRITE(*,'(A15,2X, 1P10E12.4)')  '   pop_d_sd  = ', pop_d_sd,   SUM(pop_d_sd)
  WRITE(*,'(A15,2X, 1P10E12.4)')  '   pop_d_pf  = ', pop_d_pf,   SUM(pop_d_pf)
  WRITE(*,'(A15,2X, 1P10E12.4)')  '   pop_d_sf  = ', pop_d_sf,   SUM(pop_d_sf)
  WRITE(*,'(A15,2X, 1P10E12.4)')  '   pop_sum   = ', pop_di_0,   SUM(pop_di_0)
  PRINT*, 'x'  

  
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   (l1,l2) channels: 1=present, 0=absent'
  DO l =0,lmax
     DO le1 = 0, l1e_max
        DO le2 = 0, l1e_max           
           PRINT*, l,le1,le2,w2e(l)%l12(le1,le2)        
        ENDDO
     ENDDO
  ENDDO
  
  
  PRINT*, 'xxxx DI: uncorrelated basis projection'
  PRINT*, 'xxxx  l, l1, l2, pop_d(l,l1,l2), exists(1)'

  DO l = 0, lmax           
     DO le1 = 0, l1e_max
        DO le2 = 0, l1e_max
           
           IF(w2e(l)%l12(le1,le2).EQ.1) THEN
              WRITE(*,'(3I3,1PE20.6,4X,I3)') &
                   & l, le1,le2, pop_d_channels_0(le1,le2,l), &
                   & w2e(l)%l12(le1,le2)
           ENDIF
        ENDDO        
     ENDDO
  ENDDO
  


END PROGRAM tdse_pop_sdi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!EOF
