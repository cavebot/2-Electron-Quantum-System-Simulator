!

!
! [1] PRA 68,013409, Laulan and Bachau,
! 'Correlation effects in two-photon single and DOUBLE ionization of helium'
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
     INTEGER, ALLOCATABLE, DIMENSION(:)     :: n1       !  cfg electron '1'  radial q.n
     INTEGER, ALLOCATABLE, DIMENSION(:)     :: l1       !      electron '1' angular q.n
     INTEGER, ALLOCATABLE, DIMENSION(:)     :: n2       !  cfg electron '2'  radial q.n
     INTEGER, ALLOCATABLE, DIMENSION(:)     :: l2       !      electron '2' angular q.n
     INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: l12      !  existing(1) or non-existing (0) combination of l1,l2
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: ev       !  Unperturbed energies  ev = e1 + e2
     REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: cv       !  CI  coefficients
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt      !>  \f$ P(ic) =  \sum_{ie} | cv(ie,ic) * ct(ie) |^2 \f$
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt_b    !  
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt_s    !   
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt_d    !
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt_0     !  brief time-dependent probability pt = ctr^2 + crti^2
     !
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pop_c2e  ! population in cfg of  L
     
     !
  END TYPE symmetry_LS2e

CONTAINS

  !     
  ! 1-e configurations included     ! ev(n),      n = 1, 2,  ..., N_L
  ! cv(ie,ic)   i = 1, 2, ...., N_L
  !
  !

  SUBROUTINE read_w2e_ci(this,li)   ! coulombic interaction
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
    INTEGER                            :: ie, ic, l1_max, l2_max, l1e_max
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

        

        
        ic = 0
        channels_index_to_cfg_index:DO k1 = 1, this%nch
           
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
        ENDDO channels_index_to_cfg_index

        l1_max = MAXVAL(this%l1)
        l2_max = MAXVAL(this%l2)
        l1e_max = MAX(l1_max,l2_max)
        

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
     REAL(dpk)                              :: ee, y2,yr,yi
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
PROGRAM tdse_pop_sdi
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
  ! coefficients
  !field! 
  REAL(dpk)                                   :: i0,omega,tau
  !
  !grid!

  !field-free atomic configurations
  TYPE(symmetry_ls2e), ALLOCATABLE, DIMENSION(:) :: w2e

  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: nhf, lhf, ll          !ncs
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: nllmin, nllmax, nd    !ncs nd_2 = nl2_max - nl2_min + 1 (same n1,l1,l2) 
  INTEGER                                     :: l2e, n2e, ls
  INTEGER                                     :: l1_max, l2_max,l1e_max 

  !
  REAL(dpk),   DIMENSION(:,:), POINTER        :: e1e                 !nl,ns
  REAL(dpk)                                   :: sum_b, sum_s, sum_d
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:)    :: sum_d_ll, sum_d_ll_0
  REAL(dpk)                                   :: sum_d_ss, sum_d_pp, sum_d_dd, sum_d_ff
  REAL(dpk)                                   :: sum_d_sp, sum_d_pd, sum_d_df
  REAL(dpk)                                   :: sum_d_sd, sum_d_pf
  REAL(dpk)                                   :: sum_d_ss_0, sum_d_pp_0, sum_d_dd_0, sum_d_ff_0
  REAL(dpk)                                   :: sum_d_sp_0, sum_d_pd_0, sum_d_df_0
  REAL(dpk)                                   :: sum_d_sd_0, sum_d_pf_0
  REAL(dpk)                                   :: pe, pe_b, pe_s
  REAL(dpk)                                   :: pe_d, pe_d_0
  COMPLEX(dpk)                                :: zsum
  REAL(dpk)                                   :: pop_b, pop_s, pop_d
  REAL(dpk)                                   :: e1,e2               !e1 = e1e(l1,n1), e2 = e1e(l2,n2)
  INTEGER                                     :: le1,le2
  !
  REAL(dpk)                                   :: pop_all, pop_all_cv, pop_all_0, pop_all_00 
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_l2e, pop_l2e_cv
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_l2e_0, pop_l2e_00
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_d_l2e
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_bi, pop_bi_0, pop_bb        ! E < E+
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_si, pop_si_0, pop_ss        ! E+  < E < E++
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_di, pop_di_0, pop_dd        ! E++ < E
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_d_ss, pop_d_pp, pop_d_dd, pop_d_ff
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_d_sp, pop_d_pd, pop_d_df
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_d_sd, pop_d_pf
  !
  !
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: ndi_l           ! pop with  E(ndi_l,l) = E++
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: nsi_l           ! pop with  E(ndi_l,l) = E+
  INTEGER                                     :: n_l2e           ! nof CI states of symmetry l2e
  !I/O
  CHARACTER(LEN=25)                           :: coe_file         !i
  CHARACTER(LEN=25)                           :: ion_file         !o
  CHARACTER(LEN=25)                           :: inp_file         !i
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


  ion_file = 'tdat/pop_di.dat'
  coe_file = 'tdat/coe.dat'

  CALL input_tdse_fxd              ! tinp/tdse_bs_fxd_2e.inp
  CALL output_tdse_fxd             ! tout/tdse_bs_fxd_2e.out
  CALL read_target_energies(e1e)   ! read 1e energies (in Ryd)
  CALL read_w2e_ct(w2e,lmax)       ! read TDSE wavefunction $\Psi(t) = \sum_i c_{i}(t) \Phi_{i}$
  DO l = 0, lmax
     CALL read_w2e_ci(w2e(l),l)    !> read coefs $\Phi_j = \sum_i v_{ij} \Phi^{(0)}_{i}$
  ENDDO



  
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


     
     e1e = e1e * 0.5_dpk              ! convert energies in A.U.
     DO l = 0, lmax
        w2e(l)%e2e = 0.5_dpk * w2e(l)%e2e 
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
        !        IF( ( w2e(l)%e2e(ie) + en_ion_1) > 0.0_dpk ) EXIT
        IF( ( w2e(l)%e2e(ie) ) > 0.0_dpk ) EXIT        
     ENDDO di
   
     WRITE(*,*) "tdse_sdi_pop:: si threshold  E(", nsi_l(l),  ") = ", w2e(l)%e2e(nsi_l(l))
     WRITE(*,*) "tdse_sdi_pop:: di threshold  E(", ndi_l(l),  ") = ", w2e(l)%e2e(ndi_l(l))
     WRITE(*,*) "tdse_sdi_pop:: max energy    E(", n_l2e,     ") = ", w2e(l)%e2e(n_l2e)

  ENDDO find_si_di_thresholds

  
  !
  !  
  !

  ! population based on  Eq.(11) of Ref [1]   (Laulan and Bachau)
  
  ALLOCATE( pop_bi_0(0:lmax) )   ! E   < E+
  ALLOCATE( pop_si_0(0:lmax) )   ! E+  < E < E++
  ALLOCATE( pop_di_0(0:lmax) )   ! E++ < E  
  ALLOCATE(pop_l2e_0(0:lmax) )
  ALLOCATE(pop_l2e_00(0:lmax) )
  
  pop_bi_0 = 0.0_dpk
  pop_si_0 = 0.0_dpk
  pop_di_0 = 0.0_dpk
  pop_l2e_0 = 0.0_dpk
  pop_l2e_00 = 0.0_dpk
  !
  partial_probability_Lv:DO  l = 0, lmax

     PRINT*, '#          l2e = ',l
     PRINT*, '#      ncf = ',w2e(l)%ncf

     ALLOCATE( w2e(l)%pvt(  1:w2e(l)%ncf ))         ! all E
     ALLOCATE( w2e(l)%pvt_b(1:w2e(l)%ncf ))         !  E   < E+
     ALLOCATE( w2e(l)%pvt_s(1:w2e(l)%ncf ))         !  E+  < E < E++
     ALLOCATE( w2e(l)%pvt_d(1:w2e(l)%ncf ))         !  E++ < E
     ALLOCATE( w2e(l)%pvt_0(1:w2e(l)%ncf ))        

     
     w2e(l)%pvt   = 0.0_dpk
     w2e(l)%pvt_b = 0.0_dpk
     w2e(l)%pvt_s = 0.0_dpk
     w2e(l)%pvt_d = 0.0_dpk

     ie_min = 1    ;
     IF(l.EQ.0) ie_min = 2   ! do not include population of GROUND state (L=0)
                             ! exclude from projection (see Eq.(11) of Ref [1])
     
     ie_max = w2e(l)%net        
     
     population_of_cfg_a:  DO ic = 1, w2e(l)%ncf  ! cfg index  

       
        pop_b = 0.0_dpk        
        bb:DO ie = ie_min, nsi_l(l) - 1              ! E < E+
           
           pop_b = pop_b + ABS( w2e(l)%cv(ie,ic))**2 * w2e(l)%pt(ie)
        ENDDO bb
        !
        pop_s = 0.0_dpk 
        ss:DO ie = nsi_l(l), ndi_l(l)-1        ! E+ < E < E++
           
           pop_s = pop_s + ABS(w2e(l)%cv(ie,ic))**2 * w2e(l)%pt(ie)           
        ENDDO ss        
        !        
        pop_d = 0.0_dpk
        dd:DO ie = ndi_l(l), ie_max           !  E++ < E
           
           pop_d = pop_d + ABS(w2e(l)%cv(ie,ic))**2 * w2e(l)%pt(ie)                      
        ENDDO dd

        w2e(l)%pvt_b(ic) = pop_b       !so this pop_b do not include pop of ground state
        w2e(l)%pvt_s(ic) = pop_s
        w2e(l)%pvt_d(ic) = pop_d
        
        w2e(l)%pvt(ic) = pop_s + pop_d !+ pop_b
        
        zsum = 0.0_dpk
           DO ie = ie_min, ie_max                              ! j = (n1l1;n2l2;LS)
              zsum = zsum + w2e(l)%cv(ie,ic) * w2e(l)%ct(ie)   !v_j(t) = S_i c_i(t)v_ij
           ENDDO
        
        w2e(l)%pvt_0(ic) = ABS(zsum)**2                  !p_j(t) = |v_j(t)|^2  
        
     ENDDO population_of_cfg_a
        
     pop_bi_0(l)   = SUM( w2e(l)%pvt_b )    ! sum over cfgs  ! E+  > E
     pop_si_0(l)   = SUM( w2e(l)%pvt_s )                     ! E+  < E < E++ 
     pop_di_0(l)   = SUM( w2e(l)%pvt_d )                     ! E++ < E
     !
     pop_l2e_0(l)  = SUM( w2e(l)%pvt   )   ! sum over all cfgs (n1l1;n2l2) ionization in L channel
     pop_l2e_00(l) = SUM( w2e(l)%pvt_0 )   
     
  ENDDO partial_probability_Lv

  
  !
  pop_all_00 = SUM(pop_l2e_00)
  pop_all_0  = SUM(pop_l2e_0 )                   ! sum over l2e

  
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

  ALLOCATE(pop_bb(0:lmax) )     !bound
  ALLOCATE(pop_ss(0:lmax) )     !single ionization
  ALLOCATE(pop_dd(0:lmax) )     !double ionization
  ALLOCATE(pop_l2e_cv(0:lmax))
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
  ALLOCATE(pop_d_l2e(0:lmax))
  !
  pop_bb = 0.0_dpk
  pop_ss = 0.0_dpk
  pop_dd = 0.0_dpk
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

  ! find the largest angular momentum included in the cfg channels.

  l1e_max = 0
  DO l = 0, lmax   
     IF(SIZE(w2e(l)%l12,dim=1).GT.l1e_max) l1e_max = SIZE(w2e(l)%l12,dim=1)
  ENDDO
  l1e_max = l1e_max - 1   ! since l1e_max inside the loop returns the size of l12
  
  PRINT*,"l1e_max = ", l1e_max
  
  ALLOCATE( sum_d_ll  ( 0:lmax, 0:l1e_max, 0:l1e_max) )
  ALLOCATE( sum_d_ll_0( 0:lmax, 0:l1e_max, 0:l1e_max) )

       !
  sum_d_ll   = 0.0_dpk
  sum_d_ll_0 = 0.0_dpk

  PRINT*, "zero-order projections:"
  angular_symmetry_2e_zero_order:DO l = 0, lmax

     PRINT*, "l2e = ", l


     !di E > E++
     pop_b = 0.0_dpk
     pop_s = 0.0_dpk
     pop_d = 0.0_dpk     
     sum_b = 0.0_dpk
     sum_s = 0.0_dpk
     sum_d = 0.0_dpk
     
     sum_d_ss = 0.0_dpk
     sum_d_pp = 0.0_dpk
     sum_d_dd = 0.0_dpk
     sum_d_ff = 0.0_dpk
     sum_d_sp = 0.0_dpk
     sum_d_pd = 0.0_dpk
     sum_d_df = 0.0_dpk
     sum_d_sd = 0.0_dpk
     sum_d_pf = 0.0_dpk
     !
     sum_d_ss = 0.0_dpk
     sum_d_pp = 0.0_dpk
     sum_d_dd = 0.0_dpk
     sum_d_ff = 0.0_dpk
     sum_d_sp = 0.0_dpk
     sum_d_pd = 0.0_dpk
     sum_d_df = 0.0_dpk
     sum_d_sd = 0.0_dpk
     sum_d_pf = 0.0_dpk

     
     sum_all_configurations: DO ic = 1, w2e(l)%ncf

        !
        ! based on (L,ic) --->  e1,l1 and e2,l2
        !
        
        e1  = e1e( w2e(l)%l1(ic), w2e(l)%n1(ic) )
        e2  = e1e( w2e(l)%l2(ic), w2e(l)%n2(ic) )
        le1 = w2e(l)%l1(ic)
        le2 = w2e(l)%l2(ic)

        !
        ! population of (L,ic) from all states excluding ALL bound states
        !
        ! pop(L,j) = S_i |c_i(t)|^2 |v_ij|^2
        !
        
        pe   = w2e(l)%pvt(ic)


        ! 
        ! population of (L,ic) from all states excluding ONLY the GROUND state
        !
        ! pop(L,j) = |S_i c_i(t) * v_ij|^2
        !
        ! most likely it should be
        !
        !        pe   = w2e(l)%pvt_0(ic) 
        


        
        pick_channels:IF( (e1>0).AND.(e2>0) )  THEN
           
           pe_d =   w2e(l)%pvt_d(ic)  
                      
           pop_d = pop_d  + pe       !  (e1,e2 > 0)           
           sum_d = sum_d +  pe_d     !  (e1,e2 > 0) and (E > E++)

           sum_d_ll(l, le1, le2) =   sum_d_ll(l, le1, le2 ) + pe_d


           
           ! according Eq. (11) of ref [1]
           
           pe_d_0 = w2e(l)%pvt_0(ic)  !population p(j: e_1,e2>0) = |S_i c_i(t)v_ij|^2  
           sum_d_ll_0(l, le1, le2) = sum_d_ll_0(l, le1, le2 ) + pe_d_0

           
           IF((le1.EQ.0).AND.(le2.EQ.0)) sum_d_ss = sum_d_ss + pe_d_0
           IF((le1.EQ.1).AND.(le2.EQ.1)) sum_d_pp = sum_d_pp + pe_d_0
           IF((le1.EQ.2).AND.(le2.EQ.2)) sum_d_dd = sum_d_dd + pe_d_0
           IF((le1.EQ.3).AND.(le2.EQ.3)) sum_d_ff = sum_d_ff + pe_d_0
           IF((le1.EQ.0).AND.(le2.EQ.1)) sum_d_sp = sum_d_sp + pe_d_0
           IF((le1.EQ.1).AND.(le2.EQ.0)) sum_d_sp = sum_d_sp + pe_d_0
           IF((le1.EQ.1).AND.(le2.EQ.2)) sum_d_pd = sum_d_pd + pe_d_0
           IF((le1.EQ.2).AND.(le2.EQ.1)) sum_d_pd = sum_d_pd + pe_d_0
           IF((le1.EQ.2).AND.(le2.EQ.3)) sum_d_df = sum_d_df + pe_d_0
           IF((le1.EQ.3).AND.(le2.EQ.2)) sum_d_df = sum_d_df + pe_d_0
           IF((le1.EQ.0).AND.(le2.EQ.2)) sum_d_sd = sum_d_sd + pe_d_0
           IF((le1.EQ.2).AND.(le2.EQ.0)) sum_d_sd = sum_d_sd + pe_d_0
           IF((le1.EQ.1).AND.(le2.EQ.3)) sum_d_pf = sum_d_pf + pe_d_0
           IF((le1.EQ.3).AND.(le2.EQ.1)) sum_d_pf = sum_d_pf + pe_d_0

           
        ELSE IF ( (e1<0).AND.(e2<0) ) THEN

           !i think it should be
           !pe_b  = w2e(l)%pvt_0(ic)

           pe_b  = w2e(l)%pvt_b(ic)
           
           pop_b = pop_b + pe        ! (e1,e2 < 0)   
           sum_b = sum_b + pe_b      ! (e1,e2 < 0)  and   (E < E+)

           
        ELSE

           ! i think it should be
           !pe_s  = w2e(l)%pvt_0(ic)
           
           pe_s  = w2e(l)%pvt_s(ic)
           
           pop_s = pop_s + pe        ! (e1<0 or e2 < 0)
           sum_s = sum_s + pe_s      ! (e1<0 or e2<0) and (E+ < E < E++)
           
        ENDIF pick_channels

        
     ENDDO sum_all_configurations


     !      
     pop_bb(l) = sum_b
     pop_ss(l) = sum_s
     pop_dd(l) = sum_d
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
     
     
     
  ENDDO angular_symmetry_2e_zero_order


  pop_l2e_cv  = pop_bb + pop_ss + pop_dd   ! sum over b,s,d
  pop_all_cv  = SUM(pop_l2e_cv)            ! sum over l2e



  pop_d_l2e = pop_d_ss + pop_d_pp + pop_d_dd + pop_d_ff + pop_d_sp + pop_d_pd  + pop_d_df + pop_d_sd + pop_d_pf 

  
  !
  !
  !> probabilities in terms of the CI 2e-functions  psi(t) = S_NL C_NL(t) Phi_NL(r_1,r_2)
  !
  !

  
  ALLOCATE(  pop_l2e(0:lmax) )  
  ALLOCATE(   pop_bi(0:lmax) ) !bound
  ALLOCATE(   pop_si(0:lmax) ) !single ionization
  ALLOCATE(   pop_di(0:lmax) ) !double ionization

  
  pop_bi = 0.0_dpk
  pop_si = 0.0_dpk
  pop_di = 0.0_dpk     
  populations_relative_to_ci_energies:DO l = 0, lmax
               
     pop_bi(l) = SUM(w2e(l)%pt(        1: nsi_l(l)-1) )   ! E_0 < E < E+
     pop_si(l) = SUM(w2e(l)%pt( nsi_l(l): ndi_l(l)-1) )   ! E+  < E < E++
     pop_di(l) = SUM(w2e(l)%pt( ndi_l(l): w2e(l)%net) )   ! E++ < E

     
  ENDDO populations_relative_to_ci_energies

  pop_l2e = pop_bi + pop_si + pop_di     ! sum over b,s,d
  pop_all = SUM(pop_l2e)                 ! sum over l2e

   
  
  !
  !
  ! printout results
  !
  !
  
  
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   bound population'
  PRINT*, 'x'
  PRINT*, '                        E  < E+     pop_bb = ',   pop_bb,  SUM(pop_bb)
  PRINT*, '                        E  < E+     pop_bi = ',   pop_bi,  SUM(pop_bi)
  PRINT*, '                        E  < E+   pop_bi_0 = ', pop_bi_0,  SUM(pop_bi_0)
  PRINT*, 'x'
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   single ionization'
  PRINT*, 'x'
  PRINT*, ' (e1>0 or e2>0) and  E+  < E < E++    pop_ss = ', pop_ss,   SUM(pop_ss)
  PRINT*, '                     E+  < E < E++    pop_si = ', pop_si,   SUM(pop_si)
  PRINT*, ' (e1>0 or e2>0)                     pop_si_0 = ', pop_si_0, SUM(pop_si_0)
  PRINT*, 'x'
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   double ionization'
  PRINT*, 'x'
  PRINT*, ' (e1,e2 > 0) and           E++ < E    pop_dd = ', pop_dd,   SUM(pop_dd)
  PRINT*, '                           E++ < E    pop_di = ', pop_di,   SUM(pop_di)
  PRINT*, ' (e1,e2 > 0) and                    pop_di_0 = ', pop_di_0, SUM(pop_di_0)
  PRINT*, 'x'  
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   total population'
  PRINT*, 'x'
  PRINT*, '  E,e1,2                            pop_cv(l) = ', pop_l2e_cv, pop_all_cv
  PRINT*, '  E                                    pop(l) = ', pop_l2e,   pop_all
  PRINT*, '  e1,e2                              pop_0(l) = ', pop_l2e_0, pop_all_0
  PRINT*, '  e1,e2  (coherent sum)             pop_00(l) = ', pop_l2e_00, pop_all_00
  PRINT*, 'x'
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx '
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   double ionization'
  PRINT*, 'x'
  PRINT*, ' (e1,e2 > 0) and           E++ < E    pop_dd = ', pop_dd,   SUM(pop_dd)
  PRINT*, 'x          L = 0,1,...,Lmax,       sum over L'
  PRINT*, '   pop_d_ss  = ', pop_d_ss,   SUM(pop_d_ss)
  PRINT*, '   pop_d_pp  = ', pop_d_pp,   SUM(pop_d_pp)
  PRINT*, '   pop_d_dd  = ', pop_d_dd,   SUM(pop_d_dd)
  PRINT*, '   pop_d_ff  = ', pop_d_ff,   SUM(pop_d_ff)
  PRINT*, '   pop_d_sp  = ', pop_d_sp,   SUM(pop_d_sp)
  PRINT*, '   pop_d_pd  = ', pop_d_pd,   SUM(pop_d_pd)
  PRINT*, '   pop_d_df  = ', pop_d_df,   SUM(pop_d_df)
  PRINT*, '   pop_d_sd  = ', pop_d_sd,   SUM(pop_d_sd)
  PRINT*, '   pop_d_pf  = ', pop_d_pf,   SUM(pop_d_pf)
  PRINT*, '   pop_sum   = ',  pop_d_l2e, SUM(pop_d_l2e)
  PRINT*, 'x'  



  
 
  
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   (l1,l2) channels: 1=present, 0=absent'
  DO l =0,lmax
     DO k1 = 0, l1e_max
        DO k2 = 0, l1e_max
           
           PRINT*, l,k1,k2,w2e(l)%l12(k1,k2)
        
        ENDDO
     ENDDO
  ENDDO
  
  
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   DI: uncorrelated basis projection'
  DO k1 = 0, l1e_max
     DO k2 = k1, l1e_max

        print_only_existing_channels_0:DO l = 0, lmax
           IF(w2e(l)%l12(k1,k2).EQ.0) THEN
              CYCLE
           ELSE
              WRITE(*,'(2I3,10E20.6)') k1,k2,(sum_d_ll_0(l2e,k1,k2),l2e=0,lmax)
              EXIT
           ENDIF
        ENDDO print_only_existing_channels_0
        
     ENDDO
  ENDDO

  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   DI: correlated basis projection'    
  DO k1 = 0, l1e_max
     DO k2 = k1, l1e_max
        
        print_only_existing_channels:DO l = 0, lmax
           IF(w2e(l)%l12(k1,k2).EQ.0) THEN
              CYCLE
           ELSE
              WRITE(*,'(2I3,10E20.6)') k1,k2,(sum_d_ll(l2e,k1,k2),l2e=0,lmax)
              EXIT
           ENDIF
        ENDDO print_only_existing_channels
        
     ENDDO
  ENDDO

END PROGRAM tdse_pop_sdi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!EOF
