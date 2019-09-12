
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
     !
     REAL(dpk), ALLOCATABLE, DIMENSION(:,:)   :: popd_e1_Lll  ! DI (icl,k1)     for icl -> (l1,l2) channels 
     REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: popd_Lll_1   ! DI (icl), sum over k2, for icl -> (l1,l2) channels
     REAL(dpk)                                :: popd_L_1     ! DI,
     !
     REAL(dpk), ALLOCATABLE, DIMENSION(:,:)   :: popd_e2_Lll  ! DI (icl,k2)     for icl -> (l1,l2) channels
     REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: popd_Lll_2   ! DI(icl), sum over k1, for icl -> (l1,l2) channels
     REAL(dpk)                                :: popd_L_2     ! DI,     
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
    INTEGER                            :: n2e 
    INTEGER                            :: nhx, ncs
    !
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ch_l1_tmp, ch_l2_tmp  
    
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
     INTEGER                                :: l!, ls
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
