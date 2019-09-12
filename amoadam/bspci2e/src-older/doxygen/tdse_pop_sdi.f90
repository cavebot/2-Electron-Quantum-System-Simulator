!
! calculates spatial density matrix, rho(r1,r2,t)
!
!> This program calculates bound (b), single-ionization (si) and double-ionization (di)
!! probabilities based on the tdse.nc output of the tdse calculation of a 2e-system
!!
!!
!!  \f$ P_{di} = \sum_{E>0,nl>0} |\langle \phi_0|\psi(t) \rangle|^2\f$
!!

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
     REAL(dpk),    ALLOCATABLE, DIMENSION(:):: ent      !  energies in symmetry L
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
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: ev       !  Unperturbed energies  ev = e1 + e2
     REAL(dpk), ALLOCATABLE, DIMENSION(:,:) :: cv       !  CI  coefficients
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt      !> @var \f$ P(ic) =  \sum_{ie} | cv(ie,ic) * ct(ie) |^2 \f$
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt_b    !  
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt_s    !   
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pvt_d    !   
     !
     REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: pop_c2e  ! population in cfg of  L
     
     !
  END TYPE symmetry_LS2e

CONTAINS

  
  SUBROUTINE init_symmetry_ls2e(this, l_2e, s_2e, n_e2e, n_c2e)
    !
    IMPLICIT NONE
    !ARG
    TYPE(symmetry_ls2e), INTENT(inout) :: this    
    INTEGER                            :: l_2e 
    INTEGER                            :: s_2e
    INTEGER                            :: n_e2e
    INTEGER                            :: n_c2e

    !loc
    !EXE!

    this%l  = l_2e           ! angular number
    this%s  = s_2e           ! spin number
    this%net = n_e2e          ! nof energy states included in l_2e
    this%nch = n_c2e          ! nof configurations included in l_2e


    !
    ALLOCATE(this%ent(this%net))
    ALLOCATE(this%pop_c2e(this%nch))

    !
  END SUBROUTINE init_symmetry_LS2e
  !     
  ! 1-e configurations included     ! ev(n),      n = 1, 2,  ..., N_L
  ! cv(ie,ic)   i = 1, 2, ...., N_L
  !
  !
  SUBROUTINE init_w2e_configurations(this,li)
    !
    IMPLICIT NONE
    !
    TYPE(symmetry_ls2e), INTENT(inout) :: this
    !CFG
    INTEGER                            :: l, li, ls
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nhf, lhf, ll          !ncs
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nllmin, nllmax, nd    !ncs nd_2 = nl2_max - nl2_min + 1 (same n1,l1,l2) 
    INTEGER                            :: l2e, n2e 
    INTEGER                            :: nhx, ncs
    !
    INTEGER                            :: ie, ic
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
       


       ! read energies/coefficients from CI calculation


       
        ALLOCATE( this%ev(1:n2e)  )          ! ev = ev1 + ev2 (non-interacting e-)
        ALLOCATE( this%cv(1:n2e,1:nhx))     !'n' energy eigenstate index, 'c' --> configuration index
     
    
        read_energies_ci_2e:DO ie = 1, n2e
           READ(nwf2e) this%ev(ie)                      
           READ(nwf2e) this%cv(ie, :)    ! v(ie,jc)
        ENDDO read_energies_ci_2e
        
        
        ! re_index cfg(s) according to ie,ic   (n1,l1,n2,l2)

        
        ALLOCATE( this%n1(1:nhx)   )
        ALLOCATE( this%l1(1:nhx)   )
        ALLOCATE( this%n2(1:nhx)   )
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



        
        this%l   = l
        this%s   = ls
        this%nch = ncs   !channels
        this%ncf = nhx   !cfg ..
        this%nev = n2e
        
        
        
        WRITE(*,'(a1,a60)') '&','----------------------------------'
        WRITE(*,'(a1,a60,I6)') '&','        ang.  symmetry      l2e = ', this%l
        WRITE(*,'(a1,a60,I6)') '&','        spin  symmetry      s2e = ', this%s
        WRITE(*,'(a1,a60,I6)') '&',' nof zero-2e-states         n2e = ', this%nev
        WRITE(*,'(a1,a60,I6)') '&',' nof config. series         ncs = ', this%nch
        WRITE(*,'(a1,a60,I6)') '&',' matrix dimension           nhx = ', nhx
        WRITE(*,'(a1,a60,I6)') '&',' matrix dimension            ic = ', ic
        WRITE(*,'(a1,a60)'   ) '&',' read configuration finished.'
        
        
!        PRINT*, w2e(l)%n1(1),    w2e(l)%l1(1),    w2e(l)%l2(1),    w2e(l)%n2(1)
!        PRINT*, w2e(l)%n1(1831), w2e(l)%l1(1831), w2e(l)%l2(1831), w2e(l)%n2(1831)
!        PRINT*, w2e(l)%n1(3661), w2e(l)%l1(3661), w2e(l)%l2(3661), w2e(l)%n2(3661)
!        PRINT*, w2e(l)%n1(5491), w2e(l)%l1(5491), w2e(l)%l2(5491), w2e(l)%n2(5491)
!        PRINT*, w2e(l)%n1(7320), w2e(l)%l1(7320), w2e(l)%l2(7320), w2e(l)%n2(7320)
        
        !
        !
        
        DEALLOCATE(nhf,lhf,ll,nllmin,nllmax,nd)
        
        
        !    END DO loop_L_symmetries
  
     RETURN
     
   END SUBROUTINE init_w2e_configurations

  
  
END MODULE atom_2e
!
!
!
PROGRAM tdse_pop_sdi
  !
  USE PRECISION
  USE units, only: zi
  USE parameter_tdse_fxd
  USE bs_frb_2e, ONLY: read_target_energies
  !
  USE io
  USE netcdf
  USE ncFile_new
  !
  USE atom_2e, ONLY: symmetry_ls2e
  USE utils,   ONLY: asciifile,datafile
  !
  IMPLICIT NONE
  !
  !
  INTEGER, PARAMETER                          :: l1e_max = 5
  INTEGER                                     :: i, j, k
  REAL(dpk)                                   ::itmp
  INTEGER                                     :: li, l
  INTEGER                                     :: ie, ic, k1, k2
  ! coefficients
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        ::  cct        ! time-dep endent coeff
  !field! 
  REAL(dpk)                                   :: i0,omega,tau
  !
  !grid!

  INTEGER                                     :: it, ir1, ir2, n_lab !
  INTEGER                                     :: la, lb, ie_min, ie_max, ill !

  !field-free atomic configurations
  INTEGER                                     :: ncs, nhx
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: lstart, lend
  TYPE(symmetry_ls2e), ALLOCATABLE, DIMENSION(:)  :: w2e

  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: nhf, lhf, ll          !ncs
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: nllmin, nllmax, nd    !ncs nd_2 = nl2_max - nl2_min + 1 (same n1,l1,l2) 
  INTEGER                                     :: l2e, n2e, ls
  INTEGER                                     :: nwf2e
  !
  INTEGER                                     :: ns,np               ! l1e_max
  REAL(dpk),   DIMENSION(:,:), POINTER        :: e1e                 !nl,ns
  INTEGER                                     :: it1, it2, it3, it4

  REAL(dpk)                                   :: sum_b, sum_s, sum_d
  REAL(dpk)                                   :: sum_d_ss, sum_d_pp, sum_d_dd
  REAL(dpk)                                   :: sum_d_sp, sum_d_pd, sum_d_df
  REAL(dpk)                                   :: sum_d_sd, sum_d_pf
  REAL(dpk)                                   :: pe, pe_b, pe_s, pe_d
  REAL(dpk)                                   :: pop_ic
  REAL(dpk)                                   :: pop_b, pop_s, pop_d
  REAL(dpk)                                   :: e1,e2               !e1 = e1e(l1,n1), e2 = e1e(l2,n2)
  INTEGER                                     :: le1,le2
  !
  REAL(dpk)                                   :: pop_all, pop_all_0, pop_all_cv      
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_l2e, pop_l2e_0, pop_l2e_cv 
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_bi, pop_bi_0, pop_bb        ! E < E+
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_si, pop_si_0, pop_ss        ! E+  < E < E++
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_di, pop_di_0, pop_dd        ! E++ < E
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: pop_d_ss, pop_d_pp, pop_d_dd
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
  CHARACTER(LEN=25)                           :: filetype         !f90,pes,cdf
  !
  !ncFile
  INTEGER                                     :: dim_1, dim_2
  CHARACTER( LEN = 100 )                      :: ncFilename
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: iv_nc
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: rm_nc
  !
  !f90
  INTEGER                                     :: ncoe  
  LOGICAL                                     :: coe_file_exists
  INTEGER                                     :: it_index, nt_out
  INTEGER                                     :: nof_1e_states, nof_r_points
  ! read input partial pes file
  INTEGER                                     :: npes
  CHARACTER(len=1)                            :: ctmp
  INTEGER                                     :: nl, ne_tot_l
  INTEGER                                     :: l_pes, ntmp
  REAL(dpk)                                   :: ee, y2,yr,yi
  !
  
  !
  !exe!
  
      !
 !  CALL getarg(1, argv)             ! nof partial waves
 !  READ(argv,*)   lmax              !

  ! define i/o id's and file names


  nwf2e = 9
  ncoe  = 12
  
  !
  ! now read
  !


  ion_file = 'tdat/pop_di.dat'
  coe_file = 'tdat/coe.dat'

  CALL input_tdse_fxd              ! tinp/tdse_bs_fxd_2e.inp
  CALL output_tdse_fxd             ! tout/tdse_bs_fxd_2e.out

  !start h1e
  
  ! read 1e-orbitals energies 
  !

  
  CALL read_target_energies(e1e)   ! read 1e energies (in Ryd)
  e1e = e1e * 0.5_dpk              ! e1e in a.u. now
    
  WRITE(*,'(a1,a60)') '&',' 1e partial waves read.'

  !end h1e

  !start h2e TD coefficients
  
  !
  !
  !  read 2e time-dependent coefficients data         C_nl(t)
  !
  ! intensity, frequency, total pulse duration, maximum_L_2e, maximum_energy_2e
  ! en(N,L)
  ! ntot
  ! cct(N,L)
  !


  filetype = "cdf"
  !filetype = "ascii"


!  IF(fileType.EQ."f90") THEN

!!%
!!%  
!!%  !check if it is present
!!%
!!%  
!!%  WRITE(*,'(a1,a60,a20)') '&',' open coefficient file:', coe_file
!!%  INQUIRE(file=coe_file,exist=coe_file_exists)  
!!%  IF(.NOT.coe_file_exists) THEN
!!%     WRITE(*,'(a1,a60)') '&',' coefficient file is not present. Aborting'
!!%     STOP
!!%  ENDIF
!!%  WRITE(*,'(a1,a60)') '&',' read from coefficient file'  
!!%
!!%  
!!%
!!%
!!%  ! Use time-dependent coefficient file 
!!%
!!%  OPEN(ncoe,file=coe_file,form='unformatted',access='sequential')
!!%
!!%  !read coefficients file and  read dimensions
!!%  
!!%  READ(ncoe)  i0, omega, tau, lmax, en_cut
!!%  ALLOCATE( w2e(0:lmax) )
!!%  DO  l = 0, lmax
!!%     READ(ncoe)  w2e(l)%net, (rtmp, ie = 1, w2e(l)%net )      !     READ(ncoe)  i, (rcoe, ie = 1, ne_l(l) )
!!%  ENDDO
!!%  READ(ncoe) ntot
!!%  CLOSE(ncoe)
!!%
!!%  !........... allocate dimensions
!!%
!!%  DO l = 0, lmax
!!%     ALLOCATE( w2e(l)%ent(1:w2e(l)%net ) )
!!%     ALLOCATE(  w2e(l)%ct(1:w2e(l)%net ) )
!!%  ENDDO
!!%  ALLOCATE(cct(1:2*ntot))
!!%  
!!%
!!%  ! read again coefficients file and load variables
!!%  
!!%  OPEN(ncoe,file=coe_file,form='unformatted',access='sequential') 
!!%  READ(ncoe)  i0, omega, tau, lmax, en_cut
!!%  read_energies_l:DO  l = 0, lmax
!!%     READ(ncoe)  i, w2e(l)%ent
!!%  ENDDO read_energies_l
!!%  READ(ncoe) ntot
!!%  READ(ncoe) time,cct
!!%  CLOSE(ncoe)

  !
  !
  !
  !  end of coefficient file
  !
  !
     !
  IF(fileType.EQ."cdf") THEN
     
     !  ELSE IF(filetype.EQ."cdf") THEN

     !locate the netCDF file = 'atomic_name' + 'gauge' + '.nc'

     ! 
     
     ncFilename="tdat/hev.nc"!//atomic_name//TRIM(gauge)//ncFile_end

     !  CALL create_ncFile(ncFilename)          ! create file
     CALL read_rs_ncFile(ncFilename,    i0,  "intensity")           !read 'intensity'
     CALL read_rs_ncFile(ncFilename, omega,  "frequency")           !read 'frequency'
     CALL read_rs_ncFile(ncFilename,   tau,  "duration")            !read 'duration'
     CALL read_is_ncFile(ncFilename,  lmax, "maximum_partial_wave") !read 'maximum_partial_wave'

     PRINT*, ' lmax = ', lmax

     
     ALLOCATE(   w2e(0:lmax) )

     
     !
     !sof: read energies for each 'L':  w2e(l)%ent(ie)
     !
     
     CALL read_iv_dim_ncFile(ncFilename, iv_nc, "energy_index","l_dim")  !read energy states for each 'L'  
     CALL read_rm_dim_ncFile(ncFilename, rm_nc, "energy_matrix","e_dim_1","e_dim_2")  !read energy states for each 'L'

  
     dim_1 = SIZE(rm_nc,dim=1)
     dim_2 = SIZE(rm_nc,dim=2)

     PRINT*, "energy_matrix:  dim_1 =  ", dim_1
     PRINT*, "energy_matrix:  dim_2 =  ", dim_2
     
 
     DO l = 0, lmax

        w2e(l)%net = iv_nc(l+1)
        
        PRINT*, ' l, ne_l = ',  l, w2e(l)%net 
        ALLOCATE( w2e(l)%ent( 1:w2e(l)%net ) )
        
        w2e(l)%ent(:) = 0.5_dpk * rm_nc(:,l+1)      ! assign energy in A.U. array for each 'L' symmetry 
        
     ENDDO
     
     DEALLOCATE(rm_nc)    !de-allocate to re-used
     
     !
     !eof: read energies
     !
     
     
     !
     !sof: read coefficients for each 'L':  w2e(l)%ct(ie) 
     !
     
     
     ! first calculate indices for the block symmetries 'L' 
     
     ALLOCATE( lstart( 0:lmax) )
     ALLOCATE(   lend( 0:lmax) )
     
     
     lstart(0) = 1
     DO l = 1, lmax
        lstart(l) = lstart(l-1) + w2e(l)%net 
     END DO
     
     lend(0:lmax) = lstart(0:lmax) + w2e(0:lmax)%net - 1
     
     
     !
     ! read real part of the coefficients array
     !
     
     CALL read_rm_dim_ncFile(ncFilename, rm_nc, "re","re_dim_1","re_dim_2") 
     
     dim_1 = SIZE(rm_nc,dim=1)
     dim_2 = SIZE(rm_nc,dim=2)
     
     PRINT*, "dim_1 =  ", dim_1
     PRINT*, "dim_2 =  ", dim_2
     
     DO l = 0, lmax
        
        ALLOCATE( w2e(l)%ct( 1:w2e(l)%net ) )
        ALLOCATE( w2e(l)%pt( 1:w2e(l)%net ) )
        
        re_part_of_ct:DO ie = 1, w2e(l)%net
           
           yr = rm_nc( dim_1, ie + lstart(l) - 1)  ! assign the real part of ct(ie)
           
           w2e(l)%ct(ie) = yr
           w2e(l)%pt(ie) = yr*yr
           
        ENDDO re_part_of_ct
        
        
     ENDDO
     
     DEALLOCATE(rm_nc)   !de-allocate to re-used
     
     !
     ! read real part of the coefficients array
     !
     CALL read_rm_dim_ncFile(ncFilename, rm_nc, "im","im_dim_1","im_dim_2")  
     
     DO l = 0, lmax
        im_part_of_ct:DO ie = 1, w2e(l)%net
           
           yi =  rm_nc( dim_1, ie + lstart(l) - 1)
           
           w2e(l)%ct(ie) = w2e(l)%ct(ie) + zi * yi  ! add now the imag part of ct(ie)
           w2e(l)%pt(ie) = w2e(l)%pt(ie) + yi*yi
           
        ENDDO im_part_of_ct
     ENDDO
     
     
     PRINT*, i0,omega,tau
     
     
     
     DEALLOCATE(iv_nc)
     DEALLOCATE(rm_nc)
     
     DO l = 0 , lmax
        WRITE(*,'(a1,a60)') '&','................................'
        WRITE(*,'(a1,a60, I6)') '&','    angular symmetry block     l2e = ', l
        WRITE(*,'(a1,a60, I6)') '&',' nof 2e states at l2e         ne_l = ', w2e(l)%net 
        WRITE(*,'(a1,a60, I6)') '&','                            lstart = ', lstart(l)
        WRITE(*,'(a1,a60, I6)') '&','                              lend = ', lend(l) 
     ENDDO
     
     
  ELSE IF(filetype.EQ."ascii") THEN   ! read from 'tdat/pesl.dat' ascii file 
     
     !

  
     WRITE(*,'(a1,a60)')'&',' pes::  read from tdat/pes.dat.'
     npes   = 10  
     OPEN(npes,FILE='tdat/pesl.dat')              !open spectra file
     READ(npes,'(a1,1x,2i6)') ctmp, nl, ntot    !'&', nof partial waves, nof states
     
     lmax = nl - 1  
     
     ALLOCATE( w2e(0:lmax) )
     
     ne_tot_l = INT(ntot/nl)                  ! fails if nl = 1,3,5,...
     
     WRITE(*,*) '& pes           nl = ', nl
     WRITE(*,*) '& pes         lmax = ', lmax
     WRITE(*,*) '& pes     ne_total = ', ntot
     WRITE(*,*) '& pes     ne_tot_l = ', ne_tot_l
     
     
     partial_wave_l:DO  l = 0, lmax                  ! read for each partial wave
     

        w2e(l)%net = ne_tot_l
     
        ALLOCATE( w2e(l)%ent(1:w2e(l)%net ) )
        ALLOCATE(  w2e(l)%ct(1:w2e(l)%net ) )
        ALLOCATE(  w2e(l)%pt(1:w2e(l)%net ) )
        
        
        READ(npes,'(a1,1x,2i6)') ctmp, l_pes, ntmp    
        
        WRITE(*,*) '& pes:     l = ', l_pes
        
        eigenstate_of_partial_wave_l: DO  ie = 1,  w2e(l)%net 
           
           READ(npes,'(4e14.6)') ee, y2, yr, yi
           
           w2e(l)%ent(ie) = ee
           w2e(l)%ct(ie)  = yr + zi* yi
           w2e(l)%pt(ie)  = y2
           
        ENDDO eigenstate_of_partial_wave_l
     ENDDO partial_wave_l
     
     CLOSE(npes)
     !
     !end of h2e TD coefficients
     
     
     
     ! start off TD coeff checks
     
     IF( ntot.NE.SUM(w2e(:)%net) ) THEN
        WRITE(*,'(a1,a60)') '&','inconsistent data in coefficients file: ntot.ne.sum(ne_l)'
        !     PRINT*, '&', '   nof states        ntot = ',       SIZE(cct)/2 
        PRINT*, '&', '                 sum ne_l = ', SUM(w2e(:)%net) 
        PRINT*, '&    exit the calculation.'
        STOP
     ENDIF
     

          
     DO l = 0 , lmax
        WRITE(*,'(a1,a60)') '&','................................'
        WRITE(*,'(a1,a60, I6)') '&','    angular symmetry block     l2e = ', l
        WRITE(*,'(a1,a60, I6)') '&',' nof 2e states at l2e         ne_l = ', w2e(l)%net 
     ENDDO

     
  ELSE
     
     PRINT*, "#  Unknown file type = ", filetype
     PRINT*, "#  Exit "
     STOP
  ENDIF


     
         
     !check
     
  DO  l = 0, lmax
     PRINT*,"---------------- l = ", l
     DO ie = 1, 10
        !        PRINT*,ie, w2e(l)%ent(ie), cct(ie+lstart(l)-1), cct(ie + lstart(l) - 1 + ntot)
        PRINT*,ie, w2e(l)%ent(ie), REAL(w2e(l)%ct(ie)), AIMAG(w2e(l)%ct(ie))
     ENDDO
  ENDDO
  
     
     
     WRITE(*,'(a1,a60)'      ) '&','       pulse (a.u.):'
     WRITE(*,'(a1,a60,E10.3)') '&','       peak intensity       i0 = ', i0
     WRITE(*,'(a1,a60,G15.4)') '&','            frequency    omega = ', omega
     WRITE(*,'(a1,a60,G15.4)') '&','                pulse duration = ', tau
     WRITE(*,'(a1,a60)'      ) '&','    atomic system (a.u.):'
     WRITE(*,'(a1,a60,I5)   ') '&','  saved 2e states per symmetry = ', nmax
     WRITE(*,'(a1,a60,I5)   ') '&',' 2e angular symmetries l2e_max = ', lmax
     WRITE(*,'(a1,a60,E10.3)') '&','                       e2e_cut = ', en_cut
     WRITE(*,'(a1,a60,3I8)  ') '&','  used  2e states         ntot = ', ntot, SUM(w2e(:)%net)
     WRITE(*,'(a1)')'&'
     WRITE(*,'(a1,a60,G15.3)')'&',' input read finished ' 
     WRITE(*,'(a1)')'&'
     WRITE(*,'(a1)')'&'
  
  
  
     ! end of TD coeff checks

     !DO l = 0, lmax
     !CALL init_w2e_ct( w2e(l),l)
     !CALL init_w2e_v12(w2e(l),l)
     !ENDDO
     
  
     !
     ! start CI coefficients file
     !     ! read 2e-wavefunction  | phi_nL > 
     !

  
     loop_L_symmetries: DO li = 0, lmax
     
        !
        ! 1-e configurations included     ! ev(n),      n = 1, 2,  ..., N_L
        !                                 ! cv(ie,ic)   i = 1, 2, ...., N_L

        CALL hfile(nwf2e,"dat","w2eb","bin",li)                       !L symmetry
        
        READ(unit=nwf2e)    l, ls
  
        IF(li.NE.l) THEN
           WRITE(*,*) '& incosistent files'
           WRITE(*,*) '&               li = ', li
           WRITE(*,*) '&                l = ', l
           STOP
        ENDIF
     
        !
        !  read configurations data ( cfg-L.inp)
        !
        
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




        w2e(l)%l   = l
        w2e(l)%s   = ls
        w2e(l)%nch = ncs   !channels
        w2e(l)%ncf = nhx   !cfg ..
        w2e(l)%nev = n2e
     

     
        ALLOCATE( w2e(l)%ev(1:n2e)  )          ! ev = ev1 + ev2 (non-interacting e-)
        ALLOCATE( w2e(l)%cv(1:n2e,1:nhx))     !'n' energy eigenstate index, 'c' --> configuration index        
        ALLOCATE( w2e(l)%n1(1:nhx)   )
        ALLOCATE( w2e(l)%l1(1:nhx)   )
        ALLOCATE( w2e(l)%n2(1:nhx)   )
        ALLOCATE( w2e(l)%l2(1:nhx)   )
        !
        !ALLOCATE( w2e(l)%pvt(  1:nhx))           ! all E
        !ALLOCATE( w2e(l)%pvt_b(1:nhx))         !  E   < E+
        !ALLOCATE( w2e(l)%pvt_s(1:nhx))         !  E+  < E < E++
        !ALLOCATE( w2e(l)%pvt_d(1:nhx))         !  E++ < E

        
        !
        ! read energies/coefficients from CI calculation
        !
     
    
        read_energies_ci_2e:DO ie = 1, n2e
           READ(nwf2e) w2e(l)%ev(ie)                      
           READ(nwf2e) w2e(l)%cv(ie, :)    ! v(ie,jc)
        ENDDO read_energies_ci_2e
        




        !
        ! re_index cfg(s) according to ie,ic
        !


     
        ic = 0
        channels_index_to_cfg_index:DO k1 = 1, w2e(l)%nch
           
           it1 = nhf(k1)
           it2 = lhf(k1)
           it3 =  ll(k1)
           it4 = nllmin(k1)
           
           DO k2 = 1, nd(k1)
              
              ic = ic + 1
              
              w2e(l)%n1(ic) = it1
              w2e(l)%l1(ic) = it2
              w2e(l)%l2(ic) = it3
              w2e(l)%n2(ic) = it4 + k2  - 1
              
           ENDDO
        ENDDO channels_index_to_cfg_index
        !
        !
        
        
        WRITE(*,'(a1,a60)') '&','----------------------------------'
        WRITE(*,'(a1,a60,I6)') '&','        ang.  symmetry      l2e = ', w2e(l)%l
        WRITE(*,'(a1,a60,I6)') '&','        spin  symmetry      s2e = ', w2e(l)%s
        WRITE(*,'(a1,a60,I6)') '&',' nof zero-2e-states         n2e = ', w2e(l)%nev
        WRITE(*,'(a1,a60,I6)') '&',' nof config. series         ncs = ', w2e(l)%nch
        WRITE(*,'(a1,a60,I6)') '&',' matrix dimension           nhx = ', nhx
        WRITE(*,'(a1,a60,I6)') '&',' matrix dimension            ic = ', ic
        WRITE(*,'(a1,a60)'   ) '&',' read configuration finished.'
        
        
        PRINT*, w2e(l)%n1(1),    w2e(l)%l1(1),    w2e(l)%l2(1),    w2e(l)%n2(1)
        PRINT*, w2e(l)%n1(1831), w2e(l)%l1(1831), w2e(l)%l2(1831), w2e(l)%n2(1831)
        PRINT*, w2e(l)%n1(3661), w2e(l)%l1(3661), w2e(l)%l2(3661), w2e(l)%n2(3661)
        PRINT*, w2e(l)%n1(5491), w2e(l)%l1(5491), w2e(l)%l2(5491), w2e(l)%n2(5491)
        PRINT*, w2e(l)%n1(7320), w2e(l)%l1(7320), w2e(l)%l2(7320), w2e(l)%n2(7320)
        
        !
        !
        
        DEALLOCATE(nhf,lhf,ll,nllmin,nllmax,nd)
        
        
     END DO loop_L_symmetries
     
     !
     !
     !
     !
     !
   !!! FINISH READING  C(T) AND V12 coefficients
     !
     !
     !
     
     PRINT*, "tst:: lmax = ", lmax
 
     !         
     !   P(a,L) =  | Sum_{n} C_(nL)(t) * v^(nL)_(a) |^2  == |<Phi_(a;L)|Psi(t) >|^2
     !         
     !
     ! a  = (n1 l1, n2 l2)  
     ! C_(nL)(t)           : laser interaction coupling coefficient
     ! v^(nL)_(a)          : e-e  interaction  coupling coefficient

  
     ALLOCATE(    ndi_l(0:lmax) )
     ALLOCATE(    nsi_l(0:lmax) )


     ndi_l =  0
     nsi_l =  0
     find_si_di_thresholds:DO l = 0, lmax
        

        n_l2e = w2e(l)%net

        si:DO  ie = 1, n_l2e     ! find SI threshold energy  for symmetry L        
           nsi_l(l) = ie
           !        IF(  w2e(l)%ent(ie) > 0.0_dpk ) EXIT
           IF(  w2e(l)%ent(ie) > en_ion_1 ) EXIT        
        ENDDO si
        
     
     di:DO  ie = 1, n_l2e     ! find DI threshold energy  for symmetry L       
        ndi_l(l) = ie
        !        IF( ( w2e(l)%ent(ie) + en_ion_1) > 0.0_dpk ) EXIT
        IF( ( w2e(l)%ent(ie) ) > 0.0_dpk ) EXIT        
     ENDDO di
   
     WRITE(*,*) "tdse_sdi_pop::   E(", nsi_l(l),  ") = ", w2e(l)%ent(nsi_l(l))
     WRITE(*,*) "tdse_sdi_pop::   E(", ndi_l(l),  ") = ", w2e(l)%ent(ndi_l(l))

  ENDDO find_si_di_thresholds

  !
  !  
  !
  
  ALLOCATE( pop_bi_0(0:lmax) )   ! E   < E+
  ALLOCATE( pop_si_0(0:lmax) )   ! E+  < E < E++
  ALLOCATE( pop_di_0(0:lmax) )   ! E++ < E  
  ALLOCATE(pop_l2e_0(0:lmax) )
  !
  pop_bi_0 = 0.0_dpk
  pop_si_0 = 0.0_dpk
  pop_di_0 = 0.0_dpk
  pop_l2e_0 = 0.0_dpk
  !
  partial_probability_Lv:DO  l = 0, lmax

     PRINT*, '#          l2e = ',l
     PRINT*, '#      ncf = ',w2e(l)%ncf

     ALLOCATE( w2e(l)%pvt(  1:w2e(l)%ncf ))           ! all E
     ALLOCATE( w2e(l)%pvt_b(1:w2e(l)%ncf ))         !  E   < E+
     ALLOCATE( w2e(l)%pvt_s(1:w2e(l)%ncf ))         !  E+  < E < E++
     ALLOCATE( w2e(l)%pvt_d(1:w2e(l)%ncf ))         !  E++ < E
        

     
     w2e(l)%pvt   = 0.0_dpk
     w2e(l)%pvt_b = 0.0_dpk
     w2e(l)%pvt_s = 0.0_dpk
     w2e(l)%pvt_d = 0.0_dpk
     
     
     population_of_cfg_a:  DO ic = 1, w2e(l)%ncf  ! cfg index  
        
        pop_b = 0.0_dpk        
        bb:DO ie = 1, nsi_l(l) - 1              ! E < E+
           
           pop_ic = ABS( w2e(l)%cv(ie,ic) * w2e(l)%ct(ie) )**2 
           pop_b = pop_b + pop_ic           
        ENDDO bb
        !
        pop_s = 0.0_dpk 
        ss:DO ie = nsi_l(l), ndi_l(l)-1        ! E+ < E < E++
           
           pop_ic = ABS( w2e(l)%cv(ie,ic) * w2e(l)%ct(ie) )**2
           pop_s = pop_s + pop_ic
        ENDDO ss        
        !        
        pop_d = 0.0_dpk
        dd:DO ie = ndi_l(l), w2e(l)%net         !  E++ < E
           
           pop_ic = ABS( w2e(l)%cv(ie,ic) * w2e(l)%ct(ie) )**2
           pop_d = pop_d + pop_ic           
        ENDDO dd

        w2e(l)%pvt_b(ic) = pop_b
        w2e(l)%pvt_s(ic) = pop_s
        w2e(l)%pvt_d(ic) = pop_d
          w2e(l)%pvt(ic) = pop_b + pop_s + pop_d 
        
     ENDDO population_of_cfg_a

      
     pop_bi_0(l) =  SUM( w2e(l)%pvt_b )    ! sum over cfgs
     pop_si_0(l) =  SUM( w2e(l)%pvt_s )
     pop_di_0(l) =  SUM( w2e(l)%pvt_d )

     
  ENDDO partial_probability_Lv

  !
  
  pop_l2e_0 = pop_bi_0 + pop_si_0 + pop_di_0   ! sum over b,s,d
  pop_all_0 = SUM(pop_l2e_0)                   ! sum over l2e
  
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
  ALLOCATE(pop_d_sp(0:lmax))
  ALLOCATE(pop_d_pd(0:lmax))
  ALLOCATE(pop_d_df(0:lmax))     
  ALLOCATE(pop_d_sd(0:lmax))
  ALLOCATE(pop_d_pf(0:lmax))
  !
  pop_bb = 0.0_dpk
  pop_ss = 0.0_dpk
  pop_dd = 0.0_dpk
  !
  !  
  pop_d_ss = 0.0_dpk
  pop_d_pp = 0.0_dpk
  pop_d_dd = 0.0_dpk
  pop_d_sp = 0.0_dpk
  pop_d_pd = 0.0_dpk
  pop_d_df = 0.0_dpk
  pop_d_sd = 0.0_dpk
  pop_d_pf = 0.0_dpk
  
  
  angular_symmetry_2e_zero_order:DO l = 0, lmax

     PRINT*, "l2eee = ", l


     !di E > E++
     pop_b = 0.0_dpk
     pop_s = 0.0_dpk
     pop_d = 0.0_dpk     
     sum_b = 0.0_dpk
     sum_s = 0.0_dpk
     sum_d = 0.0_dpk
     !
     sum_d_ss = 0.0_dpk
     sum_d_pp = 0.0_dpk
     sum_d_dd = 0.0_dpk
     sum_d_sp = 0.0_dpk
     sum_d_pd = 0.0_dpk
     sum_d_df = 0.0_dpk
     sum_d_sd = 0.0_dpk
     sum_d_pf = 0.0_dpk

     
     sum_all_configurations: DO ic = 1, w2e(l)%ncf

  
        e1  = e1e( w2e(l)%l1(ic), w2e(l)%n1(ic) )
        e2  = e1e( w2e(l)%l2(ic), w2e(l)%n2(ic) )
        le1 = w2e(l)%l1(ic)
        le2 = w2e(l)%l2(ic)
        
        pe   = w2e(l)%pvt(ic)   !  (e1,e2 > 0)

        
        pick_channels:IF( ( e1>0).AND.(e2>0) )  THEN
           
           pe_d = w2e(l)%pvt_d(ic)  !  
           
           pop_d = pop_d  + pe       !  (e1,e2 > 0)           
           sum_d = sum_d +  pe_d     !  (e1,e2 > 0) and (E > E++)

           
           IF((le1.EQ.0).AND.(le2.EQ.0)) sum_d_ss = sum_d_ss + pe_d
           IF((le1.EQ.1).AND.(le2.EQ.1)) sum_d_pp = sum_d_pp + pe_d
           IF((le1.EQ.2).AND.(le2.EQ.2)) sum_d_dd = sum_d_dd + pe_d
           IF((le1.EQ.0).AND.(le2.EQ.1)) sum_d_sp = sum_d_sp + pe_d
           IF((le1.EQ.1).AND.(le2.EQ.0)) sum_d_sp = sum_d_sp + pe_d
           IF((le1.EQ.1).AND.(le2.EQ.2)) sum_d_pd = sum_d_pd + pe_d
           IF((le1.EQ.2).AND.(le2.EQ.1)) sum_d_pd = sum_d_pd + pe_d
           IF((le1.EQ.2).AND.(le2.EQ.3)) sum_d_df = sum_d_df + pe_d
           IF((le1.EQ.3).AND.(le2.EQ.2)) sum_d_df = sum_d_df + pe_d
           IF((le1.EQ.0).AND.(le2.EQ.2)) sum_d_sd = sum_d_sd + pe_d
           IF((le1.EQ.2).AND.(le2.EQ.0)) sum_d_sd = sum_d_sd + pe_d
           IF((le1.EQ.1).AND.(le2.EQ.3)) sum_d_pf = sum_d_pf + pe_d
           IF((le1.EQ.3).AND.(le2.EQ.1)) sum_d_pf = sum_d_pf + pe_d
           
        ELSE IF ( (e1<0).AND.(e2<0) ) THEN
           
           pe_b  = w2e(l)%pvt_b(ic)
           
           pop_b = pop_b + pe        ! (e1,e2 < 0)   
           sum_b = sum_b + pe_b      ! (e1,e2 < 0)  and   (E < E+)
           
        ELSE

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
     pop_d_sp(l) = sum_d_sp
     pop_d_pd(l) = sum_d_pd
     pop_d_df(l) = sum_d_df
     pop_d_sd(l) = sum_d_sd
     pop_d_pd(l) = sum_d_pf
     
     
     
  ENDDO angular_symmetry_2e_zero_order


  pop_l2e_cv  = pop_bb + pop_ss + pop_dd   ! sum over b,s,d
  pop_all_cv  = SUM(pop_l2e_cv)            ! sum over l2e

  
  !
  !
  ! probabilities in terms of the CI 2e-functions  psi(t) = S_NL C_NL(t) Phi_NL(r_1,r_2)
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
     
     pop_bi(l) = SUM( ABS(  w2e(l)%ct(        1: nsi_l(l)-1) )**2 )  ! E_0 < E < E+
     pop_si(l) = SUM( ABS(  w2e(l)%ct( nsi_l(l): ndi_l(l)-1) )**2 )  ! E+  < E < E++
     pop_di(l) = SUM( ABS(  w2e(l)%ct( ndi_l(l): w2e(l)%net) )**2 )  ! E++ < E
         
     pop_l2e(l) = SUM( ABS( w2e(l)%ct(1:w2e(l)%net) )**2 )                ! all over E

     
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
  PRINT*, '  E,e1,2                           pop_cv(l) = ', pop_l2e_cv, pop_all_cv
  PRINT*, '  E                                   pop(l) = ', pop_l2e,   pop_all
  PRINT*, '  e1,e                              pop_0(l) = ', pop_l2e_0, pop_all_0
  PRINT*, 'x'
  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx '


  PRINT*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxx   double ionization'
  PRINT*, 'x'
  PRINT*, ' (e1,e2 > 0) and           E++ < E    pop_dd = ', pop_dd,   SUM(pop_dd)
  PRINT*, 'x          L = 0,1,...,Lmax,       sum over L'
  PRINT*, '   pop_d_ss  = ', pop_d_ss,   SUM(pop_d_ss)
  PRINT*, '   pop_d_pp  = ', pop_d_pp,   SUM(pop_d_pp)
  PRINT*, '   pop_d_dd  = ', pop_d_dd,   SUM(pop_d_dd)
  PRINT*, '   pop_d_sp  = ', pop_d_sp,   SUM(pop_d_sp)
  PRINT*, '   pop_d_pd  = ', pop_d_pd,   SUM(pop_d_pd)
  PRINT*, '   pop_d_df  = ', pop_d_df,   SUM(pop_d_df)
  PRINT*, '   pop_d_sd  = ', pop_d_sd,   SUM(pop_d_sd)
  PRINT*, '   pop_d_pf  = ', pop_d_pf,   SUM(pop_d_pf)
  PRINT*, 'x'  
  
 

    


END PROGRAM tdse_pop_sdi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!EOF
