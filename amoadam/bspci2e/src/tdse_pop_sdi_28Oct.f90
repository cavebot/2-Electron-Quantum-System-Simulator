!
! calculates spatial density matrix, rho(r1,r2,t)
!

PROGRAM tdse_pop_sdi
  !
  USE PRECISION
  USE units, only: zi
  USE parameter_tdse_fxd
  USE bs_frb_2e, ONLY: read_target_energies
  USE io
  use netcdf
  USE ncFile
  !
  USE wf_1e
  USE utils,        ONLY: asciifile,datafile
  !
  IMPLICIT NONE
  !
  !
  INTEGER, PARAMETER                     :: l1e_max = 5
  INTEGER                                :: i, j, p, q, k
  INTEGER                                :: li, l
  INTEGER                                :: n_plot, l_plot
  INTEGER                                :: ie, ic, k1, k2
  ! coefficients
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   ::  ct        ! time-dependent coeff
  COMPLEX(dpk)                           ::  ct_nl 
  !field! 
  REAL(dpk)                              :: i0,omega,tau
  !
  REAL(dpk)                              :: time
  !grid!
  INTEGER,  ALLOCATABLE,DIMENSION(:)          :: ipw
  INTEGER,  ALLOCATABLE,DIMENSION(:)          :: pw

  INTEGER                                     :: it, ir1, ir2, n_lab !
  INTEGER                                     :: la, lb, ie_min, ie_max, ill !
  REAL(dpk),      ALLOCATABLE, DIMENSION(:,:) :: rho_t              ! rho_t(ir1,ir2)
  INTEGER                                     :: il, i1, i2, i1_start
  INTEGER                                     :: n_1, l_1, n_2, l_2
  REAL(dpk),    ALLOCATABLE,   DIMENSION(:,:) :: f2_ln        !(lalb, r2)
  REAL(dpk),    ALLOCATABLE, DIMENSION(:,:,:) :: f12_l        !(lalb, r1,r2)
  COMPLEX(dpk), ALLOCATABLE, DIMENSION(:,:,:) :: f12_ln       !(lalb, r1,r2)
  REAL(dpk), ALLOCATABLE,      DIMENSION(:,:) :: rho_12       !(r1,r2)

  !field-free atomic configurations
  INTEGER                                     :: ncs, nhx
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: ne_l, lstart, lend
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: en_l
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: ev                    !nhx
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: cv_nc                 !nhx
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: nhf, lhf, ll          !ncs
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: nll_min, nll_max, nd  !ncs nd_2 = nl2_max - nl2_min + 1 (same n1,l1,l2) 
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: n1, l1                !
  INTEGER,   ALLOCATABLE, DIMENSION(:)        :: n2, l2                !
  INTEGER                                     :: l2e, n2e, ls
  INTEGER                                     :: nwf2e
  !
  REAL(DPK), DIMENSION(:,:), POINTER          :: pl                  ! p_l(r), p'_l(r)
  REAL(DPK), DIMENSION(:),   POINTER          :: ri                  ! r_i, dr_i
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:,:)    :: w1e                 ! p_nl(r_i)
  INTEGER                                     :: ns,np               ! l1e_max
  INTEGER,   DIMENSION(:,:), POINTER          :: e1e                 !nl,ns
  !
  CHARACTER(len=6)                            :: form
  REAL(dpk), ALLOCATABLE, DIMENSION(:)        :: tm 
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)      :: ct_re, ct_im
  !io
  CHARACTER(LEN=25)                           :: atomic_file         !i
  CHARACTER(LEN=25)                           ::    coe_file         !i
  CHARACTER(LEN=25)                           ::    ion_file         !o
  CHARACTER(LEN=25)                           ::    inp_file         !i
  CHARACTER(LEN=25)                           ::    out_file         !o
  !
  INTEGER                                     :: nion, ncoe
  !ncFile
  INTEGER                                     :: ncid 
  LOGICAL                                     :: coe_file_exists
  INTEGER                                     :: dim_t 
  INTEGER                                     :: it_index, nt_out, i_out
  INTEGER                                     :: nwf1e
  INTEGER                                     :: nof_1e_states, nof_r_points
  CHARACTER(LEN=6)                            ::   sn, sl


  INTERFACE
     SUBROUTINE get_partial_wave_info(ncs, l1, l2, ipw, pw)
       USE DynamicIntegerArray
       IMPLICIT NONE
       INTEGER,                             INTENT(in)  :: ncs
       INTEGER,             DIMENSION(ncs), INTENT(in)  :: l1
       INTEGER,             DIMENSION(ncs), INTENT(in)  :: l2
       INTEGER,ALLOCATABLE,   DIMENSION(:), INTENT(out) :: ipw
       INTEGER,ALLOCATABLE,   DIMENSION(:), INTENT(out) ::  pw
     END SUBROUTINE get_partial_wave_info
  END INTERFACE
  
  !
  !exe!
  
      !
 !  CALL getarg(1, argv)             ! nof partial waves
 !  READ(argv,*)   lmax              !

  ! define i/o id's and file names


  nwf2e = 9
  nwf1e = 10
  nion  = 11
  ncoe  = 12


  !in
!  atomic_file = "dmx2ebb"      
  !out

  inp_file = 'tinp/tdse_pop_sdi.inp'
  ion_file = 'tdat/pop_di.dat'
  coe_file = 'tdat/coe.dat'

  
  OPEN( ninp, file=inp_file, status = 'UNKNOWN' )
  READ( ninp, * )               n_plot, l_plot
  READ( ninp, * )               ie_min, ie_max
  CLOSE(ninp)



  WRITE(*,'(a1,a60,2I5)') '&','test 1e orbital   (n_plot,  l_plot) = ', l_plot, n_plot
  WRITE(*,'(a1,a60,2I5)') '&', 'energy range:    (ie_min,  ie_max) = ', ie_min, ie_max

  
  WRITE(sl,'(i6)') l_plot
  WRITE(sn,'(i6)') n_plot

  !
  ! now read
  !

  CALL input_tdse_fxd              ! tinp/tdse_bs_fxd_2e.inp
  CALL output_tdse_fxd             ! tout/tdse_bs_fxd_2e.out




  CALL read_target_energies(e1e)             ! read 1e energies (in Ryd)
  e1e = e1e * 0.5_dpk                        ! e1e in a.u. now
  
  
  !
  !
  !  read 1-e partial waves    P_nl(r), l = 0, 1,...., l1e_max
  !
  !
!!%  !l = 0
!!%  !CALL datafile(nwf1e, l,"w1e-")
!!%  !READ(nwf1e) nof_1e_states, nof_r_points
!!%  !CLOSE(nwf1e)
!!%  ! WRITE(*,'(a1,a60,I5)') '&','l1e_max = ', l1e_max
!!%  !WRITE(*,'(a1,a60,I5)') '&','n1e_max = ', nof_1e_states
!!%  !WRITE(*,'(a1,a60,I5)') '&',' nr_max = ', nof_r_points
!!%  !ALLOCATE( w1e(0:l1e_max, nof_1e_states, nof_r_points) )
  !
  !     w1e = w1e(l,n,r_i)
  !
!!%
!!%  w1e = 0.0_dpk 
!!%  read_1e_wfs:DO l = 0, l1e_max
!!%
!!%  WRITE(*,'(a1,a60,G15.2)') '&',' l = ', l
!!%
!!%     CALL read_all_wf1e(l, pl, ri)
!!%
!!%
!!%     IF(l == l_plot ) THEN
!!%        out_file = "out/pl-l"//TRIM(ADJUSTL(sl))//"_n"//TRIM(ADJUSTL(sn))//".out"
!!%        OPEN(nout, file=out_file, status = 'UNKNOWN' )
!!%        DO i = 1, nof_r_points
!!%           WRITE(nout,*) ri(i), pl(n_plot,i)
!!%        ENDDO
!!%        CLOSE(nout)
!!%     ENDIF
!!%
!!%     store_in_matrix:IF(nof_1e_states.EQ.SIZE(pl,dim=1).AND.nof_r_points.EQ.SIZE(pl,dim=2)) THEN
!!%        w1e(l,:,:) = pl
!!%     ENDIF store_in_matrix
!!%
!!%  ENDDO read_1e_wfs
!!%
!!%  DEALLOCATE(pl)
!!%
!!%  WRITE(*,'(a1,a60,G15.2)') '&',' r_min = ', ri(1)
!!%  WRITE(*,'(a1,a60,G15.2)') '&',' r_max = ', ri(nof_r_points)


  !
  ! check size of w1e
  !
!!%
!!%  WRITE(*,'(a1,a60,3I5)') '&','w1e sizes =  ', SIZE(w1e,  dim=1), SIZE(w1e,  dim=2), SIZE(w1e,  dim=3)
!!%  WRITE(*,'(a1,a60,I5)')   '&','ri   size =  ', SIZE(ri)

  !
  ! plot P(r)_(nplot,l_plot)     r_min < r < r_max
  !
!!%  
!!%  out_file = "out/w1e-l"//TRIM(ADJUSTL(sl))//"_n"//TRIM(ADJUSTL(sn))//".out"
!!%  OPEN(nout, file=out_file, status = 'UNKNOWN' )
!!%  DO i = 1, nof_r_points
!!%     WRITE(nout,*) ri(i), w1e(l_plot,n_plot,i)
!!%  ENDDO
!!%  close(nout)
!!%  

  WRITE(*,'(a1,a60)') '&',' 1e partial waves read.'

  !
  !  read time-dependent coefficients        C_nl(t)
  !


  WRITE(*,'(a1,a60,a20)') '&',' open coefficient file:', coe_file

  INQUIRE(file=coe_file,exist=coe_file_exists)  
  IF(coe_file_exists) THEN 
     WRITE(*,'(a1,a60)') '&',' read from coefficient file'
  ELSE
     WRITE(*,'(a1,a60)') '&',' coefficient file is not present.aborting'
     STOP
  ENDIF
  


  !
  !
  !  read 2e time-dependent coefficients data
  !
  !
  ! coef(t) file
  !
  ! intensity, frequency, total pulse duration, maximum_L_2e, maximum_energy_2e
  ! en(N,L)
  ! ntot
  ! ct(N,L)
  !
  

  OPEN(ncoe,file=coe_file,form='unformatted',access='sequential') 
  READ(ncoe)  i0, omega, tau, lmax, en_cut

  WRITE(*,'(a1,a60)'      ) '&','       pulse (a.u.):'
  WRITE(*,'(a1,a60,E10.3)') '&','       peak intensity       i0 = ', i0
  WRITE(*,'(a1,a60,G15.4)') '&','            frequency    omega = ', omega
  WRITE(*,'(a1,a60,G15.4)') '&','                pulse duration = ', tau
  WRITE(*,'(a1,a60)'      ) '&','    atomic system (a.u.):'
  WRITE(*,'(a1,a60,I5)   ') '&','    maximum 2e angular l2e_max = ', lmax
  WRITE(*,'(a1,a60,E10.3)') '&','    maximum 2e energy  e2e_max = ', en_cut


  !
  ! nmax to be removed
  !
  
  ALLOCATE(   en_l( 1:nmax, 1:lmax))
  ALLOCATE(   ne_l(         0:lmax))

  read_energies_l:DO  l = 0, lmax
     READ(ncoe)  ne_l(l), ( en_l(ie,l), ie = 1, ne_l(l))
  ENDDO read_energies_l 
  READ(ncoe) ntot
  ALLOCATE( ct(2*ntot) )
  READ(ncoe) time, (ct(ie), ie = 1, 2*ntot)
  CLOSE(ncoe)

  
  ! 
  ! calculate starting elements of the 'L' symmetry blocks
  !

  ALLOCATE( lstart( 0:lmax) )
  ALLOCATE( lstart( 0:lmax) )

  lstart(0) = 1
  DO l = 1, lmax
     lstart(l) = lstart(l-1) + ne_l(l)
  END DO

  lend = lstart + ne_l - 1

  
  WRITE(*,'(a1,a60,I6)') '&',' total nof 2e states    ntot = ', ntot
  DO l = 0 , lmax
     WRITE(*,'(a1,a60, I6)') '&','    angular symmetry block     l2e = ', l 
     WRITE(*,'(a1,a60, I6)') '&','                            lstart = ', ne_l_start(l)
     WRITE(*,'(a1,a60,2I6)') '&','                              lend = ', ne_l_start(l) - ne_l(l), ne_l(l+1), lend(l)
     WRITE(*,'(a1,a60, I6)') '&',' nof 2e states at l2e         ne_l = ',       ne_l(l)
  ENDDO
  WRITE(*,'(a1)')'&'
  WRITE(*,'(a1,a60,G15.3)')'&',' input read finished ' 
  WRITE(*,'(a1)')'&'
  WRITE(*,'(a1)')'&'


  OPEN(nion,file=ion_file)                            ! ionization file
  
  loop_L_symmetries: DO li = 0, lmax

     ! 
     ! read 2e-wavefunction  | phi_nL > 
     !
     ! 1-e configurations included     ! ev(n),     n = 1, 2,  ..., N_L
     !                                 ! C(ie,ic)   i = 1, 2, ...., N_L
     !
     !
     
     WRITE(*,'(a1,a60,G15.3,a2)')'&',' open file ' 
     CALL hfile(nwf2e,"dat","w2eb","bin",li)                       !L symmetry
     WRITE(*,'(a1,a60,G15.3,a2)')'&',' file opened ', time, 'a.u.' 
     
     READ(unit=nwf2e)    l, ls
     IF((li-1).NE.l) THEN
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
     
     READ(nwf2e)  (    nhf(k1),    k1 = 1, ncs)
     READ(nwf2e)  (    lhf(k1),    k1 = 1, ncs)
     READ(nwf2e)  (     ll(k1),    k1 = 1, ncs)
     READ(nwf2e)  ( nllmin(k1),    k1 = 1, ncs)
     READ(nwf2e)  ( nllmax(k1),    k1 = 1, ncs)
     READ(nwf2e)  (     nd(k1),    k1 = 1, ncs)
     READ(nwf2e)          nhx                ! nof configurations (n1,l1;n2,l2)
     READ(nwf2e)          n2e


     !
     
     ALLOCATE(   ev(n2e)    )
     ALLOCATE(   cv_nc(n2e,nhx))       !'n' energy eigenstate index, 'c' --> configuration index


     !
     ! read energies/coefficients from CI calculation
     !
     read_energies_ci_2e:DO      ie = 1, SIZE(cv,dim=1)
        READ(nwf2e)  ev(ie)
        READ(nwf2e) (cv_nc(ie, ic), ic = 1, SIZE(cv,dim=2))
     ENDDO read_energies_ci_2e
     WRITE(*,*) "tdse_pop_sdi::              read ev2e and cv2e done"

     ! re_index cfg(s) according to ie,ic

     
     ALLOCATE(    n1(nhx) )
     ALLOCATE(    n2(nhx) ) 
     ALLOCATE(    l1(nhx) )
     ALLOCATE(    l2(nhx) )

     
     ic = 0
     channels_index_to_cfg_index:DO k1 = 1, ncs
        DO k2 = 1, nd(k1)

           ic = ic + 1
           
           n1(ic) = nhf(k1)
           l1(ic) = lhf(k1)
           l2(ic) =  ll(k1)        
           n2(ic) = nllmin(k1) + k2  

        ENDDO
     ENDDO single_index_channels
     
     
     WRITE(*,'(a1,a60,I6)') '&','        ang.  symmetry      l2e = ', l
     WRITE(*,'(a1,a60,I6)') '&','        spin  symmetry      s2e = ', ls
     WRITE(*,'(a1,a60,I6)') '&',' nof zero-2e-states         n2e = ', n2e
     WRITE(*,'(a1,a60,I6)') '&',' nof config. series         ncs = ', ncs
     WRITE(*,'(a1,a60,I6)') '&',' matrix dimension           nhx = ', nhx
     WRITE(*,'(a1,a60,I6)') '&',' matrix dimension            ic = ', ic
     WRITE(*,'(a1,a60)'   ) '&',' read configuration finished.'


 
     !         
     !   P(a,L) = Sum_n( | S_{n} C_(nL)(t) * C^(nL)_(a) |^2 )
     !         
     !
     ! a  = (n1 l1, n2 l2)  
     ! C_(nL)(t)           : laser interaction coupling coefficient
     ! C^(nL)_(a)          : e-e  interaction  coupling coefficient
     !
!
! time-dependent coefficient  ct
!        

     population_of_cfg_a:DO ic = 1, nhx
        
        ct_nl = 0.0_dpk
        sum_over_all_2e_energies:DO ie = 1, n2e
           
           ct_nl = ct_nl + (ct(lstart(li) + ie ) +  zi* ct(lstart(li) + ie + ntot) ) * cv_nc(ie,ic)   !*note
           
        ENDDO sum_over_all_2e_energies

     
        pop_a(ic,li) = REAL(ct_nl)**2 + AIMAG(ct_nl)**2
  ENDDO population_of_cfg_a
     

     ENDDO loop_L_symmetries



     DEALLOCATE(ev,cv) !ok
     DEALLOCATE(ipw,pw)!
     DEALLOCATE(f2_ln,f12_ln, f12_l)
     DEALLOCATE(nhf,lhf,ll,nllmin,nllmax,nd)

  ENDDO states_2e_Ln
  

  
  !  rho_12 = 0.0D+00
  !  DO li = 1,  lmax 
  !  ENDDO
  
  !     rho_12 = rho_12 + SUM( ABS(f12_l)**2, dim = 1 )     ! sum over l1,l2




  l = 12
  CALL asciifile(nwf1e, l,"rho_")
!  WRITE(nwf1e,*) SIZE(rho_12, dim = 1)
!  save_density_matrix: DO i = 1, SIZE(rho_12, dim = 1)     
  save_density_matrix: DO i = 1, 100      
!     WRITE(nwf1e,'(50E14.5)') ( rho_12(i,j), j = 1, SIZE(rho_12, dim = 2))
     WRITE(nwf1e,'(50E14.5)') ( rho_12(i,j), j = 1, 100)
!     WRITE(nwf1e,'(<nof_r_points>E14.5)') ( rho_12(i,j), j = 1, SIZE(rho_12, dim = 2))
  ENDDO save_density_matrix
  CLOSE(nwf1e)
  

  
  !  DO i = 1, SIZE(rho_12, dim = 1)
  !     DO j = 1, SIZE(rho_12, dim = 2)
  !        WRITE(nwf1e,*) i,j, rho_12(i,j)
  !     ENDDO
  !  ENDDO
  !  CLOSE(nwf1e)


END PROGRAM tdse_rho
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!EOF
