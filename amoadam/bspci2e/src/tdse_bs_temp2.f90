!!!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                                         written by xian tang, jan 1990
!                                          changes made by  hr, feb 1990
!   changes made 3/31/90 by hr,xtang,jzhang,laan:
!   changed to remove the rapid  oscillating part of the matrix
!   pulse shape changed to gaussian
!   pulse duration changed from -0.5* tau to 0.5*tau
!
!   this program uses the nag routine d02baf to integrate the diff. equation
!   and it can be used for restarting an old job with a unit 33, with time and
!   vector information
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! propagates the coefficient vector  C(t)
!
! y(t) expanded on a eigenstate basis:

! y(t) = Sum_(nl) C_(nl) * y_(nl)
!
! H_0 * y_nl = E_nl y_nl,     <y_nl|y_n'l> = d_nn'
!
!
! [ H_0 + V(t) ] y(t) = i dy/dt ==> dC/dt = |H+V(t)| * C(t)
!
! input: tinp/tdse_bs.inp    (some integration parameters)
!        tinp/pulse.inp      (pulse parameters)
!
!         dat/heLL+1.nc      (field-free data (E_nl, V(nl,n'l+1) in netCDF format
!        or
!         dat/dmx2e-v.LL+1.dat (as above in fortran binary)
!
!

PROGRAM tdse_bs_fxd_2e
  !
  USE PRECISION
  USE units
  USE parameter_tdse_fxd
  USE deriv_tdse_fxd
  USE io
  USE utils, ONLY:tafile
  USE pulse
  use netcdf
  USE ncFile
  USE DynamicIntegerArray

  !
  IMPLICIT NONE
  !
  INTEGER,                     PARAMETER :: nfft  = 16384
  !
  !
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   ::  y
  REAL(dpk), ALLOCATABLE, DIMENSION(:)   ::  yp
  !field!
  REAL(dpk)                              :: i0   ! e0 = sqrt(i0): electric field
  REAL(dpk)                              :: tp   !
  INTEGER                                :: nout_T
  !
  REAL(dpk)                              :: t, t1, tim, ti
  REAL(dpk)                              :: pop
  REAL(dpk)                              :: pop_g, pop_b, pop_k, pop_l
  !pes!
  CHARACTER(len=25)                      :: pes_file
  INTEGER                                :: npes
  INTEGER                                :: ntot_l, ne
  !
  INTEGER                                :: i, k, k1, nb, iu, lstart, kt,lk
  INTEGER                                :: nout1, nout2, nt, li, l, ikt, ik, ii
  !hohg!
  REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: dpr, dpi
  REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: time
  REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: rt
  REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: at
  REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: vt
  CHARACTER(len=6)                         :: form

  REAL(dpk), ALLOCATABLE, DIMENSION(:)     :: tm
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)   :: ct_re, ct_im
  INTEGER                                  :: ndip
  CHARACTER(len=25)                        :: dip_file
  INTEGER                                  :: tSteps
  !
  !nag
  REAL(dpk), ALLOCATABLE, DIMENSION(:,:)   :: w       !nag, work array
  INTEGER                                  :: ifail     !nag/d02baf
  !  INTEGER                                :: iwork(5)  !ode
  !  REAL(dpk), ALLOCATABLE, DIMENSION(:)   :: work      !ode nag, work array
  !io
  !  CHARACTER(LEN=25)                      :: atomic_file     !i
  CHARACTER(LEN=25)                      ::    coe_file     !o
  CHARACTER(LEN=25)                      ::    ion_file     !o
  CHARACTER(LEN=25)                      ::    pop_file     !o
  !
  INTEGER                                :: nion, npop, ncoe
  !ncFile
  INTEGER                                :: ncid
  !
  LOGICAL                                :: coe_file_exists
  !
  INTEGER                                :: dim_t
  integer                                :: it_index, nt_out, i_out
  !
  !nag
  EXTERNAL                        d02baf        ! nag
  EXTERNAL                 tdse_v,tdse_v1, tdse_l, tdse_lv        ! user defined
  !exe!

  !
  !  CALL getarg(1, argv)             ! nof partial waves
  !  READ(argv,*)   lmax              !

  ! define i/o id's and file names

  nion  = 10
  npop  = 11
  ncoe  = 12
  ndip  = 14
  npes  = 16

  !in
!  atomic_file = "dmx2ebb"
  !out
  ion_file='tdat/pgt.dat'
  pop_file='tdat/pop.dat'
  coe_file='tdat/coe.dat'
  dip_file='tdat/dip.dat'
  pes_file='tdat/pesl.dat'

  ! now read
  CALL input_tdse_fxd              ! tinp/tdse_bs_fxd_2e.inp
  CALL output_tdse_fxd             ! tout/tdse_bs_fxd_2e.out
  CALL make_space                  ! allocate space for arrays en,dmx


!xxxxxxxxxxx                    hohg        xxxxxxxxxxxxxxxxxxxxxxxxx


  IF(hohg) THEN  !form ='a' is also available.
     form = "v"
     IF(gauge=="l")  form = "l"
     CALL read_system_ncFile(form) ! read en/dmx from netCDF (dat/systemLL+1gauge.nc)
  ENDIF

  ! read en/dmx binaries

  IF(bin_type == "cdf")  THEN
     CALL read_system_ncFile(gauge)        !netCDF (dat/systemLL+1gauge.nc)
  ELSE IF(bin_type =="f90") THEN
     CALL read_system(gauge)               !f90    (dat/systemLL+1gauge.dat)
  ELSE
     WRITE(*,*)"& check for (cdf) or (f90) entries in dat/tdse_bs_fxd_2e.dat."
     STOP
  ENDIF

  CALL read_field                        ! read field parameters (tinp/pulse.inp)

  i0     = e0**2 * i0_td
  tp     = 2.0_dpk * m_pi / omega            ! T_w :         period (a.u.)
  nout_T = 50
  ti     = cycles * tp * start_time

  !  ti        = 0.0D+00

  CALL save_field       (ti, nout_T)              !  out/pulse.out  !
  CALL open_tdse_ncFile (atomic_name, gauge, ncid)
  CALL save_pulse_ncFile(ncid, i0, omega, tau)
  CALL save_en_ncFile   (ncid, n, en)



  !  OPEN(20,FILE='tout/tdse.out')
  WRITE(*,*)'#'
  WRITE(*,*)'#                                  light:'
  WRITE(*,*)'#'
  WRITE(*,*)'#        peak electric field       Eo  = ', e0,        ' a.u.'
  WRITE(*,*)'#             photon energy        Wph = ', omega ,    ' a.u.'
  WRITE(*,*)'#             pulse duration       tau = ', tau,       ' a.u.'
  WRITE(*,*)'#             pulse type     pulseType = ', pulseType
  WRITE(*,*)'#             peak intensity       Io  = ', i0,        ' a.u.'
  WRITE(*,*)'#             field period         tp  = ', tp,        ' a.u.'
  WRITE(*,*)'#             initial time         ti  = ', ti,        ' a.u.'
  WRITE(*,*)'#             nof cycles       cycles  = ', cycles,    ' a.u.'
  WRITE(*,*)'#'
  WRITE(*,*)'#                                   atom:'
  WRITE(*,*)'#'
  WRITE(*,*)'# ionic ground  state   energy    E(+) = ', en_ion_1,  ' a.u.'
  WRITE(*,*)'#       ground  state   energy    E_g  = ', en_ground, ' a.u.'
  WRITE(*,*)'#       maximum kinetic energy  En_cut = ', en_cut,    ' a.u.'
  WRITE(*,*)'#'
  WRITE(*,*)'#                             interaction'
  WRITE(*,*)'#'
  WRITE(*,*)'#      gauge interaction         gauge = ', gauge


  ! restarted or a new calculation?
  WRITE(*,*) '#        coefficient file:', coe_file
  INQUIRE(file=coe_file,exist=coe_file_exists)
  restarted:IF(irestart.EQ.1) THEN
     WRITE(*,*) '# this is a restarted run:'
     IF(.NOT.coe_file_exists) THEN
        WRITE(*,*) '# coefficient file is not present.aborting'
        STOP
     ENDIF
     WRITE(*,*) '# read from coefficient file'

  ELSE
     WRITE(*,*)' This is not a restarted run.'
     IF(coe_file_exists) THEN
        WRITE(*,*) '# coefficient file is not empty. aborting.'
        STOP
     END IF
  ENDIF restarted


  OPEN(ncoe,file=coe_file,form='unformatted',access='sequential') ! coef(t) file


  !   if restart get stored information from the previous run
  IF(irestart.EQ.1) THEN
     READ(ncoe)  i0, omega, tau, lmax, en_cut
  ELSE
     WRITE(ncoe) i0, omega, tau, lmax, en_cut
  ENDIF

  !   Set as Zero energy the 1st ionization threshold of the atom en_ion_1.
  !    en(ne, l) > 0  continuum
  !    en(ne, l) < 0  bound
  !

  ntot = 0
  change_energy_axis: DO  l = 1, lmax+1

     iu = 0
     exclude_high_energies: DO  ne = 1, n(l)

        IF(atomic_name=="he") THEN
           en(ne,l) = 0.5D+00 * en(ne,l) - en_ion_1         ! (in a.u.) !standard 2e-data
        ENDIF

        IF(en(ne,l).GT.en_cut) CYCLE                     ! drop  E > en_cut
        iu = iu + 1
     ENDDO exclude_high_energies

     n(l) = iu
     ntot = ntot + n(l)            ! total nof states (actual size of psi(t))


     check_restart_1:IF(irestart.EQ.1) THEN
        READ(ncoe)  n(l), ( en(ne,l), ne = 1, n(l))
     ELSE                                            !else write in coe file
        WRITE(ncoe) n(l), ( en(ne,l), ne = 1, n(l))
     ENDIF check_restart_1


     WRITE(*,*)'#'
     WRITE(*,*)'#      partial wave     l = ', l - 1
     WRITE(*,*)'#        nof states  n(l) = ', n(l)
     WRITE(*,*)'#'
  ENDDO change_energy_axis


  neq = 2*ntot
  ALLOCATE( y(neq), yp(neq))

  WRITE(*,*)    '# tdse_bs_mch::                              Emax = ', en_cut
  WRITE(*, *)   '# tdse_bs_mch::   nof channels propagated   n_tot = ', ntot


  check_restart_2:IF(irestart.EQ.1) THEN
     READ(ncoe) ntot
  ELSE
     WRITE(ncoe) ntot
  ENDIF check_restart_2

  WRITE(*,*)'# total nof states             ntot = ', ntot
  WRITE(*,*)'# total nof equations           neq = ', neq

  !
  !  Make the diagonal part of matrix  E(nl,nl)
  !

  nb = 0
  assign_diagonal_hamiltonian:   DO   k = 1, lmax+1
     loop_over_states_for_each_l:DO  k1 = 1, n(k)

        nb      = nb + 1
        dag(nb) = en(k1,k)              ! a.u.

     ENDDO loop_over_states_for_each_l
  ENDDO assign_diagonal_hamiltonian

  ! set by hand
  !       dag(1) = en_ground - en_ion_1

  WRITE(*,*) '#  ground state energy E_g = ', dag(1), ' a.u. '

  !
  !   if restart get information from stored values
  !

  check_restart_3:IF(irestart.EQ.1) THEN

     t = -0.5_dpk * tau

215  READ(ncoe,END=210) t1, ( yp(ii), ii = 1, neq)

     t = t1

     !          y = yp

     DO ik = 1, neq
        y(ik) = yp(ik)
     ENDDO


     GOTO 215

210  lstart = INT(0.5_dpk * ( t + 0.5_dpk * tau ) /Tp ) + 1

     t = t + 0.5_dpk * tau

  ELSE

     lstart = 1

  ENDIF check_restart_3


  initial_conditions: IF(lstart.EQ.1) THEN
     y    = 0.0_dpk
     y(1) = 1.0_dpk
     t    = ti
  ENDIF initial_conditions


  nsum(1) = 1
  DO i = 2, lmax+1
     nsum(i) = nsum(i-1) + n(i-1)
     WRITE(*,*) 'i, nsum(i) = ',  i , nsum(i)
  END DO

          !hhg
  OPEN(nion,file=ion_file)          ! ionization file
  DO i = 1, lmax+1
     CALL tafile(i,i-1,"pop") ! partial population files
  ENDDO
  OPEN(npop,file=pop_file)          ! population   file
  IF(hohg) THEN

     OPEN(ndip, file=dip_file)          ! rt(t) file

     WRITE(ndip,'(a1,a60,a1)') '&',' gauge = ', gauge
     WRITE(ndip,'(a1,8a25)')     &
     & '&',                      &
     & 't',                      &
     & 'rt = r(t)',              &
     & 'vt = p(t) + Z*A(t)',     &
     & 'at = d/dt p(t) - Z*E(t)',&
     & 'p(t) = vt-Z*A(t)'       ,&
     & 'd/dt p(t) = at + Z*E_t' ,&
     & 'A(t)'                   ,&
     & 'E(t) = -dA(t)/dt)'

  ENDIF

  WRITE(*,'(a15,2x,a15)') 't','norm'
  PRINT*, nmax, ntot, neq
     ALLOCATE(    dpr(nmax) )
     ALLOCATE(    dpi(nmax) )

  tSteps = 20
  IF(hohg) THEN
     tSteps = 128
     ALLOCATE(   time(nfft) )
     ALLOCATE(     rt(nfft) )
     ALLOCATE(     at(nfft) )
     ALLOCATE(     vt(nfft) )
  ENDIF



  ALLOCATE( yderiv(neq)  )            ! dy/dt
  ALLOCATE( w(neq, 7) )              ! work array for nag
!  ALLOCATE(work(100+21*neq))         ! work array for ode

  pop_g = 0.0_dpk
  pop_b = 0.0_dpk
  pop_k = 0.0_dpk

  dim_t = cycles
  ALLOCATE(    tm( dim_t)       )
  ALLOCATE( ct_re( dim_t, ntot) )
  ALLOCATE( ct_im( dim_t, ntot) )


  it_index = 1
  nt_out = cycles*tSteps/dim_t      ! 2
  i_out  = 1

  WRITE(*,*) it_index, i_out, nt_out


  kt     = 0
  time_loop:DO  l = lstart, cycles !* tSteps


     ! break more steps to calculate the time-dependent rt
     nout1 = kt + 1
     time_loop_for_hhg: DO  lk = 1, tSteps

        kt  = kt + 1
        tim = t + DBLE(tp/tsteps)

        !               WRITE(*,*)  t, tim, DBLE(Tp/tSteps)

        neval = 0
        tol   = 1.0D-10

        !ifail = 1
        !CALL ode ( fcn, 2*ntot, y, t, tim, tol, tol, ifail, work, iwork )

        ifail = 0

        IF(gauge.EQ.'v') THEN
           CALL d02baf(t, tim, neq, y, tol, tdse_v, w, ifail)
        ELSE IF(gauge.EQ.'l') THEN
           CALL d02baf(t, tim, neq, y, tol, tdse_l, w, ifail)
        ELSE IF(gauge.EQ.'lv') THEN
           CALL d02baf(t, tim, neq, y, tol, tdse_lv, w, ifail)
        ENDIF

        !                WRITE(*,*)  t, tim
        !                stop
        !                t   = tim

        !   O(t) = < psi(t) | O | psi(t) > = 2 * Re ( Sum_L <C_L(t) | O | C_{L+1}(t) > )
        !
        !        = Sum_L [ < C^R_L(t) | O(r) | C^R_(L+1)(t) >  + < C^I_L(t) | O | C^I_(L+1)(t) > ]
        !

        harmonics:IF(hohg) THEN



           time(kt)   = tim

           !
           !   O(r) = z == z_1 + z_2
           !
           !   z(t) = < psi(t) | z | psi(t) >   (psi(t) = length or velocity )
           !

           !
           ! position vector             (dr matrix)
           !

           !note for the 'lv' option dr(n,n') --> dr(n,n')/( w_n - w_n')


           rt(kt) = 0.0_dpk
           loop_rt:DO li = 1, lmax

              ! D_(LL+1) * C^R_{L+1}
              ! D_(LL+1) * C^I_{L+1}

              CALL dgemv('t', n(li+1), n(li), 1.0d0, dr(1,1,li), nmax, y(  0  + nsum(li+1) ), 1, 0.0d0, dpr, 1 )
              CALL dgemv('t', n(li+1), n(li), 1.0d0, dr(1,1,li), nmax, y(ntot + nsum(li+1) ), 1, 0.0d0, dpi, 1 )

              rt(kt) =   rt(kt) + &
                          &   DOT_PRODUCT( y(        nsum(li):        nsum(li) + n(li) - 1 ), dpr) & ! <C_R|D|C_R>
                          & + DOT_PRODUCT( y( ntot + nsum(li): ntot + nsum(li) + n(li) - 1 ), dpi)   ! <C_I|D|C_I>

           END DO loop_rt

           rt(kt) = 2.0_dpk * rt(kt)

!
!         pz = p_z_1 + p_z_2
!
!                  d                       dH
!      < vz(t) > =  --<z> = -i [z, H(t)] = --
!                  dt                      dpz
!
!              = < pz - Z*q*A(t) >  = < psi(t) | pz | psi(t) > - Z*q*A(t)
!

           !
           !mechanical momentum          (dzr matrix)
           !

           vt(kt) = 0.0_dpk
           loop_vt:DO li = 1, lmax

              ! P_(LL+1) * C^R_{L+1}
              ! P_(LL+1) * C^I_{L+1}

              CALL dgemv('t', n(li+1), n(li), 1.0d0, dzr(1,1,li), nmax, y(  0  + nsum(li+1) ), 1, 0.0d0, dpr, 1 )
              CALL dgemv('t', n(li+1), n(li), 1.0d0, dzr(1,1,li), nmax, y(ntot + nsum(li+1) ), 1, 0.0d0, dpi, 1 )

              vt(kt) =  vt(kt)     &
                   &      + DOT_PRODUCT( y(    0 + nsum(li):   0  + nsum(li) + n(li) - 1 ), dpi)  & ! <C_R|P|C_I>
                   &      - DOT_PRODUCT( y( ntot + nsum(li): ntot + nsum(li) + n(li) - 1 ), dpr)    ! <C_I|P|C_R>


           END DO loop_vt


           vt(kt) = 2.0_dpk * vt(kt) + A_t( time(kt) ) !! + 2.0D+00*A_t(time(ikt))


!
!
!             d                       d
!      a(t) = --  [ p - Z*q*A(t) ]  = --- < psi(t) | p | psi(t) > + Z*q*E(t)
!             dt                      dt
!
!           = < dpsi/dt | p | dpsi/dt >  + < psi(t) | p | dpsi/dt > + Z*q*E(t)
!
!
!          or
!             d                                               dH
!      a(t) = --<p> + Z*q*E(t) = -i [p, H] + Z*q*E(t) = -i(-) -- + Z*q*E(t)
!             dt                                              dr
!                dV
!           = - --- + Z * q * E(t)
!                dr
!


!
! 'mechanical' acceleration         (dpr matrix)
!

           at(kt) = 0.0_dpk
           loop_at:DO li = 1, lmax

              !<f|p|df/dt>!
              CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,yderiv(    nsum(li+1)), 1,0.0d0, dpr, 1)
              CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,yderiv(ntot+nsum(li+1)),1,0.0d0, dpi, 1)

!              at(kt) =   & at(kt)
!                         & + DOT_PRODUCT( y(nsum(li):nsum(li) + n(li) - 1 ), dpr)           &
!                         & + DOT_PRODUCT( y(ntot+nsum(li):ntot+nsum(li) + n(li) - 1 ), dpi)

              at(kt) =   at(kt) &
                         & - DOT_PRODUCT( y(ntot+nsum(li):ntot+nsum(li) + n(li) - 1 ), dpr)  &
                         & + DOT_PRODUCT( y(nsum(li):+nsum(li) + n(li) - 1 ),          dpi)

              !<df/dt|p|f>!
              CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,y(nsum(li+1)),1,0.0d0,dpr,1)
              CALL dgemv('t',n(li+1),n(li),1.0d0,dzr(1,1,li),nmax,y(ntot+nsum(li+1)),1,0.0d0,dpi,1)

              at(kt) = at(kt) &
                   & - DOT_PRODUCT( yderiv( ntot+nsum(li): ntot + nsum(li) + n(li) - 1 ), dpr) &
                   & + DOT_PRODUCT( yderiv(   0 +nsum(li):   0  + nsum(li) + n(li) - 1 ), dpi)

           END DO loop_at

              at(kt) = at(kt) - E_t(time(kt))               !- 2.0D+00*E_t(time(kt))

           ENDIF harmonics


             ! ground state population

             pop_g = y(1)**2 + y(ntot+1)**2

             ! bound state population

             pop_b = 0.0_dpk
             nt = 0
             bound_population_t: DO  k = 1, lmax+1
                DO  k1 = 1, n(k)
                   get_bound_states:IF( en(k1,k).LE.0.0D+00) THEN
                      pop_b = pop_b + ( y(k1 + nt)**2 + y( ntot+nt+k1)**2 )
                   ENDIF get_bound_states
                ENDDO
                nt = nt + n(k)
             ENDDO bound_population_t

             ! ionization
             pop_k = 1.0_dpk - pop_b

             WRITE(nion,'(5E25.14)') tim,  pop_g, pop_b, pop_k, 1-pop_g

             nt = 0
             partial_population:DO  k = 1, lmax+1
                pop_b = 0.0_dpk
                partial_bound_population:DO  k1 = 1, n(k)
                   pop_b = pop_b + ( y(k1 + nt)**2 + y( ntot+nt+k1)**2 )
                ENDDO partial_bound_population

                WRITE(k,'(6e25.14)') tim, pop_b, ( y(k1+nt)**2 + y(ntot+nt+k1)**2 , k1=1,4)
                nt = nt + n(k)
             ENDDO partial_population


          ENDDO time_loop_for_hhg
          !
          !
          !
          !                         rt = < psi_V(t) | r | psi_V(t) >
          ! note that:              vt = d/dt(rt) = p - Z * q_e * A(t)
          !                         at = d/dt(vt) = d/dt ( p )  + Z * q_e * E(t)
          !                                       = - Z * dV/dr + Z * q_e * E(t)
          !    (helium Z = 2 , q = -1)
          !
          print_harmonics: IF(hohg) THEN

             ! check the below:
             ! above nout1=kt+1 and below: nout2 = kt. So, the below will not be executed?
             !

             nout2 = kt
             DO ikt = nout1, nout2
                !                WRITE(ndip,'(2x,4e25.12)') time(ikt), rt(ikt), at(ikt), E_t(time(ikt))
!                WRITE(ndip,'(2x,5e25.12)') time(ikt), at(ikt), vt(ikt),rt(ikt), E_t(time(ikt))
                WRITE(ndip,'(2x,8e25.12)') time(ikt),          rt(ikt), &                  !r(t)
                     &                                                     vt(ikt), &                  !v(t) (in vel gauge)
                     &                                                     at(ikt), &                  !a(t) (in vel gauge)
                     &                                                     vt(ikt) - A_t(time(ikt)), & !v(t) (in len gauge)
                     &                                                     at(ikt) + E_t(time(ikt)), & !a(t) (in len gauge)
                     &                                                     A_t(time(ikt)), &           !A(t)
         &                                                     E_t(time(ikt))              !E(t)
!            WRITE(ndip,'(2x,5e25.12)') time(ikt), at(ikt), vt(ikt), rt(ikt), E_t(time(ikt))
             END DO
          ENDIF print_harmonics

          !
          !   save for later restart
          !

          save_for_restart:IF(irestart.EQ.1) THEN

             REWIND(ncoe)
             WRITE (ncoe) i0, omega, tau, lmax+1, en_cut

             DO k = 1, lmax+1
                WRITE(ncoe) n(k), ( en(k1,k), k1 = 1, n(k))
             ENDDO

             WRITE(ncoe) ntot

             irestart = 0
          ENDIF save_for_restart


          WRITE(*,*) '# cycle = ', l, tim

          tm(l) = tim
          ct_re(l, 1:ntot) = y(     1: ntot)
          ct_im(l, 1:ntot) = y(ntot+1: 2*ntot)

          WRITE(ncoe) tim, ( y(ii), ii = 1, neq)

          !  write population for partial waves

          WRITE(npop,'(i8, 1X,e15.6)') kt, tim

          pop = 0.0_dpk
          nt = 0
          population_t: DO  k = 1, lmax+1

             !  |Y(nL)|^2  , n = 1-5, L = 0, Lmax

             pop_l = 0.0_dpk
             DO  k1 = 1, n(k)
                pop_l = pop_l + y(k1 + nt)**2 + y( ntot + nt + k1 )**2
             END DO

             pop = pop + pop_l

             WRITE(npop,'(I6,1X,5e25.14)') k-1, pop_l, ( y(k1+nt)**2 + y(ntot+nt+k1)**2 , k1=1,4)

             !1             DO k1 = 1, 5
             !1               WRITE(npop,'(3e15.6)') en(k1,k),
             !1             ENDDO


             nt = nt + n(k)

          ENDDO population_t


          WRITE(*,*) tim, pop


          !          WRITE(npop,'(/)')

          nt = 0
          pop = 0.0_dpk

          !                   t   = tim

       ENDDO time_loop

       CLOSE(ncoe)
       CLOSE(npop)


       CALL save_ct_ncFile(ncid, tm, ct_re, "tr", "re")
       CALL save_ct_ncFile(ncid, tm, ct_im, "ti", "im")
       CALL close_tdse_ncFile(ncid)




!store for spectra analysis


       OPEN(npes,file=pes_file)
       WRITE(*,*) '# pes for partial waves l = 0 - ', lmax+1
       WRITE(*,*) '#                     stored in ', pes_file


       WRITE(npes,'(a1,1x,2i6)') '&', lmax+1, ntot
       ntot_l = 0
       partial_wave_l:DO  l = 0, lmax

          WRITE(npes,'(a1,1x,2i6)') '&',l, ntot_l   ! state 'n' of partial wave 'l'

          eigenstate_ne:DO  ne = 1, n(l+1)
             WRITE(npes,'(4e14.6)')&
                  & en( ne, l+1)                                    , &
                  &  y( ne + ntot_l )**2 + y(ne + ntot_l + ntot)**2 , &
                  &  y( ne + ntot_l )                               , &
                  &  y( ne + ntot_l + ntot )
          ENDDO eigenstate_ne

          ntot_l = ntot_l + n(l+1)

       ENDDO partial_wave_l

       CLOSE(npes)

       ! format statements

     END PROGRAM tdse_bs_fxd_2e

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!EOF
