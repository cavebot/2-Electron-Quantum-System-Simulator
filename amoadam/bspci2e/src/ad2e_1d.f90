!---25/01/2019---
!Program calculates the 2-e Radial Distribution in terms of
!the 1-e radial wavefunctions, configuaration interaction coefficients,
!and TDSE coefficients according (Reference to be added)

program ad2e_test


    use parameter_tdse_fxd
    use bs_frb_2e,      only: read_target_energies
    use io
    use param,          only: np  !Number of wavefunction grid points
    use precision
    use atom_2e                 !Provides routines for reading CI/TDSE coefficients
    use anglib                  !module for clebsch-gordon Coefficients
    use spherical_harmonics     !module for spherical harmonics
    use core_phase_shift

    implicit none

    real(DPK) :: pi
    real(DPK) :: dummie1, dummie2, dummie3, dummie4 ,dummie5

    integer   :: ll1,ll2, nn1, nn2, nof_l_ncf, ll1check, ll2check
    real(DPK) :: e1, e2 !single electron energies
    INTEGER   :: ne_1e, nl_1e
    INTEGER   :: l1e_max, dummie6
    real(dpk), dimension(:,:), pointer :: en1e        !nl,ns
    REAL(dpk), ALLOCATABLE, DIMENSION(:,:)         :: e1e                 !nl,ns
    integer,   allocatable, dimension(:)        :: ndi_l           ! pop with  E(ndi_l,l) = E++
    integer,   allocatable, dimension(:)        :: nsi_l           ! pop with  E(ndi_l,l) = E+
    integer                                     :: n_l2e           ! nof CI states of symmetry l2e
    integer                                     :: cycle  !cycle number (time of evaluation of distribution)

!Redundant variables now after reformulating loops. Fix this.
    integer                 :: counter !for checking progress
    integer                 :: i, ijmax                !Index for energies
    integer                 :: l, l1max                !Index for angular momenta
    integer                 :: ic, ie, ie_max    !,ie_p        !Index for reformulated summations in terms of configuration number
    integer                 :: nof_points, nof_points2, nof_points_angular
    real(dpk)               :: angle_step
    integer :: mod_dummy

    !Wavefunction variables
    type(symmetry_ls2e), allocatable, dimension(:) :: w2e
    integer                 :: ir, ir1, ir2!, ir_skip            Wavefunction point index, i.e piont 1, 2, etc
    real(DPK), dimension(:), allocatable :: r          !Wavfefunction grid values
    real(DPK), dimension(:,:,:), allocatable:: p       !Wavefunction values, P_nl(r)
    complex(DPK), dimension(:,:)  , allocatable:: psi
    complex(DPK), dimension(:), allocatable:: work_A !bipolar spherical harmonic
    real(DPK), dimension(:,:), allocatable:: y !bipolar spherical harmonic
    real(DPK), dimension(:,:)  , allocatable:: psi2
    real(DPK) :: phi1, phi2, theta1r, theta2r
    integer :: m1, m2, theta1, theta2
    real(DPK) :: pop_ic
    complex(DPK):: pop_t_ie
    complex(DPK), dimension(:,:,:), allocatable:: f, work !2e radial distribution


    character(len=32) :: wf1efilename, phasefilename, ad2efilename          !Wavefunction data file
    logical :: file_exists
    integer:: t1, t2, clock_rate, clock_max

    !command line arguments
    character(len=100):: partial_l1, partial_l2, num_cycles ,rd2efilename
    integer:: k1_int, k2_int, num_cycles_int, k1_select, k2_select, k_val

    !For phase shifts
    real(DPK), dimension(:,:,:), allocatable :: phase !delta_kl array for phase shift values
    real(DPK), dimension(:,:), allocatable :: k_value !delta_kl array for phase shift values
    integer :: ki, R1, R2, l2,k2, j, av
    integer :: reduced_ijmax, k !only energies up to roughly n=40 have valid calculated radial functions. Use this instead.

    pi = 4.D0*DATAN(1.D0)



!get cmd arguments
!    call getarg(1, L_total)
    call getarg(1, partial_l1)
    call getarg(2, partial_l2)
    call getarg(3, num_cycles)
    read(partial_l1(1:2), "(i2)") k1_int
    read(partial_l2(1:2), "(i2)") k2_int
    read(num_cycles(1:2), "(i2)") num_cycles_int
    write(*,*) num_cycles_int, "Field cycle number:"
   !read(L_total(1:1), "(i1)") L_int


!----------------------------Get number of states 1e from inp/h1e.inp -------------------------
    open(unit = 8, file = "inp/h1e.inp", status = 'old')
        !Dummie file read
        read(8, '(2f15.5,10x,a10)')           dummie1, dummie2
        read(8, '(a10,2f10.5,10x,a25)')       dummie3, dummie4, dummie5
        !read number of B-splines
        read(8, '(2I15,10x,a8)') ijmax, dummie6
        write(*,*) "#Number of b-splines = ", ijmax
        ijmax = ijmax - 2
        write(*,*) "#Number of eigenstates = ", ijmax
    close(8)



!------------------------Check wavefunction files exist. Retrive maximum angular momentum-----------
    l = 0
    write (wf1efilename, '(a,I1,a)') 'out/1ewavefunctions-',l,'.out'
    inquire(FILE = wf1efilename, EXIST = file_exists)
        if (file_exists) then
        write(*,*) "#Wavefunction files found"
        write(*,*) "#Finding number of partial wave files"
        write(*,*)
    else
        write(*,*) "#Could not find wavefunction files, 1ewavefunction-l.out."
        write(*,*) "#Terminating program"
        stop
    end if

    do while (file_exists)
        write (wf1efilename, '(a,I1,a)') 'out/1ewavefunctions-',l,'.out'
        inquire(FILE = wf1efilename, EXIST = file_exists)
        if(file_exists) then
            l = l + 1
        else
            exit
        end if
    end do

    l1max = l
    write(*,*) "#Number of partial waves = l1max = ", l1max
    write(*,*) "#number of wavefunction points = ", np

!-------Read wavefunctions to p(np,ijmax,l1max). Tested 28/01/2019. Correct to this point--------



    nof_points = 2000
    nof_points2 = 100
    mod_dummy = 2000/nof_points


    write(*,*) 'taking every', mod_dummy, 'th point'

    allocate (r(1:nof_points) )
    allocate( p(1:nof_points,1:ijmax,1:l1max) )

      loop_partial_waves: do l = 0, (l1max-1)

            write (wf1efilename, '(a,I1,a)') 'out/1ewavefunctions-',l,'.out'
            open(16, file = wf1efilename, status = 'old')

            loop_eigenstates: do i = 1, ijmax

              read(16,'(3E20.10)') r(1) , p(1, i, l+1)

                loop_wavefunction_grid: do ir = 2, 2000

                   if (mod(ir, mod_dummy) .eq. 0) then
                      read(16,'(3E20.10)') r(ir/mod_dummy) , p((ir/mod_dummy), i, l+1)
                  !     read(16,'(3E20.10)') r(ir) , p((ir), i, l+1)
                  !     write(*,*) r(ir)
                   else
                      read(16,'(3E20.10)')
                   end if
                end do loop_wavefunction_grid
                !loop_skip: do ir_skip = (nof_points+1), 2000
                !  read(16,'(3E20.10)')
                !end do loop_skip
             end do loop_eigenstates
         close(16)
     end do loop_partial_waves

      !--------------------------R￼ead TDSE coeffs and CI coeffs and configuration data-------------------------

      call input_tdse_fxd              ! tinp/tdse_bs_fxd_2e.inp
      call output_tdse_fxd             ! tout/tdse_bs_fxd_2e.out
      call read_target_energies(en1e)  ! read 1e energies (in Ryd)

      cycle = num_cycles_int
      write(*,*) "Calculating field cycle", cycle


      cycle = 12
      call read_w2e_ct(w2e,lmax,cycle)

      do l = 0, lmax
         call read_w2e_ci(w2e(l),l)     !reads coeffs Phi_j = \sum_i v_{ij} Phi^{(0)}_{i}
      end do





      !!!!!!!!1-electron energy matrix initialisations

      en1e = en1e * 0.5_dpk              ! convert energies in a.u.
      nl_1e = SIZE(en1e,dim=1)  ! e1e(l1_max, n1_max)
      ne_1e = SIZE(en1e,dim=2)  ! e1e(l1_max, n1_max)

      ! for conformity reasons we rewrite
      ! en1e(1:nl_1e,1:ne_1e) --> e1e(0:nl_1e-1,1:ne_1e) -->

      !First find the lread_target_energiesargest angular momentum included in the cfg channels.
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Calculate ionisatation thresholds

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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculate phase shifts

      reduced_ijmax = 40
      av = 100

      allocate(phase(0:3, 1:reduced_ijmax, 1:(av+1)))  !phase shifts stored in this array
      allocate(k_value(0:3, 1:reduced_ijmax))  !array for storing K values corresponding to 1e energies

      phase = 0.0_dpk
      k_value = 0.0_dpk

        loop_l: do l = 0, 3
          loop_k: do ki = 1, reduced_ijmax

            !calculate K values from 1e energies. If bound state energy, assign -1
            if (e1e(l,ki) .lt. 0) then
              k_value(l,ki) = -1.0
            else
              k_value(l,ki) = sqrt(2*e1e(l,ki))
            end if

            loop_average:  do j = 1, av !Loop to calculate average phase shift using different points on P(r)

              R1 = (nof_points-(j+10)) !!Point that phase is evaluated at, i.e close to but not on box boundary
              R2 = (nof_points-(j+20))

                phase(l,ki,j) = (p(R2, ki, (l+1))*bessel_jn(l,((k_value(l,ki))*(r(R1)))) - p(R1, ki, (l+1))*&
                bessel_jn(l,((k_value(l,ki))*(r(R2)))))/ &
                (p(R1, ki, (l+1))*bessel_yn(l,((k_value(l,ki))*(r(R2)))) - p(R2, ki, (l+1))*bessel_yn(l,((k_value(l,ki))*(r(R1)))))

                phase(l,ki,j) = atan(phase(l,ki,j))
                phase(l,ki,(av+1)) = phase(l,ki,(av+1)) + phase(l,ki,j)

            end do loop_average

            phase(l,ki,(av+1)) = phase(l,ki,(av+1))/av


            !if (phase(l,ki,(av+1)) .lt. 0.0) then
            !  phase(l,ki,(av+1)) = phase(l,ki,(av+1)) + 3.14159
            !end if


          !  write(*,*) e1e(l,ki), k_value(l,ki), phase(l,ki,(av+1)) !write energies, k's and phases

          end do loop_k

        !  write(*,*)
        end do loop_l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!save phases

        loop_l2: do l = 0, 3

          write (phasefilename, '(a,I1,a)') 'phaseshift-',l,'.out'
          open(19, file = phasefilename, status = 'replace')

          loop_k2: do ki = 1, reduced_ijmax

            if ( k_value(l,ki) .gt. 0) then
              write(*,*) l, k_value(l,ki), phase(l,ki,(av+1))!¬, (phase(l,ki,(av+1))-((l*pi)/2)+(k_value(l,ki)*r(nof_points)))/pi
              write(19,'(3E20.10)')  k_value(l,ki), phase(l,ki,(av+1))
            end if

          end do loop_k2

          close(19)

        end do loop_l2


        write(*,*) "finished phases"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calcualte angular distribution


        nof_points_angular = 100
        allocate( work_A(1:nof_points_angular))
        allocate( y(1:nof_points_angular, 1:nof_points_angular))

       !allocate( work_A(1:8, 1:8, 1:4, 1:nof_points_angular, 1:nof_points_angular))

        y = 0.0_dpk
        work_A = 0.0_dpk
        phi1 = 0
        phi2 = 0
        m1 = 0
        m2 = 0
        k1_select = k1_int
        k2_select = k2_int
        write(*,*) "working k values", k1_select, k2_select


        loop_total_L: do l = 0, 3

          write(*,*) "Working total angular momentum:", l

          nof_l_ncf = w2e(l)%ncf
          ie_max = w2e(l)%net

          loop_total_energy: do ie = ndi_l(l), ie_max           !write(*,*) ie

            pop_t_ie = w2e(l)%ct(ie)                            !TDSE coefficients

            loop_l1l2: do ic = 2, nof_l_ncf         !sum from configuration 1 to max number of configurations, for the total  L

              nn1 = w2e(l)%n1(ic) !nn1 is energy index. also indexes k values
              nn2 = w2e(l)%n2(ic)
              if (((nn1 .eq. k1_select) .and. (nn2 .eq. k2_select)) .or. ((nn2 .eq. k1_select) .and. (nn1 .eq. k2_select))) then !select only specific k1 and k2

                !write(*,*) nn1, nn2, ll1, ll2
                ll1 = w2e(l)%l1(ic)
                ll2 = w2e(l)%l2(ic)

                pop_ic = w2e(l)%cv(ie,ic)

                theta1r = 0
                theta2r = 0
                angle_step = 2*pi/nof_points_angular


                loop_angle: do theta1 = 1, nof_points_angular

                  work_A(theta1) =  work_A(theta1) + &
                  sphharm(ll1, 0, theta1r, 0.0_dpk)*sphharm(ll2, 0, 0.0_dpk, 0.0_dpk)*cleb(2*ll1, 0, 2*ll2, 0, 2*l, 0)*&
                  pop_t_ie*pop_ic*&!(0.0,-1.0)**(ll1+ll2)*&
                  cdexp((0.0,1.0)*(-(k_value(ll1,nn1)*60)-(k_value(ll2,nn2)*60)))!*(0.0,-1.0)**(ll1+ll2)!*cdexp((0.0,1.0)*(phase(ll1,nn1,(av+1))+phase(ll1,nn1,(av+1))))
  !cdexp((0.0,1.0)*(k_value(ll1,nn1)+k_value(ll2,nn2))*(60))
                  theta1r = theta1r + angle_step

                end do loop_angle

              end if

            end do loop_l1l2

            !write(*,*)

          end do loop_total_energy
        end do loop_total_L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!save to file

        if ((k1_select .lt. 10) .and. (k2_select .lt. 10)) then
          write(rd2efilename,"(a,i1,a,i1,a)") "ad2e/ad2e_1d_k",k1_select, "_", k2_select, ".out"
        else if ((k1_select .ge. 10) .and. (k2_select .ge. 10)) then
          write(rd2efilename,"(a,i2,a,i2,a)") "ad2e/ad2e_1d_k",k1_select, "_", k2_select, ".out"
        else if ((k1_select .ge. 10) .and. (k2_select .lt. 10)) then
          write(rd2efilename,"(a,i2,a,i1,a)") "ad2e/ad2e_1d_k",k1_select, "_", k2_select, ".out"
        else
          write(rd2efilename,"(a,i1,a,i2,a)") "ad2e/ad2e_1d_k",k1_select, "_", k2_select, ".out"
        end if

        open(19, file = rd2efilename, status = 'replace')


        do theta1 = 1, nof_points_angular

            write(19,'(3E20.10)')  ((theta1-1)*angle_step)*(180/pi), abs(work_A(theta1))**2

        end do

       close(19)



deallocate(r)
deallocate(p)
deallocate(phase)

end program ad2e_test
