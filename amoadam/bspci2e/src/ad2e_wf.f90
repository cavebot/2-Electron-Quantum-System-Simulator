!---25/01/2019---
!Program calculates the 2-e Radial Distribution in terms of
!the 1-e radial wavefunctions, configuaration interaction coefficients,
!and TDSE coefficients according (Reference to be added)

program ad2e_radial_integral
  implicit none

    use parameter_tdse_fxd
    use bs_frb_2e,      only: read_target_energies
    use io
    use param,          only: np  !Number of wavefunction grid points
    use precision
    use atom_2e                 !Provides routines for reading CI/TDSE coefficients
    use anglib                  !module for clebsch-gordon Coefficients
    use spherical_harmonics     !module for spherical harmonics



    real(DPK) :: dummie1, dummie2, dummie3, dummie4 ,dummie5
    integer   :: dummie6, k
    real(DPK) :: pi

    integer   :: ll1,ll2, nn1, nn2, nof_l_ncf, ll1check, ll2check
    real(DPK) :: e1, e2 !single electron energies
    INTEGER   :: ne_1e, nl_1e
    INTEGER   :: l1e_max
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
    complex(DPK), dimension(:,:), allocatable:: work_A !bipolar spherical harmonic
    real(DPK), dimension(:,:), allocatable:: y !bipolar spherical harmonic
    real(DPK), dimension(:,:)  , allocatable:: psi2
    real(DPK) :: phi1, phi2, theta1r, theta2r
    integer :: m1, m2, theta1, theta2
    real(DPK) :: pop_ic
    complex(DPK):: pop_t_ie
    complex(DPK), dimension(:,:,:), allocatable:: f, work !2e radial distribution


    character(len=25) :: wf1efilename              !Wavefunction data file
    logical :: file_exists
    integer:: t1, t2, clock_rate, clock_max

!command line arguments
    character(len=100):: partial_l1, partial_l2, partial_l3, partial_l4, num_cycles ,rd2efilename
    integer:: l1_int, l2_int, l1_int_p, l2_int_p, num_cycles_int!, L_int


!get cmd arguments
!    call getarg(1, L_total)
    call getarg(1, partial_l1)
    call getarg(2, partial_l2)
    call getarg(3, partial_l3)
    call getarg(4, partial_l4)
    call getarg(5, num_cycles)
    read(partial_l1(1:1), "(i1)") l1_int
    read(partial_l2(1:1), "(i1)") l2_int
    read(partial_l3(1:1), "(i1)") l1_int_p
    read(partial_l4(1:1), "(i1)") l2_int_p
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



    nof_points = 100
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





          counter  = 0
          !ijmax = 60
          !lmax = 7
          !Nmax = 5998 Need to import this from somewhere

          write(*,*) "#Number of L=0 configurations in sum/prime sums = ", w2e(0)%ncf
          write(*,*) "#Number of L=1 configurations in sum/prime sums = ", w2e(1)%ncf
          write(*,*) "#Number of L=2 configurations in sum/prime sums = ", w2e(2)%ncf
          write(*,*) "#Number of L=3 configurations in sum/prime sums = ", w2e(3)%ncf




!Find ionisation thresholds


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


       !-----------------------------------Apply summation equation-----------------------------------



       !Initialise 2e radial distribution matrix, prr = p(r1,r2)
       allocate(  f(0:7, 0:7, 0:7) )
       allocate(  work(0:7, 0:7, 0:7) )
       f = 0.0_dpk


        call system_clock ( t1, clock_rate, clock_max )

        loop_L: do l = 0, 3

          !if ((l .eq. 1) .or. (l .eq. 2).or. (l .eq. 3)) then
          !  cycle
          !end if
          !write(*,*) "Angularnof_points Symmetry L = ", l

          !Set number of configurations in symmetry L
          nof_l_ncf = w2e(l)%ncf
          write(*,*) "#nof configs in L = ", nof_l_ncf

          !Set maximum number of 2e energies for symmetry L
          ie_max = w2e(l)%net
          write(*,*) "#ie max = ", ie_max

          !LOOP HERE
          ir1 = 50
          ir2 = 50


                work = 0.0_dpk

                loop_total_energy: do ie = ndi_l(l), ie_max           !write(*,*) ie

                  pop_t_ie = w2e(l)%ct(ie)                            !TDSE coefficients

                  loop_config_variables: do ic = 1, nof_l_ncf         !sum from configuration 1 to max number of configurations, for the total  L
                    pop_ic = w2e(l)%cv(ie,ic)                         !Configuration interaction coefficients
                    ll1 = w2e(l)%l1(ic)
                    ll2 = w2e(l)%l2(ic)

                    !if(((ll1) .ne. (l1_int)) .or. ((ll2) .ne. (l2_int))) then
                      !write(*,*) ll1, ll2
                    !  cycle loop_config_variables
                    !end if

                    nn1 = w2e(l)%n1(ic)
                    nn2 = w2e(l)%n2(ic)
                    e1  = e1e(ll1, nn1 )
                    e2  = e1e(ll2, nn2 )

                    !write(*,*)  nn1, ll1, ";", ll1, nn2

                    if((e1 > 0) .and. (e2 > 0)) then

                      if(((ll1) .eq. (ll2)) .and. ((nn1) .eq. (nn2)))then

                        work(ll1, ll2, l) = work(ll1, ll2, l) + p( ir1, nn1,  ll1+1)*p( ir2, nn2,  ll2+1)*pop_ic*pop_t_ie

                      else
                        work(ll1, ll2, l) = work(ll1, ll2, l) + (1/sqrt(2.0))*((p( ir1, nn1,  ll1+1 )*p( ir2, nn2,  ll2+1 )) + &
                        (p( ir2, nn1,  ll1+1 )*p( ir1, nn2,  ll2+1 )))*pop_ic*pop_t_ie

                      end if
                    end if

                  end do loop_config_variables
                end do loop_total_energy

                !write(*,*) f(ll1, ll2), work(ll1, ll2)
                f = f  + work
                !write(*,*) f(ll1, ll2), work(ll1, ll2)
          end do loop_L


          !!!!!!!!!!!!!!!!!!!!!Calculate bipolar spherical harmonics


          nof_points_angular = 200
          allocate( work_A(1:nof_points_angular, 1:nof_points_angular))
          allocate( y(1:nof_points_angular, 1:nof_points_angular))

         !allocate( work_A(1:8, 1:8, 1:4, 1:nof_points_angular, 1:nof_points_angular))

          y = 0.0_dpk
          work_A = 0.0_dpk

          phi1 = 0
          phi2 = 0

          m1 = 0
          m2 = 0

          pi = 4.D0*DATAN(1.D0)


          do l = 0, 3

            nof_l_ncf = w2e(l)%ncf
            k = 0

            loop_l1l2: do ic = 2, nof_l_ncf         !sum from configuration 1 to max number of configurations, for the total  L

              if (k .ge. 5) then
                exit
              end if

              ll1 = w2e(l)%l1(ic)
              ll2 = w2e(l)%l2(ic)

              ll1check = w2e(l)%l1(ic+1)
              ll2check = w2e(l)%l2(ic+1)

              if((ll1check .eq. ll1) .and. (ll2check .eq. ll2)) then
                cycle
              end if



              !!!!!!!!!!!!!!!
              write(*,*) "working", ll1, ll2
              theta1r = 0
              theta2r = 0
              angle_step = 2*pi/nof_points_angular

              do theta1 = 1, nof_points_angular
                do theta2 = 1, nof_points_angular


                  work_A(theta1, theta2) =   work_A(theta1, theta2) + &
                  sphharm(ll1, m1, theta1r, phi1)*sphharm(ll2, m2, theta2r, phi2)!*cleb(2*ll1, 2*m1, 2*ll2, 2*m2, 2*l, 0)!*&
                  !f(ll1, ll2, l)

                  theta2r = theta2r + angle_step
                end do

                theta1r = theta1r + angle_step
              end do
              write(*,*) "finished", ll1, ll2

              !!!!!!!!!!!!!!                  work_A(theta1, theta2) =   work_A(theta1, theta2) + &


              k = k + 1
              ll1 = ll1check
              ll2 = ll2check

            end do loop_l1l2

            write(*,*)

          end do

!!!!MOD SQUARE wavefucntion
          do theta1 = 1, nof_points_angular
            do theta2 = 2, nof_points_angular

              y(theta1, theta2) = abs(work_A(theta1, theta2))**2

            end do
          end do







          call system_clock ( t2, clock_rate, clock_max )
          write ( *, * ) 'Elapsed real time = ', real ( t2 - t1 ) / real ( clock_rate )


          write(rd2efilename,"(a)") "ad2e/ad2e.out"
          rd2efilename  = trim(rd2efilename)//trim(num_cycles)
          open(19, file = rd2efilename, status = 'replace')





                   !write to file

                           do theta1 = 1, nof_points_angular
                               do theta2 = 1, nof_points_angular

                                     write(19,'(3E20.10)')  pi/theta1, pi/theta2, y(theta1, theta2)

                               end do
                               write(19,'(3E20.10)')
                           end do

         !Mirror across diagonal and write to file

          !       do ir1 = 1, nof_points2
          !           do ir2 = 1, nof_points2
          !             if ((ir2 .ge. ir1)) then
          !               if (f(ir1,ir2) .lt. 0) then
          !                 write(19,'(3E20.10)')  r(ir1)*2.645, r(ir2)*2.645, 0.00000 !SCALED.ATOMIC.UNITS TO nm
          !               else               if (f(ir2,ir1) .lt. 0) then

          !                 write(19,'(3E20.10)')  r(ir1)*2.645, r(ir2)*2.645, f(ir1,ir2) !SCALED.ATOMIC.UNITS TO nm
          !               end if
          !             else
          !               if (f(ir2,ir1) .lt. 0) then
          !                 write(19,'(3E20.10)')  r(ir1)*2.645, r(ir2)*2.645, 0.00000 !SCALED.ATOMIC.UNITS TO nm
          !               else
          !                 write(19,'(3E20.10)')  r(ir1)*2.645, r(ir2)*2.645, f(ir2,ir1) !SCALED.ATOMIC.UNITS TO nm
          !               end if
          !             end if
          !           end do
          !           write(19,'(3E20.10)')
          !       end do
         close(19)



end program ad2e_radial_integral
