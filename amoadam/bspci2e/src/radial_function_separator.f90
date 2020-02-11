!---25/01/2019---
!Program calculates the 2-e Radial Distribution in terms of
!the 1-e radial wavefunctions, configuaration interaction coefficients,
!and TDSE coefficients according (Reference to be added)

program ad2e_test

    use precision


    implicit none

    real(DPK) :: dummie1, dummie2, dummie3, dummie4 ,dummie5
    INTEGER   :: ne_1e, nl_1e
    INTEGER   :: l1e_max, dummie6, k


!Redundant variables now after reformulating loops. Fix this.
    integer                 :: counter !for checking progress
    integer                 :: i, ijmax                !Index for energies
    integer                 :: l, l1max                !Index for angular momenta
    integer                 :: nof_points, nof_points2, nof_points_angular

    !Wavefunction variables
    integer                 :: ir, ir1, ir2!, ir_skip            Wavefunction point index, i.e piont 1, 2, etc
    real(DPK), dimension(:), allocatable :: r          !Wavfefunction grid values
    real(DPK), dimension(:,:,:), allocatable:: p       !Wavefunction values, P_nl(r)



    character(len=32) :: wf1efilename, wf1efilename2             !Wavefunction data file
    logical :: file_exists


    integer :: ki, R1, R2, l2,k2, j



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
    write(*,*) "#number of wavefunction points = ", 2000

!-------Read wavefunctions to p(np,ijmax,l1max). Tested 28/01/2019. Correct to this point--------



    nof_points = 2000

    allocate (r(1:nof_points) )
    allocate( p(1:nof_points,1:ijmax,1:l1max) )

      loop_partial_waves: do l = 0, 3

            write (wf1efilename, '(a,I1,a)') 'out/1ewavefunctions-',l,'.out'
            open(16, file = wf1efilename, status = 'old')

            loop_eigenstates: do i = 1, ijmax

              write (wf1efilename2, '(a,I1,a,I2,a)') 'out/radial_wfs/radial-', l, '-', i,'.out'
              open(19, file = wf1efilename2, status = 'replace')

                loop_wavefunction_grid: do ir = 1, 2000

                      read(16,'(3E20.10)') r(ir) , p(ir, i, l+1)
                      write(19,'(3E20.10)') r(ir), p(ir, i, l+1)

                end do loop_wavefunction_grid

                close(19)

             end do loop_eigenstates
         close(16)
     end do loop_partial_waves


end program ad2e_test
