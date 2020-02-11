!
!Calculates the phase shift relative to the asymptotic form of the radial functions,
!close to the box boundary. Asymptotic form of radial functions are bessel
!functions of the first and second kind.

module core_phase_shift

  use precision

contains


!Calculate phase shift of calculated radial functions by passing them as
!an array into function and using equation tan(dkl)= p2f1-p1f2/p2g1-p1g2

subroutine phase_shift(p, r, phase, energy, n, kmax, lmax)

    implicit none

    integer :: n, kmax, lmax, l, ki, R1, R2, i, l2,k2
    real(DPK) :: x, J, Y
    real(DPK) :: p(:,:,:), r(:), phase(:,:,:), energy(:,:)
    real(DPK) :: k(lmax, kmax)




      open(19, file = "phase-shift.out", status = 'replace')
      loop_l: do l = 1, lmax
        loop_k: do ki = 1, kmax

          if (energy(l,ki) .lt. 0) then
            k(l,ki) = -1.0
          else
            k(l,ki) = sqrt(2*energy(l,ki))
          end if


          loop_average:  do j = 1, 30

            R1 = (n-(j+10)) !!Points that phase is evaluated at, i.e close to but not on box boundary
            R2 = (n-(j+20))

              !phase(ki,l) = p(R1, ki, l)
              phase(ki,l,j) = (p(R2, ki, l)*bessel_jn(l,((k(l,ki))*(R1))) - p(R1, ki, l)*bessel_jn(l,((k(l,ki))*(R2))))/ &
                            (p(R1, ki, l)*bessel_yn(l,((k(l,ki))*(R2))) - p(R2, ki, l)*bessel_yn(l,((k(l,ki))*(R1))))

              phase(ki,l,j) = atan(phase(ki,l,j))
              phase(ki,l,31) = phase(ki,l,31) + phase(ki,l,j)

          end do loop_average

          if (phase(ki,l,31) .lt. 0) then
            phase(ki,l,31) = phase(ki,l,31) + 1.57079632679
          end if

          write(*,*) energy(l,ki), k(l,ki), phase(ki,l,31)/10
          write(19,'(3E20.10)')  k(l,ki), phase(ki,l,31)/10  !write phase as function of k


        end do loop_k

        write(19,'(3E20.10)')
        write(*,*)

      end do loop_l
      close(19)





      !open(19, file = "phase-shift.out", status = 'replace')
      !    loop_l2: do l2 = 1, lmax
      !      loop_k2: do k2 = 1, kmax
      !        write(19,'(3E20.10)') energy(l2,k2), k(l2,k2), phase(k2,l2)  !write phase as function of k
      !      end do loop_k2
      !      write(19,'(3E20.10)')
      !    end do loop_l2
      !close(19)



  end subroutine phase_shift


end module core_phase_shift
