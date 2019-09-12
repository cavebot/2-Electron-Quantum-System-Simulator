!
!  subroutines of general use
!
! check_real_input: checks agains input files (real inputs)
! check_inte_input: checks agains input files (integer inputs)
!
MODULE utils
  USE PRECISION, only:dpk
  IMPLICIT NONE
  public
  contains

       SUBROUTINE check_real_input(a,b)
         USE PRECISION,ONLY:dpk
         REAL(dpk) a,b
         !
         IF(ABS(a-b)>1.0d-07) THEN
            
            WRITE(*,*) '# check_input_real:: incosistent input:'
            WRITE(*,*) '# check_input_real::              a = ', a
            WRITE(*,*) '# check_input_real::              b = ', b
            WRITE(*,*) '# check_input_rela:: exiting.'
            STOP
         ENDIF
       END SUBROUTINE check_real_input

       SUBROUTINE check_int_input(i,j)
         INTEGER i,j
         !
         IF(i.ne.j) THEN
            
            WRITE(*,*) '# check_input_int:: incosistent input:'
            WRITE(*,*) '# check_input_int::              i = ', i
            WRITE(*,*) '# check_input_int::              j = ', j
            WRITE(*,*) '# check_input_int:: exiting.'
            STOP
         ENDIF         
       END SUBROUTINE check_int_input

END MODULE utils
