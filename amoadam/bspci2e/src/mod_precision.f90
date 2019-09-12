

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
MODULE PRECISION

  IMPLICIT NONE
  PUBLIC
!  INTEGER, PARAMETER :: dpk  = KIND(1.0D+00) 
!  INTEGER, PARAMETER :: Long =  KIND(1.0D+00) 
  INTEGER, PARAMETER :: dpk = SELECTED_REAL_KIND (15, 300)
  INTEGER, PARAMETER :: Long = Selected_Real_Kind (15, 300)
  INTEGER, PARAMETER :: Longest = SELECTED_REAL_KIND (15, 300) ! if 64 bit is all that's available
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(8)
  !  INTEGER, PARAMETER :: Long = KIND(1.0D+00) 



  !  INTEGER(KIND=i8) :: ia
  !  REAL(KIND=i10) :: a
  !  PRINT *, HUGE(ia), KIN
  !  PRINT *, RANGE(a), PRECISION(a), KIND(a)


END MODULE PRECISION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



