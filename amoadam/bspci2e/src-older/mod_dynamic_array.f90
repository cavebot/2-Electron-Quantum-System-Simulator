      !****
      !
      !  Since we can't use an allocatable array as an argument
      !  to a function or a subroutine, we create a module that
      !  allows this functionality.
      !
      !****

MODULE DynamicIntegerArray
  INTEGER                               IArraySize
  INTEGER, DIMENSION(:), ALLOCATABLE :: Iarray
CONTAINS
  !****
      !****
      !
      !  Increase the size of the dynamically allocated array
      !  by one. Keep all of the old values, and store a
      !  value of zero in the new element
      !
      !EOF!
  SUBROUTINE ResizeIArray    
    !    USE DynamicIntegerArray    
    IMPLICIT NONE
    !
    INTEGER, DIMENSION(:), ALLOCATABLE :: Temp
    INTEGER                            :: AllocateStatus
    INTEGER                            :: i

    
    IF (IArraySize < 0) STOP "# IArraySize < 0 #"
    
    add_one_element:IF (IArraySize == 0) THEN
       IArraySize = 1
       ALLOCATE(Iarray(1),STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "# 1: Allocate failed #"

    ELSE

       ALLOCATE(Temp(IArraySize), STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "# 2: Allocate failed #"
       
       ! temp = Iarray
       
       DO i = 1, IArraySize
          Temp(i) = Iarray(i)
       END DO

       CALL DeallocIArray

       IArraySize = SIZE(Temp) + 1

       ALLOCATE(Iarray(IArraySize), STAT = AllocateStatus )
       IF (AllocateStatus /= 0) STOP "# 3: Allocate failed  #"

       DO i=1,SIZE(Temp)
          Iarray(i) = Temp(i)
       END DO
       
    END IF add_one_element


    Iarray(IArraySize) = 0  ! init to zero last (added) element


    END SUBROUTINE ResizeIArray

      !****
      ! Subroutine to deallocate the storage for the dynamically
      ! allocated array in the DynamicIntegerArray module
      !****

    SUBROUTINE DeallocIArray
      !      USE DynamicIntegerArray
      IMPLICIT NONE
      !
      INTEGER    :: DeAllocateStatus
      !
      IF (ALLOCATED(Iarray)) THEN
         DEALLOCATE(Iarray, STAT = DeAllocateStatus)
         IF (DeAllocateStatus /= 0) &
              STOP "# DeallocIArray: trouble deallocating #"
      END IF;
      IArraySize = 0
    END SUBROUTINE DeallocIArray
    
      
  END MODULE DynamicIntegerArray

  

  !****
  ! The program calls a subroutine ResizeIarray to
  ! increase the size of a dynamically allocated array
  ! by one and retain all of its previous contents
  !****
  
!!%      PROGRAM TestDynamicAllocationFunctionArguments
!!%      USE DynamicIntegerArray
!!%      IMPLICIT NONE
!!%
!!%!      INTERFACE 
!!%!          SUBROUTINE ResizeIArray
!!%!          END SUBROUTINE ResizeIArray
!!%!      END INTERFACE
!!%
!!%
!!%!      INTERFACE 
!!%!         SUBROUTINE DeallocIArray
!!%!         END SUBROUTINE DeallocIArray
!!%!      END INTERFACE
!!%
!!%      INTEGER :: ix
!!%      INTEGER :: i
!!%
!!%      IArraySize = 0 ! Initialize dynamic array
!!%
!!%      DO i=1,5
!!%         WRITE(*,'(A)', ADVANCE = 'NO') "Enter an integer: "
!!%         READ *, ix
!!%         CALL ResizeIArray
!!%         Iarray(i) = ix
!!%         WRITE(*,'("Stored ", I0, " in Iarray(", I0, ")"/)') ix, i
!!%      END DO
!!%      DO i=1,5
!!%         WRITE(*, '("ix(",I0,") = ", I0)') i, Iarray(i)
!!%      END DO
!!%      
!!%      CALL DeallocIArray
!!%      
!!%    END PROGRAM TestDynamicAllocationFunctionArguments

