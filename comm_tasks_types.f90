MODULE mat_com_task_types

  IMPLICIT NONE

  PRIVATE
  
  ! public parameters:
  PUBLIC :: TASK_LEN,  &
            TASK_DEST, &
            TASK_SRC,  &
            TASK_IROW, &
            TASK_ICOL, &
            TASK_COST

  INTEGER, PARAMETER :: TASK_LEN = 5
  INTEGER, PARAMETER :: TASK_LEN = 5


  TYPE mat_com_task_data
     INTEGER :: id_nr, ref_count
     INTEGER, DIMENSION(:), ALLOCATABLE :: tasks
     INTEGER :: ntasks
  END type mat_com_task_data

  


CONTAINS
  

END MODULE mat_com_task_types
