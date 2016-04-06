MODULE buffer_types

  USE kinds, ONLY: dp
#include "./base/base_uses.f90"

  IMPLICIT NONE

  INTEGER, PRIVATE, SAVE :: last_buffer_id = 0

  TYPE buffer_data
     INTEGER :: id_nr, ref_count
     INTEGER :: length
     INTEGER, DIMENSION(:), ALLOCATABLE :: disps
     INTEGER, DIMENSION(:), POINTER :: data_i
     REAL(KIND=dp), DIMENSION(:), POINTER :: data_d
     COMPLEX(KIND=dp), DIMENSION(:), POINTER :: data_z
  END type buffer_data

  TYPE buffer_obj
     TYPE(buffer_data), POINTER, PRIVATE :: obj => NULL()
  END TYPE buffer_obj
  
CONTAINS
  
  SUBROUTINE buffer_retain(buffer)
    TYPE(buffer_obj), INTENT(INOUT) :: buffer
    CPASSERT(ASSOCIATED(buffer%obj))
    buffer%obj%ref_count = buffer%obj%ref_count + 1
  END SUBROUTINE buffer_retain
  
  PURE FUNCTION is_same_obj(A, B) RESULT(res)
    TYPE(buffer_type)
  END FUNCTION is_same_obj

END MODULE buffer_types
