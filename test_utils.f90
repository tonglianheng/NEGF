MODULE test_utils
  USE kinds
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: check_pass, &
            test_finalise, &
            test_init

  LOGICAL, SAVE :: overall_pass
  INTEGER, SAVE :: overall_counter, pass_counter, fail_counter
  CHARACTER(LEN=*), PARAMETER :: PASS_STRING = "[0;32mPASS[0m", &
                                 FAIL_STRING = "[0;31mFAIL[0m"
  REAL(KIND=dp) :: tolerance

  INTERFACE check_pass
     MODULE PROCEDURE check_pass_i
     MODULE PROCEDURE check_pass_d
     MODULE PROCEDURE check_pass_z
  END INTERFACE check_pass

CONTAINS

  SUBROUTINE test_init(tol)
    REAL(KIND=dp), OPTIONAL :: tol
    overall_pass = .TRUE.
    overall_counter = 0
    pass_counter = 0
    fail_counter = 0
    IF (PRESENT(tol)) THEN
       tolerance = tol
    ELSE
       tolerance = 1.E-15_dp
    END IF
  END SUBROUTINE test_init

  SUBROUTINE test_finalise()
    WRITE (*,*) "----------------------------------------------------------------------"
    IF (overall_pass) THEN
       WRITE (*,*) "OVERALL TEST: ", PASS_STRING
    ELSE
       WRITE (*,*) "OVERALL TEST: ", FAIL_STRING
    END IF
    WRITE (*,*) "Total number of tests:  ", overall_counter
    WRITE (*,*) "Number of passed tests: ", pass_counter
    WRITE (*,*) "Number of failed tests: ", fail_counter
  END SUBROUTINE test_finalise

  SUBROUTINE check_pass_i(val, correct)
    INTEGER, INTENT(IN) :: val, correct
    overall_counter = overall_counter + 1
    IF (val .EQ. correct) THEN
       WRITE (*,*) PASS_STRING
       pass_counter = pass_counter + 1
    ELSE
       WRITE (*,*) FAIL_STRING
       WRITE (*,*) "Error = ", val - correct
       overall_pass = .FALSE.
       fail_counter = fail_counter + 1
    END IF
  END SUBROUTINE check_pass_i

  SUBROUTINE check_pass_d(val, correct)
    REAL(KIND=dp), INTENT(IN) :: val, correct
    overall_counter = overall_counter + 1
    IF (ABS(val - correct) .LT. tolerance) THEN
       WRITE (*,*) PASS_STRING
       pass_counter = pass_counter + 1
    ELSE
       WRITE (*,*) FAIL_STRING
       WRITE (*,*) "Error = ", val - correct
       overall_pass = .FALSE.
       fail_counter = fail_counter + 1
    END IF
  END SUBROUTINE check_pass_d

  SUBROUTINE check_pass_z(val, correct)
    COMPLEX(KIND=dp), INTENT(IN) :: val, correct
    REAL(KIND=dp) :: error
    overall_counter = overall_counter + 1
    error = SQRT(REAL((val - correct)*CONJG(val - correct),KIND=dp))
    IF (error .LT. tolerance) THEN
       WRITE (*,*) PASS_STRING
       pass_counter = pass_counter + 1
    ELSE
       WRITE (*,*) FAIL_STRING
       WRITE (*,*) "Error = ", error
       overall_pass = .FALSE.
       fail_counter = fail_counter + 1
    END IF
  END SUBROUTINE check_pass_z

END MODULE test_utils
