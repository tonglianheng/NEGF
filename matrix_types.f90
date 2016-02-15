MODULE matrix_types

  USE kinds, ONLY: dp
  USE machine, ONLY: default_output_unit

#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! public types
  PUBLIC :: mat_d_obj, &
            mat_z_obj

  ! public methods
  PUBLIC :: mat_create,&
            mat_delete, &
            mat_mult, &
            mat_inv_cholesky, &
            mat_read, &
            mat_write, &
            mat_nrows, &
            mat_ncols

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'matrix_types'

  TYPE mat_d_obj
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: p => NULL()
  END TYPE mat_d_obj

  TYPE mat_z_obj
     COMPLEX(KIND=dp), DIMENSION(:,:), POINTER :: p => NULL()
  END TYPE mat_z_obj

  INTERFACE mat_create
     MODULE PROCEDURE mat_create_d
     MODULE PROCEDURE mat_create_z
  END INTERFACE mat_create

  INTERFACE mat_delete
     MODULE PROCEDURE mat_delete_d
     MODULE PROCEDURE mat_delete_z
  END INTERFACE mat_delete

  INTERFACE mat_mult
     MODULE PROCEDURE mat_mult_d
     MODULE PROCEDURE mat_mult_z
  END INTERFACE mat_mult

  INTERFACE mat_inv_cholesky
     MODULE PROCEDURE mat_inv_cholesky_d
     MODULE PROCEDURE mat_inv_cholesky_z
  END INTERFACE mat_inv_cholesky

  INTERFACE mat_read
     MODULE PROCEDURE mat_read_d
     MODULE PROCEDURE mat_read_z
  END INTERFACE mat_read

  INTERFACE mat_write
     MODULE PROCEDURE mat_write_d
     MODULE PROCEDURE mat_write_z
  END INTERFACE mat_write

  INTERFACE mat_nrows
     MODULE PROCEDURE mat_nrows_d
     MODULE PROCEDURE mat_nrows_z
  END INTERFACE mat_nrows

  INTERFACE mat_ncols
     MODULE PROCEDURE mat_ncols_d
     MODULE PROCEDURE mat_ncols_z
  END INTERFACE mat_ncols

CONTAINS

  SUBROUTINE mat_create_d(mat, nrows, ncols)
    TYPE(mat_d_obj), INTENT(INOUT) :: mat
    INTEGER, INTENT(IN) :: nrows, ncols
    ALLOCATE(mat%p(nrows,ncols))
    mat%p = 0.0_dp
  END SUBROUTINE mat_create_d

  SUBROUTINE mat_create_z(mat, nrows, ncols)
    TYPE(mat_z_obj), INTENT(INOUT) :: mat
    INTEGER, INTENT(IN) :: nrows, ncols
    ALLOCATE(mat%p(nrows,ncols))
    mat%p = (0.0_dp,0.0_dp)
  END SUBROUTINE mat_create_z

  SUBROUTINE mat_delete_d(mat)
    TYPE(mat_d_obj), INTENT(INOUT) :: mat
    IF (ASSOCIATED(mat%p)) DEALLOCATE(mat%p)
  END SUBROUTINE mat_delete_d

  SUBROUTINE mat_delete_z(mat)
    TYPE(mat_z_obj), INTENT(INOUT) :: mat
    IF (ASSOCIATED(mat%p)) DEALLOCATE(mat%p)
  END SUBROUTINE mat_delete_z

  SUBROUTINE mat_real_to_complex(mat_d, mat_z)
    TYPE(mat_d_obj), INTENT(IN) :: mat_d
    TYPE(mat_z_obj), INTENT(INOUT) :: mat_z

    INTEGER :: nrows, ncols

    CPASSERT(ASSOCIATED(mat_d%p))
    IF (ASSOCIATED(mat_z%p)) DEALLOCATE(mat_z%p)
    ALLOCATE(mat_z%p(nrows,ncols))
    mat_z%p(:,:) = CMPLX(mat_d%p(:,:), 0.0, KIND=dp)
  END SUBROUTINE mat_real_to_complex

  SUBROUTINE mat_mult_d(transA, transB, alpha, A, B, beta, C)
    TYPE(mat_d_obj), INTENT(IN) :: A, B
    TYPE(mat_d_obj), INTENT(INOUT) :: C
    CHARACTER, INTENT(IN) :: transA, transB
    REAL(KIND=dp), INTENT(IN) :: alpha, beta

    INTEGER :: nrows_A, ncols_A, ncols_B

    nrows_A = SIZE(A%p, 1)
    ncols_A = SIZE(A%p, 2)
    ncols_B = SIZE(B%p, 2)

    CPASSERT(SIZE(B%p,1) .EQ. ncols_A)
    IF (ASSOCIATED(C%p)) THEN
       CPASSERT(SIZE(C%p,1) .EQ. nrows_A)
       CPASSERT(SIZE(C%p,2) .EQ. ncols_B)
    ELSE
       CALL mat_create(C, nrows_A, ncols_B)
    END IF
    CALL DGEMM(transA, transB, &
               nrows_A, ncols_B, ncols_A, &
               alpha, A%p, nrows_A, B%p, ncols_A, &
               beta, C%p, nrows_A)
  END SUBROUTINE mat_mult_d

  SUBROUTINE mat_mult_z(transA, transB, alpha, A, B, beta, C)
    TYPE(mat_z_obj), INTENT(IN) :: A, B
    TYPE(mat_z_obj), INTENT(INOUT) :: C
    CHARACTER, INTENT(IN) :: transA, transB
    COMPLEX(KIND=dp), INTENT(IN) :: alpha, beta

    INTEGER :: nrows_A, ncols_A, ncols_B

    nrows_A = SIZE(A%p, 1)
    ncols_A = SIZE(A%p, 2)
    ncols_B = SIZE(B%p, 2)

    CPASSERT(SIZE(B%p,1) .EQ. ncols_A)
    IF (ASSOCIATED(C%p)) THEN
       CPASSERT(SIZE(C%p,1) .EQ. nrows_A)
       CPASSERT(SIZE(C%p,2) .EQ. ncols_B)
    ELSE
       CALL mat_create(C, nrows_A, ncols_B)
    END IF
    CALL ZGEMM(transA, transB, &
               nrows_A, ncols_B, ncols_A, &
               alpha, A%p, nrows_A, B%p, ncols_A, &
               beta, C%p, nrows_A)
  END SUBROUTINE mat_mult_z

  SUBROUTINE mat_inv_cholesky_d(mat, inv)
    ! uses DPOTRF and DPOTRI for inversion using Cholesky
    ! decomposition
    TYPE(mat_d_obj), INTENT(IN) :: mat
    TYPE(mat_d_obj), INTENT(INOUT) :: inv

    INTEGER :: ii, jj, nrows, info
    CHARACTER(LEN=6) :: info_string

    nrows = SIZE(mat%p, 1)

    IF (ASSOCIATED(inv%p)) THEN
       CPASSERT(SIZE(inv%p,1) .EQ. nrows)
       CPASSERT(SIZE(inv%p,2) .EQ. nrows)
    ELSE
       CALL mat_create(inv, nrows, nrows)
    END IF

    ! copy matrix
    CALL DLACPY('U', nrows, nrows, mat%p, nrows, inv%p, nrows)
    ! do cholesky decomposition
    CALL DPOTRF('U', nrows, inv%p, nrows, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, "(I6)") info
       CPABORT("DPOTRF failed with info = "//info_string)
    END IF
    ! do inversion
    CALL DPOTRI('U', nrows, inv%p, nrows, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, "(I6)") info
       CPABORT("DPOTRI failed with info = "//info_string)
    END IF
    ! copy the upper triangular part to lower to get the full matrix
    DO ii = 1, nrows
       DO jj = 1, ii-1
          inv%p(ii,jj) = inv%p(jj,ii)
       END DO
    END DO
  END SUBROUTINE mat_inv_cholesky_d

  SUBROUTINE mat_inv_cholesky_z(mat, inv)
    ! uses ZPOTRF and ZPOTRI for inversion using Cholesky
    ! decomposition
    TYPE(mat_z_obj), INTENT(IN) :: mat
    TYPE(mat_z_obj), INTENT(INOUT) :: inv

    INTEGER :: ii, jj, nrows, info
    CHARACTER(LEN=6) :: info_string

    nrows = SIZE(mat%p, 1)

    IF (ASSOCIATED(inv%p)) THEN
       CPASSERT(SIZE(inv%p,1) .EQ. nrows)
       CPASSERT(SIZE(inv%p,2) .EQ. nrows)
    ELSE
       CALL mat_create(inv, nrows, nrows)
    END IF

    ! copy matrix
    CALL ZLACPY('U', nrows, nrows, mat%p, nrows, inv%p, nrows)
    ! do cholesky decomposition
    CALL ZPOTRF('U', nrows, inv%p, nrows, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, "(I6)") info
       CPABORT("ZPOTRF failed with info = "//info_string)
    END IF
    ! do inversion
    CALL ZPOTRI('U', nrows, inv%p, nrows, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, "(I6)") info
       CPABORT("ZPOTRI failed with info = "//info_string)
    END IF
    ! copy the upper triangular part to lower to get the full matrix
    DO ii = 1, nrows
       DO jj = 1, ii-1
          inv%p(ii,jj) = inv%p(jj,ii)
       END DO
    END DO
  END SUBROUTINE mat_inv_cholesky_z

  SUBROUTINE mat_read_d(mat, filename)
    TYPE(mat_d_obj), INTENT(INOUT) :: mat
    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER, PARAMETER :: UNIT_NR = 100
    INTEGER :: nrows, ncols, ii
    ! open file
    OPEN(UNIT_NR, file=filename)
    ! read into matrix
    READ (UNIT_NR, FMT=*) nrows, ncols
    IF (ASSOCIATED(mat%p)) DEALLOCATE(mat%p)
    ALLOCATE(mat%p(nrows,ncols))
    mat%p = 0.0_dp
    DO ii = 1, nrows
       READ (UNIT_NR, FMT=*) mat%p(ii,:)
    END DO
    ! close file
    CLOSE(UNIT_NR)
  END SUBROUTINE mat_read_d

  SUBROUTINE mat_read_z(mat, filename)
    TYPE(mat_z_obj), INTENT(INOUT) :: mat
    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER, PARAMETER :: UNIT_NR = 100
    INTEGER :: nrows, ncols, ii, jj
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: tmp
    ! open file
    OPEN(UNIT_NR, file=filename)
    ! read into matrix
    READ (UNIT_NR, FMT=*) nrows, ncols
    IF (ASSOCIATED(mat%p)) DEALLOCATE(mat%p)
    ALLOCATE(mat%p(nrows,ncols))
    mat%p = (0.0_dp,0.0_dp)
    ALLOCATE(tmp(nrows,2*ncols))
    tmp = 0.0_dp
    DO ii = 1, nrows
       READ (UNIT_NR, FMT=*) tmp(ii,:)
    END DO
    DO jj = 1, ncols
       DO ii = 1, nrows
          mat%p(ii,jj) = CMPLX(tmp(ii,2*jj-1), tmp(ii,2*jj), KIND=dp)
       END DO
    END DO
    DEALLOCATE(tmp)
    ! close file
    CLOSE(UNIT_NR)
  END SUBROUTINE mat_read_z

  SUBROUTINE mat_write_d(mat, filename)
    TYPE(mat_d_obj), INTENT(IN) :: mat
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename

    INTEGER :: my_unit_nr, nrows, ncols, ii, jj

    CPASSERT(ASSOCIATED(mat%p))
    nrows = SIZE(mat%p, 1)
    ncols = SIZE(mat%p, 2)
    IF (PRESENT(filename)) THEN
       my_unit_nr = 100
       CPASSERT(my_unit_nr .NE. default_output_unit)
       OPEN(my_unit_nr, file=filename)
    ELSE
       my_unit_nr = default_output_unit
    END IF
    ! write matrix
    WRITE (my_unit_nr, FMT="(2I6)") nrows, ncols
    DO ii = 1, nrows
       DO jj = 1, ncols
          WRITE (my_unit_nr, FMT="(F10.7,2X)", ADVANCE="no") mat%p(ii,jj)
       END DO
       WRITE (my_unit_nr, *) ""
    END DO
    ! close file if we need to
    IF (PRESENT(filename)) THEN
       CLOSE(my_unit_nr)
    END IF
  END SUBROUTINE mat_write_d

  SUBROUTINE mat_write_z(mat, filename)
    TYPE(mat_z_obj), INTENT(IN) :: mat
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename

    INTEGER :: my_unit_nr, nrows, ncols, ii, jj

    CPASSERT(ASSOCIATED(mat%p))
    nrows = SIZE(mat%p, 1)
    ncols = SIZE(mat%p, 2)
    IF (PRESENT(filename)) THEN
       my_unit_nr = 100
       CPASSERT(my_unit_nr .NE. default_output_unit)
       OPEN(my_unit_nr, file=filename)
    ELSE
       my_unit_nr = default_output_unit
    END IF
    ! write matrix
    WRITE (my_unit_nr, FMT="(2I6)") nrows, ncols
    DO ii = 1, nrows
       DO jj = 1, ncols
          WRITE (my_unit_nr, FMT="(F10.7,1X,F10.7,3X)", ADVANCE="no") mat%p(ii,jj)
       END DO
       WRITE (my_unit_nr, *) ""
    END DO
    ! close file if we need to
    IF (PRESENT(filename)) THEN
       CLOSE(my_unit_nr)
    END IF
  END SUBROUTINE mat_write_z

  PURE FUNCTION mat_nrows_d(mat) RESULT(res)
    TYPE(mat_d_obj), INTENT(IN) :: mat
    INTEGER :: res
    res = SIZE(mat%p,1)
  END FUNCTION mat_nrows_d

  PURE FUNCTION mat_ncols_d(mat) RESULT(res)
    TYPE(mat_d_obj), INTENT(IN) :: mat
    INTEGER :: res
    res = SIZE(mat%p,2)
  END FUNCTION mat_ncols_d

  PURE FUNCTION mat_nrows_z(mat) RESULT(res)
    TYPE(mat_z_obj), INTENT(IN) :: mat
    INTEGER :: res
    res = SIZE(mat%p,1)
  END FUNCTION mat_nrows_z

  PURE FUNCTION mat_ncols_z(mat) RESULT(res)
    TYPE(mat_z_obj), INTENT(IN) :: mat
    INTEGER :: res
    res = SIZE(mat%p,2)
  END FUNCTION mat_ncols_z


END MODULE matrix_types
