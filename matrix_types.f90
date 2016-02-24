MODULE matrix_types

  USE kinds, ONLY: dp
  USE machine, ONLY: default_output_unit

#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! public types
  PUBLIC :: mat_d_obj, &
            mat_z_obj

  ! public parameters
  PUBLIC :: MAT_GENERAL, &
            MAT_SYMMETRIC, &
            MAT_HERMITIAN

  ! public methods
  PUBLIC :: mat_create,&
            mat_release, &
            mat_symmetry, &
            mat_real_to_complex, &
            mat_mult, &
            mat_inv_cholesky, &
            mat_inv_lu, &
            mat_read, &
            mat_write, &
            mat_nrows, &
            mat_ncols, &
            mat_axpy, &
            mat_zero, &
            mat_copy

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'matrix_types'
  INTEGER, PRIVATE, SAVE :: last_mat_id = 0

  ! symmetry types
  INTEGER, PARAMETER :: MAT_GENERAL = 0, &
                        MAT_SYMMETRIC = 1, &
                        MAT_HERMITIAN = 2

  ! lapack function declarations
  REAL(KIND=dp), EXTERNAL :: DLANGE

  TYPE mat_d_data
     REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: p
     INTEGER :: symmetry
     INTEGER :: id_nr, ref_count
  END type mat_d_data

  TYPE mat_d_obj
     TYPE(mat_d_data), POINTER, PRIVATE :: obj => NULL()
  END TYPE mat_d_obj

  TYPE mat_z_data
     COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: p
     INTEGER :: symmetry
     INTEGER :: id_nr, ref_count
  END type mat_z_data

  TYPE mat_z_obj
     TYPE(mat_z_data), POINTER, PRIVATE :: obj => NULL()
  END TYPE mat_z_obj

  INTERFACE mat_create
     MODULE PROCEDURE mat_create_d
     MODULE PROCEDURE mat_create_z
  END INTERFACE mat_create

  INTERFACE mat_release
     MODULE PROCEDURE mat_release_d
     MODULE PROCEDURE mat_release_z
  END INTERFACE mat_release

  INTERFACE mat_symmetry
     MODULE PROCEDURE mat_symmetry_d
     MODULE PROCEDURE mat_symmetry_z
  END INTERFACE mat_symmetry

  INTERFACE mat_mult
     MODULE PROCEDURE mat_mult_d
     MODULE PROCEDURE mat_mult_z
  END INTERFACE mat_mult

  INTERFACE mat_inv_cholesky
     MODULE PROCEDURE mat_inv_cholesky_d
     MODULE PROCEDURE mat_inv_cholesky_z
  END INTERFACE mat_inv_cholesky

  INTERFACE mat_inv_lu
     MODULE PROCEDURE mat_inv_lu_d
     MODULE PROCEDURE mat_inv_lu_z
  END INTERFACE mat_inv_lu

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

  INTERFACE mat_axpy
     MODULE PROCEDURE mat_axpy_d
     MODULE PROCEDURE mat_axpy_z
  END INTERFACE mat_axpy

  INTERFACE mat_zero
     MODULE PROCEDURE mat_zero_d
     MODULE PROCEDURE mat_zero_z
  END INTERFACE mat_zero

  INTERFACE mat_copy
     MODULE PROCEDURE mat_copy_d
     MODULE PROCEDURE mat_copy_z
  END INTERFACE mat_copy

  INTERFACE mat_norm_d
     MODULE PROCEDURE mat_norm_d
     MODULE PROCEDURE mat_norm_z
  END INTERFACE mat_norm_d

CONTAINS

  SUBROUTINE mat_retain_d(mat)
    TYPE(mat_d_obj), INTENT(INOUT) :: mat
    CPASSERT(ASSOCIATED(mat%obj))
    mat%obj%ref_count = mat%obj%ref_count + 1
  END SUBROUTINE mat_retain_d

  SUBROUTINE mat_retain_z(mat)
    TYPE(mat_z_obj), INTENT(INOUT) :: mat
    CPASSERT(ASSOCIATED(mat%obj))
    mat%obj%ref_count = mat%obj%ref_count + 1
  END SUBROUTINE mat_retain_z

  SUBROUTINE mat_create_d(mat, nrows, ncols, symmetry)
    TYPE(mat_d_obj), INTENT(INOUT) :: mat
    INTEGER, OPTIONAL :: symmetry
    INTEGER, INTENT(IN) :: nrows, ncols
    CPASSERT(.NOT. ASSOCIATED(mat%obj))
    ALLOCATE(mat%obj)
    ALLOCATE(mat%obj%p(nrows,ncols))
    mat%obj%p = 0.0_dp
    IF (PRESENT(symmetry)) THEN
       mat%obj%symmetry = symmetry
    ELSE
       mat%obj%symmetry = MAT_GENERAL
    END IF
    mat%obj%ref_count = 1
    ! book keeping
    mat%obj%id_nr = last_mat_id + 1
    last_mat_id = mat%obj%id_nr
  END SUBROUTINE mat_create_d

  SUBROUTINE mat_create_z(mat, nrows, ncols, symmetry)
    TYPE(mat_z_obj), INTENT(INOUT) :: mat
    INTEGER, OPTIONAL :: symmetry
    INTEGER, INTENT(IN) :: nrows, ncols
    CPASSERT(.NOT. ASSOCIATED(mat%obj))
    ALLOCATE(mat%obj)
    ALLOCATE(mat%obj%p(nrows,ncols))
    mat%obj%p = (0.0_dp,0.0_dp)
    IF (PRESENT(symmetry)) THEN
       mat%obj%symmetry = symmetry
    ELSE
       mat%obj%symmetry = MAT_GENERAL
    END IF
    mat%obj%ref_count = 1
    ! book keeping
    mat%obj%id_nr = last_mat_id + 1
    last_mat_id = mat%obj%id_nr
  END SUBROUTINE mat_create_z

  SUBROUTINE mat_release_d(mat)
    TYPE(mat_d_obj), INTENT(INOUT) :: mat
    IF (ASSOCIATED(mat%obj)) THEN
       mat%obj%ref_count = mat%obj%ref_count - 1
       IF (mat%obj%ref_count .EQ. 0) THEN
          IF (ALLOCATED(mat%obj%p)) DEALLOCATE(mat%obj%p)
          DEALLOCATE(mat%obj)
       END IF
       NULLIFY(mat%obj)
    END IF
  END SUBROUTINE mat_release_d

  SUBROUTINE mat_release_z(mat)
    TYPE(mat_z_obj), INTENT(INOUT) :: mat
    IF (ASSOCIATED(mat%obj)) THEN
       mat%obj%ref_count = mat%obj%ref_count - 1
       IF (mat%obj%ref_count .EQ. 0) THEN
          IF (ALLOCATED(mat%obj%p)) DEALLOCATE(mat%obj%p)
          DEALLOCATE(mat%obj)
       END IF
       NULLIFY(mat%obj)
    END IF
  END SUBROUTINE mat_release_z

  PURE FUNCTION mat_symmetry_d(mat) RESULT(res)
    TYPE(mat_d_obj), INTENT(IN) :: mat
    INTEGER :: res
    IF (ASSOCIATED(mat%obj)) THEN
       res = mat%obj%symmetry
    ELSE
       res = MAT_GENERAL
    END IF
  END FUNCTION mat_symmetry_d

  PURE FUNCTION mat_symmetry_z(mat) RESULT(res)
    TYPE(mat_z_obj), INTENT(IN) :: mat
    INTEGER :: res
    IF (ASSOCIATED(mat%obj)) THEN
       res = mat%obj%symmetry
    ELSE
       res = MAT_GENERAL
    END IF
  END FUNCTION mat_symmetry_z

  SUBROUTINE mat_real_to_complex(mat_d, mat_z)
    TYPE(mat_d_obj), INTENT(IN) :: mat_d
    TYPE(mat_z_obj), INTENT(INOUT) :: mat_z
    CPASSERT(ASSOCIATED(mat_d%obj))
    IF (ASSOCIATED(mat_z%obj)) CALL mat_release(mat_z)
    CALL mat_create(mat_z, mat_nrows(mat_d), mat_ncols(mat_d))
    mat_z%obj%p(:,:) = CMPLX(mat_d%obj%p(:,:), 0.0, KIND=dp)
    mat_z%obj%symmetry = mat_d%obj%symmetry
  END SUBROUTINE mat_real_to_complex

  ! SUBROUTINE mat_mult_d(transA, transB, alpha, A, B, beta, C)
  !   TYPE(mat_d_obj), INTENT(IN) :: A, B
  !   TYPE(mat_d_obj), INTENT(INOUT) :: C
  !   CHARACTER, INTENT(IN) :: transA, transB
  !   REAL(KIND=dp), INTENT(IN) :: alpha, beta

  !   INTEGER :: nrows_A, ncols_A, ncols_B
  !   LOGICAL :: do_symm
  !   INTEGER :: ii, jj

  !   nrows_A = mat_nrows(A)
  !   ncols_A = mat_ncols(A)
  !   ncols_B = mat_ncols(B)

  !   CPASSERT(mat_nrows(B) .EQ. ncols_A)
  !   IF (ASSOCIATED(C%obj)) THEN
  !      IF (ALLOCATED(C%obj%p)) CPASSERT(mat_nrows(C) .EQ. nrows_A)
  !      IF (ALLOCATED(C%obj%p)) CPASSERT(mat_ncols(C) .EQ. ncols_B)
  !   ELSE
  !      CALL mat_create(C, nrows_A, ncols_B, MAT_GENERAL)
  !   END IF

  !   ! if either A or B are symmetric, then their corresponding trans
  !   ! must not be 'T'
  !   IF (mat_symmetr(A) .EQ. MAT_SYMMETRIC) CPASSERT(transA .EQ. 'N')
  !   IF (mat_symmetr(B) .EQ. MAT_SYMMETRIC) CPASSERT(transB .EQ. 'N')

  !   do_symm = .FALSE.
  !   IF ((mat_symmetry(A) .EQ. MAT_SYMMETRIC) .AND. &
  !       (transB .EQ. 'N')) THEN
  !      do_symm = .TRUE.
  !      ! need to make B full if it is symmetric
  !      IF (mat_symmetry(B) .EQ. MAT_SYMMETRIC) THEN
  !         DO ii = 1, ncols_A
  !            DO jj = 1, ii
  !               B%obj%p(ii,jj) = B%obj%p(jj,ii)
  !            END DO
  !         END DO
  !      END IF
  !   ELSE IF ((mat_symmetry(B) .EQ. MAT_SYMMETRIC) .AND. &
  !            (transA .EQ. 'N')) THEN
  !      do_symm = .TRUE.
  !      ! need to make A full if it is symmetric
  !      IF (mat_symmetry(A) .EQ. MAT_SYMMETRIC) THEN
  !         DO ii = 1, nrows_A
  !            DO jj = 1, ii
  !               A%obj%p(ii,jj) = A%obj%p(jj,ii)
  !            END DO
  !         END DO
  !      END IF
  !   END IF
  !   IF (do_symm) THEN
  !      IF (mat_symmetry(A) .EQ. MAT_SYMMETRIC) THEN
  !         CALL DSYMM('L', 'U', nrows_A, ncols_B, &
  !                    alpha, A%obj%p, nrows_A, B%obj%p, ncols_A, &
  !                    beta, C%obj%p, nrows_A)
  !      ELSE IF (mat_symmetry(B) .EQ. MAT_SYMMETRIC) THEN
  !         CALL DSYMM('R', 'U', nrows_A, ncols_B, &
  !                    alpha, B%obj%p, ncols_A, A%obj%p, nrows_A, &
  !                    beta, C%obj%p, nrows_A)
  !      END IF
  !   ELSE
  !      CALL DGEMM(transA, transB, &
  !                 nrows_A, ncols_B, ncols_A, &
  !                 alpha, A%obj%p, nrows_A, B%obj%p, ncols_A, &
  !                 beta, C%p, nrows_A)
  !   END IF
  ! END SUBROUTINE mat_mult_d

  SUBROUTINE mat_mult_d(transA, transB, alpha, A, B, beta, C)
    TYPE(mat_d_obj), INTENT(IN) :: A, B
    TYPE(mat_d_obj), INTENT(INOUT) :: C
    CHARACTER, INTENT(IN) :: transA, transB
    REAL(KIND=dp), INTENT(IN) :: alpha, beta

    INTEGER :: nrows_A, ncols_A, ncols_B

    nrows_A = mat_nrows(A)
    ncols_A = mat_ncols(A)
    ncols_B = mat_ncols(B)

    CPASSERT(mat_nrows(B) .EQ. ncols_A)
    IF (ASSOCIATED(C%obj)) THEN
       CPASSERT(mat_nrows(C) .EQ. nrows_A)
       CPASSERT(mat_ncols(C) .EQ. ncols_B)
    ELSE
       CALL mat_create(C, nrows_A, ncols_B)
    END IF
    CALL DGEMM(transA, transB, &
               nrows_A, ncols_B, ncols_A, &
               alpha, A%obj%p, nrows_A, B%obj%p, ncols_A, &
               beta, C%obj%p, nrows_A)
  END SUBROUTINE mat_mult_d

  SUBROUTINE mat_mult_z(transA, transB, alpha, A, B, beta, C)
    TYPE(mat_z_obj), INTENT(IN) :: A, B
    TYPE(mat_z_obj), INTENT(INOUT) :: C
    CHARACTER, INTENT(IN) :: transA, transB
    COMPLEX(KIND=dp), INTENT(IN) :: alpha, beta

    INTEGER :: nrows_A, ncols_A, ncols_B

    nrows_A = mat_nrows(A)
    ncols_A = mat_ncols(A)
    ncols_B = mat_ncols(B)

    CPASSERT(mat_nrows(B) .EQ. ncols_A)
    IF (ASSOCIATED(C%obj)) THEN
       CPASSERT(mat_nrows(C) .EQ. nrows_A)
       CPASSERT(mat_ncols(C) .EQ. ncols_B)
    ELSE
       CALL mat_create(C, nrows_A, ncols_B)
    END IF
    CALL ZGEMM(transA, transB, &
               nrows_A, ncols_B, ncols_A, &
               alpha, A%obj%p, nrows_A, B%obj%p, ncols_A, &
               beta, C%obj%p, nrows_A)
  END SUBROUTINE mat_mult_z


  SUBROUTINE mat_inv_cholesky_d(mat, inv)
    ! uses DPOTRF and DPOTRI for inversion using Cholesky
    ! decomposition
    TYPE(mat_d_obj), INTENT(IN) :: mat
    TYPE(mat_d_obj), INTENT(INOUT) :: inv

    INTEGER :: ii, jj, nrows, info
    CHARACTER(LEN=6) :: info_string

    nrows = mat_nrows(mat)
    CPASSERT(mat_ncols(mat) .EQ. nrows)

    IF (ASSOCIATED(inv%obj)) THEN
       CPASSERT(mat_nrows(inv) .EQ. nrows)
       CPASSERT(mat_ncols(inv) .EQ. nrows)
    ELSE
       CALL mat_create(inv, nrows, nrows)
    END IF

    ! copy matrix
    CALL DLACPY('U', nrows, nrows, mat%obj%p, nrows, inv%obj%p, nrows)
    ! do cholesky decomposition
    CALL DPOTRF('U', nrows, inv%obj%p, nrows, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, "(I6)") info
       CPABORT("DPOTRF failed with info = "//info_string)
    END IF
    ! do inversion
    CALL DPOTRI('U', nrows, inv%obj%p, nrows, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, "(I6)") info
       CPABORT("DPOTRI failed with info = "//info_string)
    END IF
    ! copy the upper triangular part to lower to get the full matrix
    DO ii = 1, nrows
       DO jj = 1, ii-1
          inv%obj%p(ii,jj) = inv%obj%p(jj,ii)
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

    nrows = mat_nrows(mat)
    CPASSERT(mat_ncols(mat) .EQ. nrows)

    IF (ASSOCIATED(inv%obj)) THEN
       CPASSERT(mat_nrows(mat) .EQ. nrows)
       CPASSERT(mat_ncols(mat) .EQ. nrows)
    ELSE
       CALL mat_create(inv, nrows, nrows)
    END IF

    ! copy matrix
    CALL ZLACPY('U', nrows, nrows, mat%obj%p, nrows, inv%obj%p, nrows)
    ! do cholesky decomposition
    CALL ZPOTRF('U', nrows, inv%obj%p, nrows, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, "(I6)") info
       CPABORT("ZPOTRF failed with info = "//info_string)
    END IF
    ! do inversion
    CALL ZPOTRI('U', nrows, inv%obj%p, nrows, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, "(I6)") info
       CPABORT("ZPOTRI failed with info = "//info_string)
    END IF
    ! copy the upper triangular part to lower to get the full matrix
    DO ii = 1, nrows
       DO jj = 1, ii-1
          inv%obj%p(ii,jj) = inv%obj%p(jj,ii)
       END DO
    END DO
  END SUBROUTINE mat_inv_cholesky_z

  SUBROUTINE mat_inv_lu_d(mat, inv)
    TYPE(mat_d_obj), INTENT(IN) :: mat
    TYPE(mat_d_obj), INTENT(INOUT) :: inv

    INTEGER :: nrows, info, lwork
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
    CHARACTER(LEN=6) :: info_string
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: work
    nrows = mat_nrows(mat)
    IF (ASSOCIATED(inv%obj)) THEN
       CPASSERT(mat_nrows(inv) .EQ. nrows)
       CPASSERT(mat_ncols(inv) .EQ. nrows)
    ELSE
       CALL mat_create(inv, nrows, nrows)
    END IF
    ! copy matrix
    CALL DLACPY('N', nrows, nrows, mat%obj%p, nrows, inv%obj%p, nrows)
    ! do LU factorisation
    ALLOCATE(ipiv(nrows))
    CALL DGETRF(nrows, nrows, inv%obj%p, nrows, ipiv, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, FMT="(I6)") info
       CPABORT("DGETRF failed with info = "//info_string)
    END IF
    ! do inversion
    ALLOCATE(work(1))
    CALL DGETRI(nrows, inv%obj%p, nrows, ipiv, work, -1, info)
    lwork = work(1)
    DEALLOCATE(work)
    ALLOCATE(work(lwork))
    CALL DGETRI(nrows, inv%obj%p, nrows, ipiv, work, lwork, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, FMT="(I6)") info
       CPABORT("DGETRI failed with info = "//info_string)
    END IF
    ! cleanup
    DEALLOCATE(ipiv)
    DEALLOCATE(work)
  END SUBROUTINE mat_inv_lu_d

  SUBROUTINE mat_inv_lu_z(mat, inv)
    TYPE(mat_z_obj), INTENT(IN) :: mat
    TYPE(mat_z_obj), INTENT(INOUT) :: inv

    INTEGER :: nrows, info, lwork
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
    CHARACTER(LEN=6) :: info_string
    COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE :: work
    nrows = mat_nrows(mat)
    IF (ASSOCIATED(inv%obj)) THEN
       CPASSERT(mat_nrows(inv) .EQ. nrows)
       CPASSERT(mat_ncols(inv) .EQ. nrows)
    ELSE
       CALL mat_create(inv, nrows, nrows)
    END IF
    ! copy matrix
    CALL ZLACPY('N', nrows, nrows, mat%obj%p, nrows, inv%obj%p, nrows)
    ! do LU factorisation
    ALLOCATE(ipiv(nrows))
    CALL ZGETRF(nrows, nrows, inv%obj%p, nrows, ipiv, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, FMT="(I6)") info
       CPABORT("ZGETRF failed with info = "//info_string)
    END IF
    ! do inversion
    ALLOCATE(work(1))
    CALL ZGETRI(nrows, inv%obj%p, nrows, ipiv, work, -1, info)
    lwork = work(1)
    DEALLOCATE(work)
    ALLOCATE(work(lwork))
    CALL ZGETRI(nrows, inv%obj%p, nrows, ipiv, work, lwork, info)
    IF (info .NE. 0) THEN
       WRITE (info_string, FMT="(I6)") info
       CPABORT("ZGETRI failed with info = "//info_string)
    END IF
    ! cleanup
    DEALLOCATE(ipiv)
    DEALLOCATE(work)
  END SUBROUTINE mat_inv_lu_z

  SUBROUTINE mat_read_d(mat, filename)
    TYPE(mat_d_obj), INTENT(INOUT) :: mat
    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER, PARAMETER :: UNIT_NR = 100
    INTEGER :: nrows, ncols, ii
    ! open file
    OPEN(UNIT_NR, file=filename)
    ! read into matrix
    READ (UNIT_NR, FMT=*) nrows, ncols
    IF (ASSOCIATED(mat%obj)) CALL mat_release(mat)
    CALL mat_create(mat, nrows, ncols)
    DO ii = 1, nrows
       READ (UNIT_NR, FMT=*) mat%obj%p(ii,:)
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
    IF (ASSOCIATED(mat%obj)) CALL mat_release(mat)
    CALL mat_create(mat, nrows, ncols)
    ALLOCATE(tmp(nrows,2*ncols))
    tmp = 0.0_dp
    DO ii = 1, nrows
       READ (UNIT_NR, FMT=*) tmp(ii,:)
    END DO
    DO jj = 1, ncols
       DO ii = 1, nrows
          mat%obj%p(ii,jj) = CMPLX(tmp(ii,2*jj-1), tmp(ii,2*jj), KIND=dp)
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

    CPASSERT(ASSOCIATED(mat%obj))
    nrows = mat_nrows(mat)
    ncols = mat_ncols(mat)
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
          WRITE (my_unit_nr, FMT="(F10.7,2X)", ADVANCE="no") mat%obj%p(ii,jj)
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

    CPASSERT(ASSOCIATED(mat%obj))
    nrows = mat_nrows(mat)
    ncols = mat_ncols(mat)
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
          WRITE (my_unit_nr, FMT="(F10.7,1X,F10.7,3X)", ADVANCE="no") mat%obj%p(ii,jj)
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
    res = SIZE(mat%obj%p,1)
  END FUNCTION mat_nrows_d

  PURE FUNCTION mat_nrows_z(mat) RESULT(res)
    TYPE(mat_z_obj), INTENT(IN) :: mat
    INTEGER :: res
    res = SIZE(mat%obj%p,1)
  END FUNCTION mat_nrows_z

  PURE FUNCTION mat_ncols_d(mat) RESULT(res)
    TYPE(mat_d_obj), INTENT(IN) :: mat
    INTEGER :: res
    res = SIZE(mat%obj%p,2)
  END FUNCTION mat_ncols_d

  PURE FUNCTION mat_ncols_z(mat) RESULT(res)
    TYPE(mat_z_obj), INTENT(IN) :: mat
    INTEGER :: res
    res = SIZE(mat%obj%p,2)
  END FUNCTION mat_ncols_z

  SUBROUTINE mat_axpy_d(a, transX, X, Y)
    ! computes Y = aX^(T) + Y
    REAL(KIND=dp), INTENT(IN) :: a
    CHARACTER(LEN=*), INTENT(IN) :: transX
    TYPE(mat_d_obj), INTENT(IN) :: X
    TYPE(mat_d_obj), INTENT(INOUT) :: Y

    INTEGER :: ii, jj, nrows, ncols
    nrows = mat_nrows(Y)
    ncols = mat_ncols(Y)
    IF (transX .EQ. 'T') THEN
       CPASSERT(mat_nrows(X) .EQ. ncols)
       CPASSERT(mat_ncols(X) .EQ. nrows)
       DO jj = 1, ncols
          DO ii = 1, nrows
             Y%obj%p(ii,jj) = a * X%obj%p(jj,ii) + Y%obj%p(ii,jj)
          END DO
       END DO
    ELSE
       CPASSERT(mat_nrows(X) .EQ. nrows)
       CPASSERT(mat_ncols(X) .EQ. ncols)
       Y%obj%p = a*X%obj%p + Y%obj%p
    END IF
  END SUBROUTINE mat_axpy_d

  SUBROUTINE mat_axpy_z(a, transX, X, Y)
    ! computes Y = aX^(T) + Y
    COMPLEX(KIND=dp), INTENT(IN) :: a
    CHARACTER(LEN=*), INTENT(IN) :: transX
    TYPE(mat_z_obj), INTENT(IN) :: X
    TYPE(mat_z_obj), INTENT(INOUT) :: Y

    INTEGER :: ii, jj, nrows, ncols
    nrows = mat_nrows(Y)
    ncols = mat_ncols(Y)
    SELECT CASE (transX)
       CASE ('T')
          CPASSERT(mat_nrows(X) .EQ. ncols)
          CPASSERT(mat_ncols(X) .EQ. nrows)
          DO jj = 1, ncols
             DO ii = 1, nrows
                Y%obj%p(ii,jj) = a * X%obj%p(jj,ii) + Y%obj%p(ii,jj)
             END DO
          END DO
       CASE ('H')
          CPASSERT(mat_nrows(X) .EQ. ncols)
          CPASSERT(mat_ncols(X) .EQ. nrows)
          DO jj = 1, ncols
             DO ii = 1, nrows
                Y%obj%p(ii,jj) = a * CONJG(X%obj%p(jj,ii)) + Y%obj%p(ii,jj)
             END DO
          END DO
       CASE DEFAULT
          CPASSERT(mat_nrows(X) .EQ. nrows)
          CPASSERT(mat_ncols(X) .EQ. ncols)
          Y%obj%p = a*X%obj%p + Y%obj%p
       END SELECT
  END SUBROUTINE mat_axpy_z

  SUBROUTINE mat_zero_d(mat)
    TYPE(mat_d_obj), INTENT(INOUT) :: mat
    ! assume assoicated obj always leads to allocated obj%p, as no
    ! module method should allocate obj but not also allocate obj%p
    IF (ASSOCIATED(mat%obj)) THEN
       mat%obj%p = 0.0_dp
    END IF
  END SUBROUTINE mat_zero_d

  SUBROUTINE mat_zero_z(mat)
    TYPE(mat_z_obj), INTENT(INOUT) :: mat
    ! assume assoicated obj always leads to allocated obj%p, as no
    ! module method should allocate obj but not also allocate obj%p
    IF (ASSOCIATED(mat%obj)) THEN
       mat%obj%p = (0.0_dp,0.0_dp)
    END IF
  END SUBROUTINE mat_zero_z

  SUBROUTINE mat_copy_d(A, B)
    TYPE(mat_d_obj), INTENT(IN) :: A
    TYPE(mat_d_obj), INTENT(INOUT) :: B
    INTEGER :: nrows, ncols
    nrows = mat_nrows(A)
    ncols = mat_ncols(A)
    IF (ASSOCIATED(B%obj)) CALL mat_release(B)
    CALL mat_create(B, mat_nrows(A), mat_ncols(A))
    CALL DLACPY('N', nrows, ncols, A%obj%p, nrows, B%obj%p, nrows)
  END SUBROUTINE mat_copy_d

  SUBROUTINE mat_copy_z(A, B)
    TYPE(mat_z_obj), INTENT(IN) :: A
    TYPE(mat_z_obj), INTENT(INOUT) :: B
    INTEGER :: nrows, ncols
    nrows = mat_nrows(A)
    ncols = mat_ncols(A)
    IF (ASSOCIATED(B%obj)) CALL mat_release(B)
    CALL mat_create(B, mat_nrows(A), mat_ncols(A))
    CALL ZLACPY('N', nrows, ncols, A%obj%p, nrows, B%obj%p, nrows)
  END SUBROUTINE mat_copy_z

  FUNCTION mat_norm_d(mat) RESULT(res)
    TYPE(mat_d_obj), INTENT(IN) :: mat
    REAL(KIND=dp) :: res
    INTEGER :: nrows, ncols
    REAL(KIND=dp), DIMENSION(1) :: work
    nrows = mat_nrows(mat)
    ncols = mat_ncols(mat)
    ! work is a dummy for Frobenius norm
    res = DLANGE('F', nrows, ncols, mat%obj%p, nrows, work)
  END FUNCTION mat_norm_d

  FUNCTION mat_norm_z(mat) RESULT(res)
    TYPE(mat_z_obj), INTENT(IN) :: mat
    COMPLEX(KIND=dp) :: res
    INTEGER :: nrows, ncols
    COMPLEX(KIND=dp), DIMENSION(1) :: work
    nrows = mat_nrows(mat)
    ncols = mat_ncols(mat)
    ! work is a dummy for Frobenius norm
    res = DLANGE('F', nrows, ncols, mat%obj%p, nrows, work)
  END FUNCTION mat_norm_z


END MODULE matrix_types
