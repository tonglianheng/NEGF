MODULE matrix_types

   USE kinds, ONLY: dp
   USE machine, ONLY: default_output_unit
   USE cp_fm_types, ONLY: cp_fm_type, &
                          cp_fm_create, &
                          cp_fm_to_fm_triangular, &
                          cp_fm_set_submatrix, &
                          cp_fm_release, &
                          cp_fm_get_info, &
                          cp_fm_to_fm_triangular
   USE cp_cfm_types, ONLY: cp_cfm_types, &
                           cp_cfm_create, &
                           cp_cfm_to_cfm_triangular, &
                           cp_cfm_set_submatrix, &
                           cp_cfm_release, &
                           cp_cfm_get_info, &
                           cp_fm_to_cfm
   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_fm_struct, ONLY: cp_fm_struct_type, &
                           cp_fm_struct_create, &
                           cp_fm_struct_get
   USE cp_fm_basic_linalg, ONLY: cp_fm_gemm, &
                                 cp_fm_cholesky_decompose, &
                                 cp_fm_cholesky_invert
   USE cp_cfm_basic_linalg, ONLY: cp_cfm_gemm, &
                                  cp_cfm_cholesky_decompose, &
                                  cp_cfm_cholesky_invert

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
   PUBLIC :: is_same_obj, &
             mat_associate, &
             mat_axpy, &
             mat_blacs_env, &
             mat_copy, &
             mat_create, &
             mat_inv_cholesky, &
             mat_inv_lu, &
             mat_mult, &
             mat_ncols, &
             mat_norm, &
             mat_nrows, &
             mat_read, &
             mat_real_to_complex, &
             mat_release, &
             mat_scale, &
             mat_struct_equiv, &
             mat_symmetry, &
             mat_trace, &
             mat_write, &
             mat_zero

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'matrix_types'
   INTEGER, PRIVATE, SAVE :: last_mat_id = 0

   ! symmetry types
   INTEGER, PARAMETER :: MAT_GENERAL = 0, &
                         MAT_SYMMETRIC = 1, &
                         MAT_HERMITIAN = 2

   ! lapack function declarations
   REAL(KIND=dp), EXTERNAL :: DLANGE, ZLANGE

   TYPE mat_d_data
      TYPE(cp_fm_type), POINTER :: p
      INTEGER :: symmetry
      INTEGER :: id_nr, ref_count
   END type mat_d_data

   TYPE mat_d_obj
      TYPE(mat_d_data), POINTER, PRIVATE :: obj => NULL()
   END TYPE mat_d_obj

   TYPE mat_z_data
      TYPE(cp_cfm_type), POINTER :: p
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

   INTERFACE mat_retain
      MODULE PROCEDURE mat_retain_d
      MODULE PROCEDURE mat_retain_z
   END INTERFACE mat_retain

   INTERFACE mat_release
      MODULE PROCEDURE mat_release_d
      MODULE PROCEDURE mat_release_z
   END INTERFACE mat_release

   INTERFACE mat_associate
      MODULE PROCEDURE mat_associate_d
      MODULE PROCEDURE mat_associate_z
   END INTERFACE mat_associate

   INTERFACE mat_copy
      MODULE PROCEDURE mat_copy_d
      MODULE PROCEDURE mat_copy_z
   END INTERFACE mat_copy

   INTERFACE mat_scale
      MODULE PROCEDURE mat_scale_d
      MODULE PROCEDURE mat_scale_z
   END INTERFACE mat_scale

   INTERFACE mat_symmetry
      MODULE PROCEDURE mat_symmetry_d
      MODULE PROCEDURE mat_symmetry_z
   END INTERFACE mat_symmetry

   INTERFACE mat_blacs_env
      MODULE PROCEDURE mat_blacs_env_d
      MODULE PROCEDURE mat_blacs_env_z
   END INTERFACE mat_blacs_env

   INTERFACE mat_struct_equiv
      MODULE PROCEDURE mat_struct_equiv_d
      MODULE PROCEDURE mat_struct_equiv_z
   END INTERFACE mat_struct_equiv

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

   INTERFACE mat_norm
      MODULE PROCEDURE mat_norm_d
      MODULE PROCEDURE mat_norm_z
   END INTERFACE mat_norm

   INTERFACE is_same_obj
      MODULE PROCEDURE is_same_obj_d
      MODULE PROCEDURE is_same_obj_z
   END INTERFACE is_same_obj

   INTERFACE mat_trace
      MODULE PROCEDURE mat_trace_d
      MODULE PROCEDURE mat_trace_z
   END INTERFACE mat_trace


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

   SUBROUTINE mat_create_d(mat, nrows, ncols, blacs_env, symmetry, content, name)
      TYPE(mat_d_obj), INTENT(INOUT) :: mat
      INTEGER, INTENT(IN) :: nrows, ncols
      TYPE(cp_blacs_env_type), POINTER :: blacs_env
      INTEGER, OPTIONAL :: symmetry
      REAL(KIND=dp), DIMENSION(:,:), INTENT(IN), OPTIONAL :: content
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: name
      TYPE(cp_fm_struct_type), POINTER :: fmstruct
      CPASSERT(.NOT. ASSOCIATED(mat%obj))
      ALLOCATE(mat%obj)
      CALL cp_fm_struct_create(fmstruct=fmstruct, &
                               context=blacs_env, &
                               nrows_global=nrows, &
                               ncols_global=ncols)
      IF (PRESENT(name)) THEN
         CALL cp_fm_create(mat%obj%p, fmstruct, name=name)
      ELSE
         CALL cp_fm_create(mat%obj%p, fmstruct)
      END IF
      IF (PRESENT(symmetry)) THEN
         mat%obj%symmetry = symmetry
      ELSE
         mat%obj%symmetry = MAT_GENERAL
      END IF
      IF (PRESENT(content)) THEN
         CALL cp_fm_set_submatrix(fm=mat%obj%p, &
                                  new_values=content, &
                                  start_row=1, &
                                  start_col=1, &
                                  n_rows=nrows, &
                                  n_cols=ncols)
      END IF
      mat%obj%ref_count = 1
      ! book keeping
      mat%obj%id_nr = last_mat_id + 1
      last_mat_id = mat%obj%id_nr
   END SUBROUTINE mat_create_d

   SUBROUTINE mat_create_z(mat, nrows, ncols, blacs_env, symmetry, content, name)
      TYPE(mat_z_obj), INTENT(INOUT) :: mat
      INTEGER, INTENT(IN) :: nrows, ncols
      TYPE(cp_blacs_env_type), POINTER :: blacs_env
      INTEGER, OPTIONAL :: symmetry
      COMPLEX(KIND=dp), DIMENSION(:,:), INTENT(IN), OPTIONAL :: content
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: name
      TYPE(cp_fm_struct_type), POINTER :: fmstruct
      CPASSERT(.NOT. ASSOCIATED(mat%obj))
      ALLOCATE(mat%obj)
      CALL cp_fm_struct_create(fmstruct=fmstruct, &
                               context=blacs_env, &
                               nrows_global=nrows, &
                               ncols_global=ncols)
      IF (PRESENT(name)) THEN
         CALL cp_cfm_create(mat%obj%p, fmstruct, name=name)
      ELSE
         CALL cp_cfm_create(mat%obj%p, fmstruct)
      END IF
      IF (PRESENT(symmetry)) THEN
         mat%obj%symmetry = symmetry
      ELSE
         mat%obj%symmetry = MAT_GENERAL
      END IF
      IF (PRESENT(content)) THEN
         CALL cp_cfm_set_submatrix(fm=mat%obj%p, &
                                   new_values=content, &
                                   start_row=1, &
                                   start_col=1, &
                                   n_rows=nrows, &
                                   n_cols=ncols)
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
            IF (ASSOCIATED(mat%obj%p)) CALL cp_fm_release(mat%obj%p)
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
            IF (ASSOCIATED(mat%obj%p)) CALL cp_cfm_release(mat%obj%p)
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

   FUNCTION mat_blacs_env_d(mat) RESULT(res)
      TYPE(mat_d_obj), INTENT(IN) :: mat
      TYPE(cp_blacs_env_type), POINTER :: res
      TYPE(cp_fm_struct_type), POINTER :: struct
      CPASSERT(ASSOCIATED(mat%obj))
      CPASSERT(ASSOCIATED(mat%obj%p))
      CALL cp_fm_get_info(matrix=mat%obj%p, &
                          matrix_struct=struct)
      CALL cp_fm_struct_get(fmstruct=struct, &
                            context=res)
   END FUNCTION mat_blacs_env_d

   FUNCTION mat_blacs_env_z(mat) RESULT(res)
      TYPE(mat_z_obj), INTENT(IN) :: mat
      TYPE(cp_blacs_env_type), POINTER :: res
      TYPE(cp_fm_struct_type), POINTER :: struct
      CPASSERT(ASSOCIATED(mat%obj))
      CPASSERT(ASSOCIATED(mat%obj%p))
      CALL cp_cfm_get_info(matrix=mat%obj%p, &
           matrix_struct=struct)
      CALL cp_fm_struct_get(fmstruct=struct, &
           context=res)
   END FUNCTION mat_blacs_env_z

   FUNCTION mat_struct_equiv_d(mat_a, mat_b) RESULT(res)
      TYPE(mat_d_obj), INTENT(IN) :: mat_a, mat_b
      LOGICAL                     :: res
      TYPE(cp_fm_struct_type), POINTER :: struct_a, struct_b
      CPASSERT(ASSOCIATED(mat_a%obj))
      CPASSERT(ASSOCIATED(mat_b%obj))
      CPASSERT(ASSOCIATED(mat_a%obj%p))
      CPASSERT(ASSOCIATED(mat_b%obj%p))
      CALL cp_fm_get_info(matrix=mat_a%obj%p, &
                          matrix_struct=struct_a)
      CALL cp_fm_get_info(matrix=mat_b%obj%p, &
                          matrix_struct=struct_b)
      res = cp_fm_struct_equivalent(struct_a, struct_b)
   END FUNCTION mat_struct_equiv_d

   FUNCTION mat_struct_equiv_z(mat_a, mat_b) RESULT(res)
      TYPE(mat_z_obj), INTENT(IN) :: mat_a, mat_b
      LOGICAL                     :: res
      TYPE(cp_fm_struct_type), POINTER :: struct_a, struct_b
      CPASSERT(ASSOCIATED(mat_a%obj))
      CPASSERT(ASSOCIATED(mat_b%obj))
      CPASSERT(ASSOCIATED(mat_a%obj%p))
      CPASSERT(ASSOCIATED(mat_b%obj%p))
      CALL cp_cfm_get_info(matrix=mat_a%obj%p, &
                           matrix_struct=struct_a)
      CALL cp_cfm_get_info(matrix=mat_b%obj%p, &
                           matrix_struct=struct_b)
      res = cp_fm_struct_equivalent(struct_a, struct_b)
   END FUNCTION mat_struct_equiv_z

   SUBROUTINE mat_real_to_complex(mat_d, mat_z)
      TYPE(mat_d_obj), INTENT(IN) :: mat_d
      TYPE(mat_z_obj), INTENT(INOUT) :: mat_z
      CPASSERT(ASSOCIATED(mat_d%obj))
      CALL mat_release(mat_z)
      CALL mat_create(mat_z, &
                      mat_nrows(mat_d), &
                      mat_ncols(mat_d), &
                      mat_blacs_env(mat_d))
      CALL cp_fm_to_cfm(msourcer=mat_d%obj%p, &
                        mtarget=mat_z%obj%p)
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

   !   CPASSERT(.NOT. is_same_obj(A,C))
   !   CPASSERT(.NOT. is_same_obj(B,C))

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
      TYPE(cp_blacs_env_type), POINTER :: struct

      CPASSERT(.NOT. is_same_obj(A,C))
      CPASSERT(.NOT. is_same_obj(B,C))

      nrows_A = mat_nrows(A)
      ncols_A = mat_ncols(A)
      ncols_B = mat_ncols(B)
      struct => mat_blacs_env(A)

      CPASSERT(mat_nrows(B) .EQ. ncols_A)
      IF (ASSOCIATED(C%obj)) THEN
         CPASSERT(mat_nrows(C) .EQ. nrows_A)
         CPASSERT(mat_ncols(C) .EQ. ncols_B)
      ELSE
         CALL mat_create(C, nrows_A, ncols_B, struct)
      END IF
      CALL cp_fm_gemm(transA, transB, &
                      nrows_A, ncols_B, ncols_A, &
                      alpha, A%obj%p, B%obj%p, &
                      beta, C%obj%p)
   END SUBROUTINE mat_mult_d

   SUBROUTINE mat_mult_z(transA, transB, alpha, A, B, beta, C)
      TYPE(mat_z_obj), INTENT(IN) :: A, B
      TYPE(mat_z_obj), INTENT(INOUT) :: C
      CHARACTER, INTENT(IN) :: transA, transB
      COMPLEX(KIND=dp), INTENT(IN) :: alpha, beta

      INTEGER :: nrows_A, ncols_A, ncols_B
      TYPE(cp_blacs_env_type), POINTER :: struct

      CPASSERT(.NOT. is_same_obj(A,C))
      CPASSERT(.NOT. is_same_obj(B,C))

      nrows_A = mat_nrows(A)
      ncols_A = mat_ncols(A)
      ncols_B = mat_ncols(B)
      struct => mat_blacs_env(A)

      CPASSERT(mat_nrows(B) .EQ. ncols_A)
      IF (ASSOCIATED(C%obj)) THEN
         CPASSERT(mat_nrows(C) .EQ. nrows_A)
         CPASSERT(mat_ncols(C) .EQ. ncols_B)
      ELSE
         CALL mat_create(C, nrows_A, ncols_B, struct)
      END IF
      CALL cp_cfm_gemm(transA, transB, &
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
      TYPE(mat_d_obj) :: tmp

      CPASSERT(.NOT. is_same_obj(mat,inv))
      nrows = mat_nrows(mat)
      CPASSERT(mat_ncols(mat) .EQ. nrows)

      IF (ASSOCIATED(inv%obj)) THEN
         CPASSERT(mat_nrows(inv) .EQ. nrows)
         CPASSERT(mat_ncols(inv) .EQ. nrows)
         CPASSERT(mat_struct_equiv(mat, inv))
      ELSE
         CALL mat_create(inv, nrows, nrows, mat_blacs_env(mat))
      END IF
      ! copy matrix
      CALL cp_fm_to_fm_triangular(mat%obj%p, inv%obj%p, 'U')
      ! do cholesky decomposition
      CALL cp_fm_cholesky_decompose(inv%obj%p, info_out=info)
      IF (info .NE. 0) THEN
         WRITE (info_string, "(I6)") info
         CPABORT("Cholesky decomposition failed with info = "//info_string)
      END IF
      ! do inversion
      CALL cp_fm_cholesky_invert(inv%obj%p, info_out=info)
      IF (info .NE. 0) THEN
         WRITE (info_string, "(I6)") info
         CPABORT("Cholesky inversion failed with info = "//info_string)
      END IF
      ! copy the upper triangular part to lower to get the full matrix
      CALL mat_create(tmp, mat_nrows(inv), mat_ncols(inv), mat_blacs_env(inv))
      CALL mat_transpose(inv, tmp)
      CALL cp_fm_to_fm_triangular(tmp%obj%p, inv%obj%p, 'L')
      CALL mat_release(tmp)
   END SUBROUTINE mat_inv_cholesky_d

   SUBROUTINE mat_inv_cholesky_z(mat, inv)
      ! uses DPOTRF and DPOTRI for inversion using Cholesky
      ! decomposition
      TYPE(mat_z_obj), INTENT(IN) :: mat
      TYPE(mat_z_obj), INTENT(INOUT) :: inv

      INTEGER :: ii, jj, nrows, info
      CHARACTER(LEN=6) :: info_string
      TYPE(mat_z_obj) :: tmp

      CPASSERT(.NOT. is_same_obj(mat,inv))
      nrows = mat_nrows(mat)
      CPASSERT(mat_ncols(mat) .EQ. nrows)

      IF (ASSOCIATED(inv%obj)) THEN
         CPASSERT(mat_nrows(inv) .EQ. nrows)
         CPASSERT(mat_ncols(inv) .EQ. nrows)
         CPASSERT(mat_struct_equiv(mat, inv))
      ELSE
         CALL mat_create(inv, nrows, nrows, mat_blacs_env(mat))
      END IF
      ! copy matrix
      CALL cp_cfm_to_cfm_triangular(mat%obj%p, inv%obj%p, 'U')
      ! do cholesky decomposition
      CALL cp_cfm_cholesky_decompose(inv%obj%p, info_out=info)
      IF (info .NE. 0) THEN
         WRITE (info_string, "(I6)") info
         CPABORT("Cholesky decomposition failed with info = "//info_string)
      END IF
      ! do inversion
      CALL cp_cfm_cholesky_invert(inv%obj%p, info_out=info)
      IF (info .NE. 0) THEN
         WRITE (info_string, "(I6)") info
         CPABORT("Cholesky inversion failed with info = "//info_string)
      END IF
      ! copy the upper triangular part to lower to get the full matrix
      CALL mat_create(tmp, mat_nrows(inv), mat_ncols(inv), mat_blacs_env(inv))
      CALL mat_transpose(inv, tmp)
      CALL cp_cfm_to_cfm_triangular(tmp%obj%p, inv%obj%p, 'L')
      CALL mat_release(tmp)
   END SUBROUTINE mat_inv_cholesky_z

   SUBROUTINE mat_inv_lu_d(mat, inv)
      TYPE(mat_d_obj), INTENT(IN) :: mat
      TYPE(mat_d_obj), INTENT(INOUT) :: inv

      INTEGER :: nrows, info, lwork
      INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
      CHARACTER(LEN=6) :: info_string
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: work
      CPASSERT(.NOT. is_same_obj(mat,inv))
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
      CPASSERT(.NOT. is_same_obj(mat,inv))
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
            WRITE (my_unit_nr, FMT="(F12.7,2X)", ADVANCE="no") mat%obj%p(ii,jj)
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
            WRITE (my_unit_nr, FMT="(F12.7,1X,F12.7,3X)", ADVANCE="no") mat%obj%p(ii,jj)
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
      CPASSERT(.NOT. is_same_obj(X,Y))
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
      CPASSERT(.NOT. is_same_obj(X,Y))
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
      CASE ('C')
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

   SUBROUTINE mat_scale_d(mat, scalar)
      TYPE(mat_d_obj), INTENT(INOUT) :: mat
      REAL(KIND=dp), INTENT(IN) :: scalar
      mat%obj%p = scalar * mat%obj%p
   END SUBROUTINE mat_scale_d

   SUBROUTINE mat_scale_z(mat, scalar)
      TYPE(mat_z_obj), INTENT(INOUT) :: mat
      COMPLEX(KIND=dp), INTENT(IN) :: scalar
      mat%obj%p = scalar * mat%obj%p
   END SUBROUTINE mat_scale_z

   SUBROUTINE mat_associate_d(A, B)
      TYPE(mat_d_obj), INTENT(IN) :: A
      TYPE(mat_d_obj), INTENT(INOUT) :: B
      CALL mat_release(B)
      B%obj => A%obj
      CALL mat_retain(B)
   END SUBROUTINE mat_associate_d

   SUBROUTINE mat_associate_z(A, B)
      TYPE(mat_z_obj), INTENT(IN) :: A
      TYPE(mat_z_obj), INTENT(INOUT) :: B
      CALL mat_release(B)
      B%obj => A%obj
      CALL mat_retain(B)
   END SUBROUTINE mat_associate_z

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
      REAL(KIND=dp) :: res
      INTEGER :: nrows, ncols
      COMPLEX(KIND=dp), DIMENSION(1) :: work
      nrows = mat_nrows(mat)
      ncols = mat_ncols(mat)
      ! work is a dummy for Frobenius norm
      res = ZLANGE('F', nrows, ncols, mat%obj%p, nrows, work)
   END FUNCTION mat_norm_z

   PURE FUNCTION is_same_obj_d(A, B) RESULT(res)
      TYPE(mat_d_obj), INTENT(IN) :: A, B
      LOGICAL :: res
      res = .FALSE.
      IF (ASSOCIATED(A%obj) .AND. ASSOCIATED(B%obj)) THEN
         IF (A%obj%id_nr .EQ. B%obj%id_nr) THEN
            res = .TRUE.
         END IF
      END IF
   END FUNCTION is_same_obj_d

   FUNCTION is_same_obj_z(A, B) RESULT(res)
      TYPE(mat_z_obj), INTENT(IN) :: A, B
      LOGICAL :: res
      res = .FALSE.
      IF (ASSOCIATED(A%obj) .AND. ASSOCIATED(B%obj)) THEN
         IF (A%obj%id_nr .EQ. B%obj%id_nr) THEN
            res = .TRUE.
         END IF
      END IF
   END FUNCTION is_same_obj_z

   FUNCTION mat_trace_d(mat) RESULT(res)
      TYPE(mat_d_obj), INTENT(IN) :: mat
      REAL(KIND=dp) :: res
      INTEGER :: ii
      CPASSERT(ASSOCIATED(mat%obj))
      CPASSERT(ALLOCATED(mat%obj%p))
      res = 0.0_dp
      DO ii = 1, mat_nrows(mat)
         res = res + mat%obj%p(ii,ii)
      END DO
   END FUNCTION mat_trace_d

   FUNCTION mat_trace_z(mat) RESULT(res)
      TYPE(mat_z_obj), INTENT(IN) :: mat
      COMPLEX(KIND=dp) :: res
      INTEGER :: ii
      CPASSERT(ASSOCIATED(mat%obj))
      CPASSERT(ALLOCATED(mat%obj%p))
      res = 0.0_dp
      DO ii = 1, mat_nrows(mat)
         res = res + mat%obj%p(ii,ii)
      END DO
   END FUNCTION mat_trace_z

END MODULE matrix_types
