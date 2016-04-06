MODULE matrix_types

   USE kinds, ONLY: dp
   USE machine, ONLY: default_output_unit
   USE cp_fm_types, ONLY: cp_fm_type, &
                          cp_fm_create, &
                          cp_fm_to_fm_triangular, &
                          cp_fm_set_submatrix, &
                          cp_fm_release, &
                          cp_fm_get_info, &
                          cp_fm_get_submatrix, &
                          cp_fm_to_fm, &
                          cp_fm_to_fm_triangular
   USE cp_cfm_types, ONLY: cp_cfm_types, &
                           cp_cfm_create, &
                           cp_cfm_to_cfm_triangular, &
                           cp_cfm_set_submatrix, &
                           cp_cfm_release, &
                           cp_cfm_get_info, &
                           cp_cfm_get_submatrix, &
                           cp_cfm_to_cfm, &
                           cp_fm_to_cfm, &
                           cp_cfm_to_fm
   USE cp_blacs_env, ONLY: cp_blacs_env_type, &
                           get_blacs_info
   USE cp_fm_struct, ONLY: cp_fm_struct_type, &
                           cp_fm_struct_create, &
                           cp_fm_struct_get, &
                           cp_fm_struct_release, &
                           cp_fm_struct_equivalent
   USE cp_fm_basic_linalg, ONLY: cp_fm_gemm, &
                                 cp_fm_cholesky_decompose, &
                                 cp_fm_cholesky_invert, &
                                 cp_fm_lu_invert, &
                                 cp_fm_transpose, &
                                 cp_fm_geadd, &
                                 cp_fm_scale, &
                                 cp_fm_norm, &
                                 cp_fm_latra
   USE cp_cfm_basic_linalg, ONLY: cp_cfm_gemm, &
                                  cp_cfm_cholesky_decompose, &
                                  cp_cfm_cholesky_invert, &
                                  cp_cfm_lu_invert, &
                                  cp_cfm_transpose, &
                                  cp_cfm_geadd, &
                                  cp_cfm_scale, &
                                  cp_cfm_norm, &
                                  cp_cfm_latra
   USE cp_para_types, ONLY: cp_para_env_type

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
             mat_complex_to_real, &
             mat_release, &
             mat_scale, &
             mat_set, &
             mat_struct_equiv, &
             mat_symmetry, &
             mat_trace, &
             mat_transpose, &
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

   INTERFACE mat_set
      MODULE PROCEDURE mat_set_d
      MODULE PROCEDURE mat_set_z
   END INTERFACE mat_set

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

   INTERFACE mat_transpose
      MODULE PROCEDURE mat_transpose_d
      MODULE PROCEDURE mat_transpose_z
   END INTERFACE mat_transpose

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
      CALL cp_fm_struct_release(fmstruct)
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
      CALL cp_fm_struct_release(fmstruct)
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
            IF (ASSOCIATED(mat%obj%p)) THEN
               CALL cp_fm_release(mat%obj%p)
            END IF
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
            IF (ASSOCIATED(mat%obj%p)) THEN
               CALL cp_cfm_release(mat%obj%p)
            END IF
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
      CPASSERT(ASSOCIATED(mat_d%obj%p))
      CALL mat_release(mat_z)
      CALL mat_create(mat_z, &
                      mat_nrows(mat_d), &
                      mat_ncols(mat_d), &
                      mat_blacs_env(mat_d))
      CALL cp_fm_to_cfm(msourcer=mat_d%obj%p, &
                        mtarget=mat_z%obj%p)
      mat_z%obj%symmetry = mat_d%obj%symmetry
   END SUBROUTINE mat_real_to_complex

   SUBROUTINE mat_complex_to_real(mat_z, mat_r, mat_i)
      TYPE(mat_z_obj), INTENT(IN) :: mat_z
      TYPE(mat_d_obj), INTENT(INOUT), OPTIONAL :: mat_r, mat_i
      CPASSERT(ASSOCIATED(mat_z%obj))
      CPASSERT(ASSOCIATED(mat_z%obj%p))
      IF (PRESENT(mat_r)) THEN
         CALL mat_release(mat_r)
         CALL mat_create(mat_r, &
                         mat_nrows(mat_z), &
                         mat_ncols(mat_z), &
                         mat_blacs_env(mat_z))
         CALL cp_cfm_to_fm(msource=mat_z%obj%p, &
                           mtargetr=mat_r%obj%p)
      END IF
      IF (PRESENT(mat_i)) THEN
         CALL mat_release(mat_i)
         CALL mat_create(mat_i, &
                         mat_nrows(mat_z), &
                         mat_ncols(mat_z), &
                         mat_blacs_env(mat_z))
         CALL cp_cfm_to_fm(msource=mat_z%obj%p, &
                           mtargeti=mat_i%obj%p)
      END IF
   END SUBROUTINE mat_complex_to_real

   SUBROUTINE mat_mult_d(transA, transB, alpha, A, B, beta, C)
      TYPE(mat_d_obj), INTENT(IN) :: A, B
      TYPE(mat_d_obj), INTENT(INOUT) :: C
      CHARACTER, INTENT(IN) :: transA, transB
      REAL(KIND=dp), INTENT(IN) :: alpha, beta

      INTEGER :: nrows_A, ncols_A, ncols_B
      TYPE(cp_blacs_env_type), POINTER :: context

      CPASSERT(.NOT. is_same_obj(A,C))
      CPASSERT(.NOT. is_same_obj(B,C))

      nrows_A = mat_nrows(A)
      ncols_A = mat_ncols(A)
      ncols_B = mat_ncols(B)
      context => mat_blacs_env(A)

      CPASSERT(mat_nrows(B) .EQ. ncols_A)
      IF (ASSOCIATED(C%obj)) THEN
         CPASSERT(mat_nrows(C) .EQ. nrows_A)
         CPASSERT(mat_ncols(C) .EQ. ncols_B)
      ELSE
         CALL mat_create(C, nrows_A, ncols_B, context)
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
      TYPE(cp_blacs_env_type), POINTER :: context

      CPASSERT(.NOT. is_same_obj(A,C))
      CPASSERT(.NOT. is_same_obj(B,C))

      nrows_A = mat_nrows(A)
      ncols_A = mat_ncols(A)
      ncols_B = mat_ncols(B)
      context => mat_blacs_env(A)

      CPASSERT(mat_nrows(B) .EQ. ncols_A)
      IF (ASSOCIATED(C%obj)) THEN
         CPASSERT(mat_nrows(C) .EQ. nrows_A)
         CPASSERT(mat_ncols(C) .EQ. ncols_B)
      ELSE
         CALL mat_create(C, nrows_A, ncols_B, context)
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

      INTEGER :: nrows, info
      CHARACTER(LEN=6) :: info_string
      CPASSERT(.NOT. is_same_obj(mat,inv))
      nrows = mat_nrows(mat)
      IF (ASSOCIATED(inv%obj)) THEN
         CPASSERT(mat_nrows(inv) .EQ. nrows)
         CPASSERT(mat_ncols(inv) .EQ. nrows)
         CPASSERT(mat_struct_equiv(mat, inv))
      ELSE
         CALL mat_create(inv, nrows, nrows, mat_blacs_env(mat))
      END IF
      ! copy matrix
      CALL cp_fm_to_fm_triangular(mat%obj%p, inv%ob%p, 'U')
      ! do inversion
      CALL cp_fm_lu_invert(inv%obj%p, info)
      IF (info .NE. 0) THEN
         WRITE (info_string, FMT="(I6)") info
         CPABORT("LU inversion failed with info = "//info_string)
      END IF
   END SUBROUTINE mat_inv_lu_d

   SUBROUTINE mat_inv_lu_z(mat, inv)
      TYPE(mat_z_obj), INTENT(IN) :: mat
      TYPE(mat_z_obj), INTENT(INOUT) :: inv

      INTEGER :: nrows, info
      CHARACTER(LEN=6) :: info_string
      CPASSERT(.NOT. is_same_obj(mat,inv))
      nrows = mat_nrows(mat)
      IF (ASSOCIATED(inv%obj)) THEN
         CPASSERT(mat_nrows(inv) .EQ. nrows)
         CPASSERT(mat_ncols(inv) .EQ. nrows)
         CPASSERT(mat_struct_equiv(mat, inv))
      ELSE
         CALL mat_create(inv, nrows, nrows, mat_blacs_env(mat))
      END IF
      ! copy matrix
      CALL cp_cfm_to_cfm_triangular(mat%obj%p, inv%ob%p, 'U')
      ! do inversion
      CALL cp_cfm_lu_invert(inv%obj%p, info)
      IF (info .NE. 0) THEN
         WRITE (info_string, FMT="(I6)") info
         CPABORT("LU inversion failed with info = "//info_string)
      END IF
   END SUBROUTINE mat_inv_lu_z

   SUBROUTINE mat_set_d(mat, &
                        values, &
                        start_row, &
                        start_col, &
                        nrows, &
                        ncols, &
                        transpose, &
                        symmetry)
      TYPE(mat_d_obj), INTENT(INOUT) :: mat
      REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: values
      INTEGER, INTENT(IN), OPTIONAL :: start_row, start_col, nrows, ncols, symmetry
      CHARACTER, INTENT(IN), OPTIONAL :: transpose

      INTEGER :: nrows_global, ncols_global, &
                 my_start_row, my_start_col, my_nrows, my_ncols, &
                 my_symmetry
      LOGICAL :: do_transpose
      ! set defaults
      my_start_row = 1
      my_start_col = 1
      my_nrows = SIZE(values, 1)
      my_ncols = SIZE(values, 2)
      do_transpose = .FALSE.
      my_symmetry = MAT_GENERAL
      ! modify with user input
      IF (PRESENT(start_row)) my_start_row = start_row
      IF (PRESENT(start_col)) my_start_col = start_col
      IF (PRESENT(nrows)) my_nrows = nrows
      IF (PRESENT(ncols)) my_ncols = ncols
      IF (PRESENT(transpose)) do_transpose = (transpose .EQ. 'T')
      IF (PRESENT(symmetry)) my_symmetry = symmetry
      ! check matrix size compatible
      CPASSERT(ASSOCIATED(mat%obj))
      CPASSERT(ASSOCIATED(mat%obj%p))
      CALL cp_fm_get_info(matrix=mat%obj%p, &
                          nrow_global=nrows_global, &
                          ncol_global=ncols_global)
      IF (my_start_row+my_nrows-1 .GT. nrow_global) THEN
         CPABORT("limits exceed total number of rows")
      END IF
      IF (my_start_col+my_ncols-1 .GT. ncol_global) THEN
         CPABORT("limits exceed total number of cols")
      END IF
      ! copy data
      CALL cp_fm_set_submatrix(fm=matrix%obj%p, &
                               new_values=values, &
                               start_row=my_start_row, &
                               start_col=my_start_col, &
                               transpose=do_transpose)
   END SUBROUTINE mat_set_d

   SUBROUTINE mat_set_z(mat, &
                        values, &
                        start_row, &
                        start_col, &
                        nrows, &
                        ncols, &
                        transpose, &
                        symmetry)
      TYPE(mat_z_obj), INTENT(INOUT) :: mat
      COMPLEX(KIND=dp), DIMENSION(:,:), INTENT(IN) :: values
      INTEGER, INTENT(IN), OPTIONAL :: start_row, start_col, nrows, ncols, symmetry
      CHARACTER, INTENT(IN), OPTIONAL :: transpose

      INTEGER :: nrows_global, ncols_global, &
                 my_start_row, my_start_col, my_nrows, my_ncols, &
                 my_symmetry
      LOGICAL :: do_transpose
      ! set defaults
      my_start_row = 1
      my_start_col = 1
      my_nrows = SIZE(values, 1)
      my_ncols = SIZE(values, 2)
      do_transpose = .FALSE.
      my_symmetry = MAT_GENERAL
      ! modify with user input
      IF (PRESENT(start_row)) my_start_row = start_row
      IF (PRESENT(start_col)) my_start_col = start_col
      IF (PRESENT(nrows)) my_nrows = nrows
      IF (PRESENT(ncols)) my_ncols = ncols
      IF (PRESENT(transpose)) do_transpose = (transpose .EQ. 'T')
      IF (PRESENT(symmetry)) my_symmetry = symmetry
      ! check matrix size compatible
      CPASSERT(ASSOCIATED(mat%obj))
      CPASSERT(ASSOCIATED(mat%obj%p))
      CALL cp_cfm_get_info(matrix=mat%obj%p, &
                           nrow_global=nrows_global, &
                           ncol_global=ncols_global)
      IF (my_start_row+my_nrows-1 .GT. nrow_global) THEN
         CPABORT("limits exceed total number of rows")
      END IF
      IF (my_start_col+my_ncols-1 .GT. ncol_global) THEN
         CPABORT("limits exceed total number of cols")
      END IF
      ! copy data
      CALL cp_cfm_set_submatrix(fm=matrix%obj%p, &
                                new_values=values, &
                                start_row=my_start_row, &
                                start_col=my_start_col, &
                                transpose=do_transpose)
   END SUBROUTINE mat_set_z

   SUBROUTINE mat_read_d(mat, filename, blacs_env)
      TYPE(mat_d_obj), INTENT(INOUT) :: mat
      CHARACTER(LEN=*), INTENT(IN) :: filename
      TYPE(blacs_env_type), POINTER :: blacs_env

      INTEGER, PARAMETER :: UNIT_NR = 100
      INTEGER :: nrows, ncols, ii
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: content
      ! read into matrix
      OPEN(UNIT_NR, file=filename)
      READ (UNIT_NR, FMT=*) nrows, ncols
      ALLOCATE (content(nrows,ncols))
      content = 0.0_dp
      DO ii = 1, nrows
         READ (UNIT_NR, FMT=*) content(ii,:)
      END DO
      CLOSE(UNIT_NR)
      ! reallocate matrix
      IF (ASSOCIATED(mat%obj)) CALL mat_release(mat)
      CALL mat_create(mat, nrows, ncols, blacs_env)
      ! copy data
      CALL mat_set(mat, content)
      ! cleanup
      DEALLOCATE (content)
   END SUBROUTINE mat_read_d

   SUBROUTINE mat_read_z(mat, filename, blacs_env)
      TYPE(mat_z_obj), INTENT(INOUT) :: mat
      CHARACTER(LEN=*), INTENT(IN) :: filename

      INTEGER, PARAMETER :: UNIT_NR = 100
      INTEGER :: nrows, ncols, ii, jj
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: tmp
      COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: content
      ! read into matrix
      OPEN(UNIT_NR, file=filename)
      READ (UNIT_NR, FMT=*) nrows, ncols
      ALLOCATE (tmp(2*ncols))
      ALLOCATE (content(nrows,ncols))
      tmp = 0.0_dp
      DO ii = 1, nrows
         READ (UNIT_NR, FMT=*) tmp(:)
         DO jj = 1, ncols
            content(ii,jj) = CMPLX(tmp(2*jj-1), tmp(2*jj), KIND=dp)
         END DO
         tmp = 0.0_dp
      END DO
      DEALLOCATE(tmp)
      CLOSE(UNIT_NR)
      ! reallocate matrix
      IF (ASSOCIATED(mat%obj)) CALL mat_release(mat)
      CALL mat_create(mat, nrows, ncols, blacs_env)
      ! copy data
      CALL mat_set(mat, content)
      ! cleanup
      DEALLOCATE(content)
   END SUBROUTINE mat_read_z

   SUBROUTINE mat_write_d(mat, filename)
      TYPE(mat_d_obj), INTENT(IN) :: mat
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename

      INTEGER :: my_unit_nr, nrows, ncols, ii, jj, &
                 mepos, ionode
      TYPE(cp_blacs_env_type), POINTER :: context
      TYPE(cp_para_env_type), POINTER :: para_env
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: content

      CPASSERT(ASSOCIATED(mat%obj))
      CPASSERT(ASSOCIATED(mat%obj%p))
      nrows = mat_nrows(mat)
      ncols = mat_ncols(mat)
      context => mat_blacs_env(mat)
      CALL get_blacs_info(blacs_env=context, &
                          para_env=para_env)
      mepos = para_env%mepos
      ionode = para_env%ionode
      CALL cp_fm_get_submatrix(fm=matrix%obj%p, &
                               target_m=context)
      ! only ionode should write the matrix
      IF (mepos .EQ. ionode) THEN
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
               WRITE (my_unit_nr, FMT="(F12.7,2X)", ADVANCE="no") content(ii,jj)
            END DO
            WRITE (my_unit_nr, *) ""
         END DO
         ! close file if we need to
         IF (PRESENT(filename)) THEN
            CLOSE(my_unit_nr)
         END IF
      END IF
   END SUBROUTINE mat_write_d

   SUBROUTINE mat_write_z(mat, filename)
      TYPE(mat_z_obj), INTENT(IN) :: mat
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: filename

      INTEGER :: my_unit_nr, nrows, ncols, ii, jj, &
                 mepos, ionode
      TYPE(cp_blacs_env_type), POINTER :: context
      TYPE(cp_para_env_type), POINTER :: para_env
      COMPLEX(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: content

      CPASSERT(ASSOCIATED(mat%obj))
      CPASSERT(ASSOCIATED(mat%obj%p))
      nrows = mat_nrows(mat)
      ncols = mat_ncols(mat)
      context => mat_blacs_env(mat)
      CALL get_blacs_info(blacs_env=context, &
                          para_env=para_env)
      mepos = para_env%mepos
      ionode = para_env%ionode
      CALL cp_cfm_get_submatrix(fm=matrix%obj%p, &
                                target_m=context)
      ! only ionode should write the matrix
      IF (mepos .EQ. ionode) THEN
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
               WRITE (my_unit_nr, FMT="(F12.7,1X,F12.7,3X)", ADVANCE="no") content(ii,jj)
            END DO
            WRITE (my_unit_nr, *) ""
         END DO
         ! close file if we need to
         IF (PRESENT(filename)) THEN
            CLOSE(my_unit_nr)
         END IF
      END IF
   END SUBROUTINE mat_write_z

   PURE FUNCTION mat_nrows_d(mat) RESULT(res)
      TYPE(mat_d_obj), INTENT(IN) :: mat
      INTEGER :: res
      CALL cp_fm_get_info(matrix=mat%obj%p, &
                          nrow_global=res)
   END FUNCTION mat_nrows_d

   PURE FUNCTION mat_nrows_z(mat) RESULT(res)
      TYPE(mat_z_obj), INTENT(IN) :: mat
      INTEGER :: res
      CALL cp_cfm_get_info(matrix=mat%obj%p, &
                           nrow_global=res)
   END FUNCTION mat_nrows_z

   PURE FUNCTION mat_ncols_d(mat) RESULT(res)
      TYPE(mat_d_obj), INTENT(IN) :: mat
      INTEGER :: res
      CALL cp_fm_get_info(matrix=mat%obj%p, &
                          ncol_global=res)
   END FUNCTION mat_ncols_d

   PURE FUNCTION mat_ncols_z(mat) RESULT(res)
      TYPE(mat_z_obj), INTENT(IN) :: mat
      INTEGER :: res
      CALL cp_cfm_get_info(matrix=mat%obj%p, &
                           ncol_global=res)
   END FUNCTION mat_ncols_z

   SUBROUTINE mat_axpy_d(a, transX, X, Y)
      ! computes Y = aX^(T) + Y
      REAL(KIND=dp), INTENT(IN) :: a
      CHARACTER, INTENT(IN) :: transX
      TYPE(mat_d_obj), INTENT(IN) :: X
      TYPE(mat_d_obj), INTENT(INOUT) :: Y
      CPASSERT(.NOT. is_same_obj(X,Y))
      CALL cp_fm_geadd(a, transX, X%obj%p, 1.0_dp, Y%obj%p)
   END SUBROUTINE mat_axpy_d

   SUBROUTINE mat_axpy_z(a, transX, X, Y)
      ! computes Y = aX^(T) + Y
      COMPLEX(KIND=dp), INTENT(IN) :: a
      CHARACTER, INTENT(IN) :: transX
      TYPE(mat_z_obj), INTENT(IN) :: X
      TYPE(mat_z_obj), INTENT(INOUT) :: Y
      CPASSERT(.NOT. is_same_obj(X,Y))
      CALL cp_cfm_geadd(a, transX, X%obj%p, (1.0_dp,0.0_dp), Y%obj%p)
   END SUBROUTINE mat_axpy_z

   SUBROUTINE mat_zero_d(mat)
      TYPE(mat_d_obj), INTENT(INOUT) :: mat
      ! assume assoicated obj always leads to allocated obj%p, as no
      ! module method should allocate obj but not also allocate obj%p
      IF (ASSOCIATED(mat%obj)) THEN
         CALL cp_fm_scale(0.0_dp, mat%obj%p)
      END IF
   END SUBROUTINE mat_zero_d

   SUBROUTINE mat_zero_z(mat)
      TYPE(mat_z_obj), INTENT(INOUT) :: mat
      ! assume assoicated obj always leads to allocated obj%p, as no
      ! module method should allocate obj but not also allocate obj%p
      IF (ASSOCIATED(mat%obj)) THEN
         CALL cp_cfm_scale((0.0_dp,0.0_dp), mat%obj%p)
      END IF
   END SUBROUTINE mat_zero_z

   SUBROUTINE mat_copy_d(A, B)
      TYPE(mat_d_obj), INTENT(IN) :: A
      TYPE(mat_d_obj), INTENT(INOUT) :: B
      IF (ASSOCIATED(B%obj)) CALL mat_release(B)
      CALL mat_create(B, mat_nrows(A), mat_ncols(A), mat_blacs_env(A))
      CALL cp_fm_to_fm(A%obj%p, B%obj%p)
   END SUBROUTINE mat_copy_d

   SUBROUTINE mat_copy_z(A, B)
      TYPE(mat_z_obj), INTENT(IN) :: A
      TYPE(mat_z_obj), INTENT(INOUT) :: B
      IF (ASSOCIATED(B%obj)) CALL mat_release(B)
      CALL mat_create(B, mat_nrows(A), mat_ncols(A), mat_blacs_env(A))
      CALL cp_cfm_to_cfm(A%obj%p, B%obj%p)
   END SUBROUTINE mat_copy_z

   SUBROUTINE mat_scale_d(mat, scalar)
      TYPE(mat_d_obj), INTENT(INOUT) :: mat
      REAL(KIND=dp), INTENT(IN) :: scalar
      IF (ASSOCIATED(mat%obj)) THEN
         CALL cp_fm_scale(scalar, mat%obj%p)
      END IF
   END SUBROUTINE mat_scale_d

   SUBROUTINE mat_scale_z(mat, scalar)
      TYPE(mat_z_obj), INTENT(INOUT) :: mat
      COMPLEX(KIND=dp), INTENT(IN) :: scalar
      IF (ASSOCIATED(mat%obj)) THEN
         CALL cp_cfm_scale(scalar, mat%obj%p)
      END IF
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
      res = cp_fm_norm(mat%obj%p, 'F')
   END FUNCTION mat_norm_d

   FUNCTION mat_norm_z(mat) RESULT(res)
      TYPE(mat_z_obj), INTENT(IN) :: mat
      REAL(KIND=dp) :: res
      res = cp_cfm_norm(mat%obj%p, 'F')
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
      CPASSERT(ASSOCIATED(mat%obj%p))
      res = cp_fm_latra(mat%obj%p)
   END FUNCTION mat_trace_d

   FUNCTION mat_trace_z(mat) RESULT(res)
      TYPE(mat_z_obj), INTENT(IN) :: mat
      COMPLEX(KIND=dp) :: res
      INTEGER :: ii
      CPASSERT(ASSOCIATED(mat%obj))
      CPASSERT(ASSOCIATED(mat%obj%p))
      res = cp_cfm_latra(mat%obj%p)
   END FUNCTION mat_trace_z

   SUBROUTINE mat_transpose_d(mat, mat_t)
      TYPE(mat_d_obj), INTENT(IN) :: mat
      TYPE(mat_d_obj), INTENT(INOUT) :: mat_t
      CPASSERT(ASSOCIATED(mat%obj))
      CPASSERT(ASSOCIATED(mat%obj%p))
      CALL mat_release(mat_t)
      CALL mat_create(mat_t, &
                      mat_nrows(mat), &
                      mat_ncols(mat), &
                      mat_blacs_env(mat))
      CALL cp_fm_transpose(mat%obj%p, mat_t%obj%p)
   END SUBROUTINE mat_transpose_d

   SUBROUTINE mat_transpose_z(mat, mat_t, trans)
      TYPE(mat_z_obj), INTENT(IN) :: mat
      TYPE(mat_z_obj), INTENT(INOUT) :: mat_t
      CHARACTER, INTENT(IN), OPTIONAL :: trans
      CHARACTER :: my_trans
      CPASSERT(ASSOCIATED(mat%obj))
      CPASSERT(ASSOCIATED(mat%obj%p))
      my_trans = 'T'
      IF (PRESENT(trans)) my_trans = trans
      CALL mat_release(mat_t)
      CALL mat_create(mat_t, &
                      mat_nrows(mat), &
                      mat_ncols(mat), &
                      mat_blacs_env(mat))
      CALL cp_cfm_transpose(mat%obj%p, my_trans, mat_t%obj%p)
   END SUBROUTINE mat_transpose_z

END MODULE matrix_types
