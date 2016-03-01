PROGRAM test_matrix
  USE test_utils, ONLY: check_pass, &
                        test_init, &
                        test_finalise
  USE kinds, ONLY: dp
  USE matrix_types

  IMPLICIT NONE

  INTEGER :: unit_nr, ii, jj
  INTEGER, PARAMETER :: N1 = 4, &
                        N2 = 4, &
                        M2 = 3
  REAL(KIND=dp) :: norm2
  REAL(KIND=dp), DIMENSION(N1,N1) :: eigvectors, mat1_diag, mat1_identity
  REAL(KIND=dp), DIMENSION(N1,N1) :: mat1_data, tmp1_data
  REAL(KIND=dp), DIMENSION(N2,M2) :: mat2_data
  REAL(KIND=dp), DIMENSION(N1,M2) :: mat3_data
  COMPLEX(KIND=dp), DIMENSION(N1,N1) :: zmat1_identity
  TYPE(mat_d_obj) :: mat1, mat1_inv, mat1_copy, mat2, mat3
  TYPE(mat_d_obj) :: identity, tmp
  TYPE(mat_z_obj) :: zmat1, zmat1_inv, zmat1_copy, zmat2, zmat3
  TYPE(mat_z_obj) :: zidentity, ztmp

  ! initialise test utilities
!  CALL test_init(tol=1.E-15_dp)
  CALL test_init()

  ! constructing a positive definite matrix for mat1_data
  CALL RANDOM_SEED()
  CALL RANDOM_NUMBER(eigvectors)
  eigvectors = 10.0_dp * eigvectors
  mat1_diag = 0.0_dp
  DO ii = 1, N1
     CALL RANDOM_NUMBER(mat1_diag(ii,ii))
  END DO
  CALL DGEMM('N', 'T', N1, N1, N1, &
             1.0_dp, mat1_diag, N1, eigvectors, N1, &
             0.0_dp, tmp1_data, N1)
  CALL DGEMM('N', 'N', N1, N1, N1, &
             1.0_dp, eigvectors, N1, tmp1_data, N1, &
             0.0_dp, mat1_data, N1)

  ! constructing an identity matrix compatible with mat1_data
  mat1_identity = 0.0_dp
  zmat1_identity = (0.0_dp,0.0_dp)
  DO ii = 1, N1
     mat1_identity(ii,ii) = 1.0_dp
     zmat1_identity(ii,ii) = (1.0_dp,0.0_dp)
  END DO

  ! constructing a random matrix for mat2_data
  CALL RANDOM_NUMBER(mat2_data)

  ! write matrices to files
  OPEN(UNIT=unit_nr, FILE="test_matrix_mat1.txt", ACTION='WRITE')
  WRITE (unit_nr, *) N1, N1
  DO ii = 1, N1
     WRITE (unit_nr, *) mat1_data(ii,:)
  END DO
  CLOSE(unit_nr)

  ! write matrices to files
  OPEN(UNIT=unit_nr, FILE="test_matrix_mat2.txt", ACTION='WRITE')
  WRITE (unit_nr, *) N2, M2
  DO ii = 1, N2
     WRITE (unit_nr, *) mat2_data(ii,:)
  END DO
  CLOSE(unit_nr)

  ! calculate Frobenius norm of mat2
  norm2 = 0.0_dp
  DO jj = 1, M2
     DO ii = 1, N2
        norm2 = norm2 + mat2_data(ii,jj)**2
     END DO
  END DO
  norm2 = SQRT(norm2)

  ! do checks
  WRITE (*,*) "==================== check read ===================="
  WRITE (*,*) "correct mat1"
  DO ii = 1, N1
     WRITE (*, *) mat1_data(ii,:)
  END DO
  CALL mat_read(mat1, "./test_matrix_mat1.txt")
  WRITE (*,*) "mat1"
  CALL mat_write(mat1)
  WRITE (*,*) "correct mat2"
  DO ii = 1, N1
     WRITE (*, *) mat2_data(ii,:)
  END DO
  CALL mat_read(mat2, "./test_matrix_mat2.txt")
  WRITE (*,*) "mat2"
  CALL mat_write(mat2)

  WRITE (*,*) "==================== check nrows ===================="
  WRITE (*,*) "correct number of rows for mat2 = ", N2
  WRITE (*,*) "number of rows for mat2 = ", mat_nrows(mat2)
  CALL check_pass(mat_nrows(mat2), N2)

  WRITE (*,*) "==================== check ncols ===================="
  WRITE (*,*) "correct number of cols for mat2 = ", M2
  WRITE (*,*) "number of cols for mat2 = ", mat_ncols(mat2)
  CALL check_pass(mat_ncols(mat2), M2)

  WRITE (*,*) "==================== check copy ===================="
  CALL mat_copy(mat1, mat1_copy)
  CALL mat_write(mat1)
  CALL mat_write(mat1_copy)

  WRITE (*,*) "==================== check norm ===================="
  WRITE (*,*) "correct norm for mat2 = ", norm2
  WRITE (*,*) "calculated norm for mat2 = ", mat_norm(mat2)
  CALL check_pass(mat_norm(mat2), norm2)

  WRITE (*,*) "==================== check axpy ===================="
  CALL mat_release(tmp)
  CALL mat_create(tmp, mat_nrows(mat2), mat_ncols(mat2))
  CALL mat_write(mat2)
  CALL mat_axpy(2.0_dp, 'N', mat2, tmp)
  CALL mat_write(tmp)
  CALL mat_axpy(-2.0_dp, 'N', mat2, tmp)
  CALL check_pass(mat_norm(tmp), 0.0_dp)

  WRITE (*,*) "==================== check associate ===================="
  CALL mat_associate(mat1, mat1_copy)
  CALL mat_write(mat1)
  CALL mat_write(mat1_copy)
  CALL mat_copy(mat1_copy, tmp)
  CALL mat_axpy(-1.0_dp, 'N', mat1, tmp)
  CALL check_pass(mat_norm(tmp), 0.0_dp)

  WRITE (*,*) "==================== check zero ===================="
  CALL mat_copy(mat1, mat1_copy)
  CALL mat_write(mat1_copy)
  CALL mat_zero(mat1_copy)
  CALL mat_write(mat1_copy)
  CALL check_pass(mat_norm(mat1_copy), 0.0_dp)

  WRITE (*,*) "==================== check multiply ===================="
  CALL mat_write(mat1)
  CAlL mat_write(mat2)
  CALL DGEMM('N', 'N', N1, M2, N1, &
             1.0_dp, mat1_data, N1, mat2_data, N2, &
             0.0_dp, mat3_data, N1)
  CALL mat_release(tmp)
  CALL mat_create(tmp, N1, M2, content=mat3_data)
  WRITE (*,*) "correct mat1 * mat2"
  CALL mat_write(tmp)
  CALL mat_mult("N", "N", 1.0_dp, mat1, mat2, 0.0_dp, mat3)
  WRITE (*,*) "calculated mat1 * mat2"
  CALL mat_write(mat3)
  CALL mat_axpy(-1.0_dp, 'N', mat3, tmp)
  CALL check_pass(mat_norm(tmp), 0.0_dp)

  WRITE (*,*) "==================== check inverse cholesky ===================="
  CALL mat_write(mat1)
  CALL mat_inv_cholesky(mat1, mat1_inv)
  WRITE (*,*) "calculated inverse"
  CALL mat_write(mat1_inv)
  CALL mat_mult("N", "N", 1.0_dp, mat1, mat1_inv, 0.0_dp, identity)
  WRITE (*,*) "mat * inverse"
  CALL mat_write(identity)
  CALL mat_release(tmp)
  CALL mat_create(tmp, N1, N1, content=mat1_identity)
  CALL mat_axpy(-1.0_dp, 'N', identity, tmp)
  CALL check_pass(mat_norm(tmp), 0.0_dp)

  WRITE (*,*) "==================== check inverse LU ===================="
  CALL mat_write(mat1)
  CALL mat_inv_lu(mat1, mat1_inv)
  WRITE (*,*) "calculated inverse"
  CALL mat_write(mat1_inv)
  CALL mat_mult("N", "N", 1.0_dp, mat1, mat1_inv, 0.0_dp, identity)
  CALL mat_write(identity)
  CALL mat_release(tmp)
  CALL mat_create(tmp, N1, N1, content=mat1_identity)
  CALL mat_axpy(-1.0_dp, 'N', identity, tmp)
  CALL check_pass(mat_norm(tmp), 0.0_dp)

  ! ------------------------------------------------------------------------
  ! Complex Tests
  ! ------------------------------------------------------------------------

  WRITE (*,*) ""
  WRITE (*,*) "========================================================================"
  WRITE (*,*) "* COMPLEX TESTS"
  WRITE (*,*) "========================================================================"
  WRITE (*,*) ""

  WRITE (*,*) "==================== check real_to_complex ===================="
  CALL mat_real_to_complex(mat1, zmat1)
  CALL mat_write(mat1)
  CALL mat_write(zmat1)
  CALL mat_real_to_complex(mat2, zmat2)
  CALL mat_write(mat2)
  CALL mat_write(zmat2)

  WRITE (*,*) "==================== check nrows ===================="
  WRITE (*,*) "correct number of rows for zmat2 = ", N2
  WRITE (*,*) "number of rows for zmat2 = ", mat_nrows(zmat2)
  CALL check_pass(mat_nrows(zmat2), N2)

  WRITE (*,*) "==================== check ncols ===================="
  WRITE (*,*) "correct number of cols for zmat2 = ", M2
  WRITE (*,*) "number of cols for zmat2 = ", mat_ncols(zmat2)
  CALL check_pass(mat_ncols(zmat2), M2)

  WRITE (*,*) "==================== check copy ===================="
  CALL mat_copy(zmat1, zmat1_copy)
  CALL mat_write(zmat1)
  CALL mat_write(zmat1_copy)

  WRITE (*,*) "==================== check norm ===================="
  WRITE (*,*) "correct norm for zmat2 = ", norm2
  WRITE (*,*) "calculated norm for zmat2 = ", mat_norm(zmat2)
  CALL check_pass(mat_norm(zmat2), norm2)

  WRITE (*,*) "==================== check axpy ===================="
  CALL mat_release(ztmp)
  CALL mat_create(ztmp, mat_nrows(zmat2), mat_ncols(zmat2))
  CALL mat_write(zmat2)
  CALL mat_axpy((2.0_dp,0.0_dp), 'N', zmat2, ztmp)
  CALL mat_write(ztmp)
  CALL mat_axpy((-2.0_dp,0.0_dp), 'N', zmat2, ztmp)
  CALL check_pass(mat_norm(ztmp), 0.0_dp)

  WRITE (*,*) "==================== check associate ===================="
  CALL mat_associate(zmat1, zmat1_copy)
  CALL mat_write(zmat1)
  CALL mat_write(zmat1_copy)
  CALL mat_copy(zmat1_copy, ztmp)
  CALL mat_axpy((-1.0_dp,0.0_dp), 'N', zmat1, ztmp)
  CALL check_pass(mat_norm(ztmp), 0.0_dp)

  WRITE (*,*) "==================== check zero ===================="
  CALL mat_copy(zmat1, zmat1_copy)
  CALL mat_write(zmat1_copy)
  CALL mat_zero(zmat1_copy)
  CALL mat_write(zmat1_copy)
  CALL check_pass(mat_norm(zmat1_copy), 0.0_dp)

  WRITE (*,*) "==================== check multiply ===================="
  CALL mat_write(zmat1)
  CAlL mat_write(zmat2)
  CALL DGEMM('N', 'N', N1, M2, N1, &
             1.0_dp, mat1_data, N1, mat2_data, N2, &
             0.0_dp, mat3_data, N1)
  CALL mat_release(tmp)
  CALL mat_create(tmp, N1, M2, content=mat3_data)
  CALL mat_real_to_complex(tmp, ztmp)
  WRITE (*,*) "correct zmat1 * zmat2"
  CALL mat_write(ztmp)
  CALL mat_mult("N", "N", (1.0_dp,0.0_dp), zmat1, zmat2, (0.0_dp,0.0_dp), zmat3)
  WRITE (*,*) "calculated zmat1 * zmat2"
  CALL mat_write(zmat3)
  CALL mat_axpy((-1.0_dp,0.0_dp), 'N', zmat3, ztmp)
  CALL check_pass(mat_norm(ztmp), 0.0_dp)

  WRITE (*,*) "==================== check inverse cholesky ===================="
  CALL mat_write(zmat1)
  CALL mat_inv_cholesky(zmat1, zmat1_inv)
  WRITE (*,*) "calculated inverse"
  CALL mat_write(zmat1_inv)
  CALL mat_mult("N", "N", (1.0_dp,0.0_dp), zmat1, zmat1_inv, (0.0_dp,0.0_dp), zidentity)
  WRITE (*,*) "mat * inverse"
  CALL mat_write(zidentity)
  CALL mat_release(ztmp)
  CALL mat_create(ztmp, N1, N1, content=zmat1_identity)
  CALL mat_axpy((-1.0_dp,0.0_dp), 'N', zidentity, ztmp)
  CALL check_pass(mat_norm(ztmp), 0.0_dp)

  WRITE (*,*) "==================== check inverse LU ===================="
  CALL mat_write(zmat1)
  CALL mat_inv_lu(zmat1, zmat1_inv)
  WRITE (*,*) "calculated inverse"
  CALL mat_write(zmat1_inv)
  CALL mat_mult("N", "N", (1.0_dp,0.0_dp), zmat1, zmat1_inv, (0.0_dp,0.0_dp), zidentity)
  CALL mat_write(zidentity)
  CALL mat_release(ztmp)
  CALL mat_create(ztmp, N1, N1, content=zmat1_identity)
  CALL mat_axpy((-1.0_dp,0.0_dp), 'N', zidentity, ztmp)
  CALL check_pass(mat_norm(ztmp), 0.0_dp)

  ! finish tests
  CALL test_finalise()

  ! clean up
  CALL mat_release(mat1)
  CALL mat_release(mat2)
  CALL mat_release(mat3)
  CALL mat_release(tmp)
  CALL mat_release(mat1_inv)
  CALL mat_release(identity)

  CALL mat_release(zmat1)
  CALL mat_release(zmat2)
  CALL mat_release(zmat3)
  CALL mat_release(ztmp)
  CALL mat_release(zmat1_inv)
  CALL mat_release(zidentity)

END PROGRAM test_matrix
