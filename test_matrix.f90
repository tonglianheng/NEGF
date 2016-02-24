PROGRAM test_matrix
  USE matrix_types
  USE kinds
  
  IMPLICIT NONE

  TYPE(mat_d_obj) :: HS, TSL, HS_TSL, HS_inv, identity
  
  CALL mat_read(HS, "../input/Daijiro/HS.dat")
  CALL mat_write(HS)

  CALL mat_read(TSL, "../input/Daijiro/TSL.dat")
  CALL mat_write(TSL)

  CALL mat_mult("N", "N", 1.0_dp, HS, TSL, 0.0_dp, HS_TSL)
  CALL mat_write(HS_TSL)

  CALL mat_inv_cholesky(HS, HS_inv)
  CALL mat_write(HS_inv)
  CALL mat_mult("N", "N", 1.0_dp, HS, HS_inv, 0.0_dp, identity)
  CALL mat_write(identity)
  
  CALL mat_inv_lu(HS, HS_inv)
  CALL mat_write(HS_inv)
  CALL mat_mult("N", "N", 1.0_dp, HS, HS_inv, 0.0_dp, identity)
  CALL mat_write(identity)

END PROGRAM test_matrix
