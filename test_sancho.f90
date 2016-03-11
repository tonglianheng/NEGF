PROGRAM test_sancho
  USE kinds, ONLY: dp
  USE matrix_types
  USE sancho_method

  IMPLICIT NONE

  REAL(KIND=dp),PARAMETER           :: Tol = 1e-4
  REAL(KIND=dp),PARAMETER           :: eta = 1e-4
  COMPLEX(KIND=dp),PARAMETER        :: E   = (0,eta)
  INTEGER, PARAMETER                :: N1  = 2

  REAL(KIND=dp), DIMENSION(N1,N1) :: hL, tL, shL, stL

  TYPE(mat_z_obj) :: gL
  TYPE(mat_d_obj) :: mat1, mat2, mat3, mat4

  ! constructting input matrices
  CALL RANDOM_SEED()
  CALL RANDOM_NUMBER(hL)
  CALL RANDOM_NUMBER(tL)
  CALL RANDOM_NUMBER(shL)
  CALL RANDOM_NUMBER(stL)

  hL = 10.0_dp * hL
  tL = 10.0_dp * tL
  ! end: contructing input matrices
  
  CALL mat_create(mat1, nrows=N1, ncols=N1, content=hL)
  CALL mat_create(mat2, nrows=N1, ncols=N1, content=tL)
  CALL mat_create(mat3, nrows=N1, ncols=N1, content=shL)
  CALL mat_create(mat4, nrows=N1, ncols=N1, content=stL)
  
  CALL surface_GR_sancho(GR_surface= gL, &
                         energy    = E, &
                         H_onsite  = mat1, &
                         H_hopping = mat2, &
                         S_onsite  = mat3, &
                         S_hopping = mat4, &
                         tolerance = Tol)

  WRITE (*,*) "hL = "
  CALL mat_write(mat1)
  WRITE (*,*) "tL = "
  CALL mat_write(mat2)
  WRITE (*,*) "shL = "
  CALL mat_write(mat3)
  WRITE (*,*) "stL = "
  CALL mat_write(mat4)
END PROGRAM test_sancho
