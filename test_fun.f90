PROGRAM test_fun
!test gL, Sigma, Gamma, G

USE kinds, ONLY: dp
USE matrix_types
USE sancho_method

IMPLICIT NONE
REAL(KIND=dp),PARAMETER           :: Tol = 1e-4
REAL(KIND=dp),PARAMETER           :: eta = 1e-4
COMPLEX(KIND=dp),PARAMETER        :: E   = (0,eta)


TYPE(mat_z_obj) :: mat_gL
TYPE(mat_d_obj) :: mat_hL,mat_tL,mat_shL,mat_stL


!initialize test utilities -- add it in the end

! constructting input matrices from random numbers
!CALL RANDOM_SEED()
!CALL RANDOM_NUMBER(hL)
!CALL RANDOM_NUMBER(tL)
!CALL RANDOM_NUMBER(shL)
!CALL RANDOM_NUMBER(stL)
!hL = 10.0_dp * hL
!tL = 10.0_dp * tL
!matrix --> object
!CALL mat_create(mat1, nrows=N1, ncols=N1, content=hL)
!CALL mat_create(mat2, nrows=N1, ncols=N1, content=tL)
!CALL mat_create(mat3, nrows=N1, ncols=N1, content=shL)
!CALL mat_create(mat4, nrows=N1, ncols=N1, content=stL)
!end: matrix --> object
! end: contructing input matrices from  random numbers


! read-out the matrices from input/XXX already in a form of a matrix obj
CALL mat_read(mat_hL, "input/simple/hL.dat")
WRITE (*,*) "mat_hL"
CALL mat_read(mat_tL, "input/simple/tL.dat")
WRITE (*,*) "mat_tL"
CALL mat_read(mat_shL, "input/simple/ShL.dat")
WRITE (*,*) "mat_shL"
CALL mat_read(mat_stL, "input/simple/StL.dat")
WRITE (*,*) "mat_stL"
! end: read-out the matrices from input/XXX ...

 CALL surface_GR_sancho(GR_surface= mat_gL, &
                        energy    = E, &
                        H_onsite  = mat_hL, &
                        H_hopping = mat_tL, &
                        S_onsite  = mat_shL, &
                        S_hopping = mat_stL, &
                        tolerance = Tol)

WRITE (*,*) "hL = "
CALL mat_write(mat_hL)
WRITE (*,*) "tL = "
CALL mat_write(mat_tL)
WRITE (*,*) "shL = "
CALL mat_write(mat_shL)
WRITE (*,*) "stL = "
CALL mat_write(mat_stL)
WRITE (*,*) "gs="
CALL mat_write(mat_gL)



END PROGRAM test_fun
