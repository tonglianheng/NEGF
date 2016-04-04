PROGRAM test_GR
!test gL, Sigma, Gamma, G

USE kinds, ONLY: dp
USE matrix_types
USE sancho_method
USE negf_env_types

IMPLICIT NONE
REAL(KIND=dp),PARAMETER           :: Tol = 1e-4
REAL(KIND=dp),PARAMETER           :: eta = 1e-4
REAL(KIND=dp),PARAMETER           :: my_energy = 0
COMPLEX(KIND=dp),PARAMETER        :: E   = (my_energy,eta)


TYPE(mat_z_obj) :: mat_gL, mat_GR
TYPE(mat_d_obj) :: mat_hL,mat_tL,mat_shL,mat_stL
TYPE(negf_env_obj) :: my_negf_env

! read-out the matrices from input/XXX already in a form of a matrix obj
CALL mat_read(mat_hL, "input/simple/hL.dat")
WRITE (*,*) 
CALL mat_write(mat_hL)
CALL mat_read(mat_tL, "input/simple/tL.dat")
WRITE (*,*) 
CALL mat_write(mat_tL)
CALL mat_read(mat_shL, "input/simple/ShL.dat")
WRITE (*,*) 
CALL mat_write(mat_shL)
CALL mat_read(mat_stL, "input/simple/StL.dat")
WRITE (*,*) 
CALL mat_write(mat_stL)
CALL mat_read(mat_stL, "input/simple/StL.dat")
WRITE (*,*) 
CALL mat_write(mat_stL)

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

!----------------------------------------------------------!
CALL negf_env_create(negf_env = my_negf_env, nterminals = 2, eps_E = eta, eps_Sancho = Tol)

!on-site leads
CALL mat_read(my_negf_env%obj%H_L_onsite(1),"input/simple/hL.dat")
CALL mat_read(my_negf_env%obj%H_L_onsite(2),"input/simple/hR.dat")
CALL mat_read(my_negf_env%obj%S_L_onsite(1),"input/simple/ShL.dat")
CALL mat_read(my_negf_env%obj%S_L_onsite(2),"input/simple/ShR.dat")
!hopping leads
CALL mat_read(my_negf_env%obj%H_L_hopping(1),"input/simple/tL.dat")
CALL mat_read(my_negf_env%obj%H_L_hopping(2),"input/simple/tR.dat")
CALL mat_read(my_negf_env%obj%S_L_hopping(1),"input/simple/StL.dat")
CALL mat_read(my_negf_env%obj%S_L_hopping(2),"input/simple/StR.dat")
!leads -- system
CALL mat_read(my_negf_env%obj%H_LS(1),"input/simple/TLS.dat")
CALL mat_read(my_negf_env%obj%H_LS(2),"input/simple/TSR.dat")
CALL mat_read(my_negf_env%obj%S_LS(1),"input/simple/STLS.dat")
CALL mat_read(my_negf_env%obj%S_LS(2),"input/simple/STSR.dat")
!system
CALL mat_read(my_negf_env%obj%H_S,"input/simple/HS.dat")
CALL mat_read(my_negf_env%obj%S_S,"input/simple/SHS.dat")

 
CALL calc_GR(negf_env = my_negf_env, energy = my_energy, GRetarded = mat_GR)

WRITE(*,*) "GR = "
CALL mat_write(mat_GR)

!CALL quad(negf_env, a, b, tol, trace)

END PROGRAM test_GR

