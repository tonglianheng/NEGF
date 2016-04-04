PROGRAM test_quad
USE test_utils, ONLY: check_pass, &
                       test_init,&
                       test_finalise
USE kinds, ONLY: dp
USE matrix_types
USE sancho_method
USE negf_env_types
USE integral

IMPLICIT NONE

REAL(KIND=dp)           :: Tol
REAL(KIND=dp)           :: eta
REAL(KIND=dp)           :: my_energy
COMPLEX(KIND=dp)        :: E
REAL(KIND=dp)           :: my_a, my_b

TYPE(mat_z_obj) :: mat_gL, mat_GR, my_integral,matlab_integral
TYPE(mat_d_obj) :: mat_hL,mat_tL,mat_shL,mat_stL,mat_temp
TYPE(negf_env_obj) :: my_negf_env
INTEGER :: fcnt
CHARACTER(:), ALLOCATABLE :: SysName


CALL test_init(tol = 1.E-6_dp)

!choose system to test
SysName = "simple"

! -------------------- test GR ---------------------!
!initialization
open (1, file =  "test/" // SysName // "/GR/Tol.dat" )
READ (1, FMT=*) Tol
open (2, file = "test/" // SysName // "/GR/eta.dat" )
READ (2, FMT=*) eta
open (3, file = "test/" // SysName // "/GR/E.dat" )
READ (3, FMT=*) my_energy
close(1) 
close(2) 
close(3)

open (1, file = "test/" // SysName // "/quad/ab.dat" )
READ (1, FMT=*) my_a, my_b
close(1)


E   = CMPLX(my_energy,eta)
!end: initialization

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
CALL mat_read(my_negf_env%obj%H_L_onsite(1),"input/" // SysName // "/hL.dat")
CALL mat_read(my_negf_env%obj%H_L_onsite(2),"input/" // SysName // "/hR.dat")
CALL mat_read(my_negf_env%obj%S_L_onsite(1),"input/" // SysName // "/ShL.dat")
CALL mat_read(my_negf_env%obj%S_L_onsite(2),"input/" // SysName // "/ShR.dat")
!hopping leads
CALL mat_read(my_negf_env%obj%H_L_hopping(1),"input/" // SysName // "/tL.dat")
CALL mat_read(my_negf_env%obj%H_L_hopping(2),"input/" // SysName // "/tR.dat")
CALL mat_read(my_negf_env%obj%S_L_hopping(1),"input/" // SysName // "/StL.dat")
CALL mat_read(my_negf_env%obj%S_L_hopping(2),"input/" // SysName // "/StR.dat")
!leads -- system
CALL mat_read(my_negf_env%obj%H_LS(1),"input/" // SysName // "/TSL.dat")
CALL mat_read(my_negf_env%obj%H_LS(2),"input/" // SysName // "/TSR.dat")
CALL mat_read(my_negf_env%obj%S_LS(1),"input/" // SysName // "/STSL.dat")
CALL mat_read(my_negf_env%obj%S_LS(2),"input/" // SysName // "/STSR.dat")
!system
CALL mat_read(my_negf_env%obj%H_S,"input/" // SysName // "/HS.dat")
CALL mat_read(my_negf_env%obj%S_S,"input/" // SysName // "/SHS.dat")

 
CALL calc_GR(my_negf_env, my_energy, mat_GR)

WRITE(*,*) "GR = "
CALL mat_write(mat_GR)


!INTEGRAL -- quad testing
WRITE (*,*) "============================ check quad  ============================="
CALL quad(my_negf_env, my_a, my_b, Tol, my_integral,fcnt,calc_GR)
WRITE(*,*) "MY_INTEGRAL = "
CALL mat_write(my_integral)
WRITE(*,*) "fcnt = ", fcnt
CALL mat_read(matlab_integral,"test/" // SysName // "/GR/GR.dat")
CALL check_pass(mat_norm(matlab_integral),mat_norm(my_integral))

! finish tests
CALL test_finalise()

!detele - start
!CALL mat_read(mat_temp,"TSR.dat")
!WRITE(*,*) 'TSR.dat'
!CALL mat_write(mat_temp)
!delete - end



END PROGRAM test_quad

