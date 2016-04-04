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
REAL(KIND=dp)           :: my_energy,ReTr_integral_matlab, ImTr_integral_matlab
COMPLEX(KIND=dp)        :: E,Tr_my_integral,Tr_integral_matlab
REAL(KIND=dp)           :: my_a, my_b

TYPE(mat_z_obj)    :: mat_gL, mat_GR, my_integral,matlab_integral,matlab_GR,matlab_gs
TYPE(mat_d_obj)    :: mat_hL,mat_tL,mat_shL,mat_stL,mat_temp
TYPE(negf_env_obj) :: my_negf_env
INTEGER :: fcnt,fcnt_matlab
CHARACTER(:), ALLOCATABLE :: SysName


CALL test_init(tol = 1.E-7_dp)

!choose system to test
SysName = "simple"
!

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

!test surface GR sancho
WRITE (*,*) "======================= check surface_GR_sancho  ======================="
CALL mat_read(mat_hL, "input/" // SysName // "/hL.dat")
CALL mat_read(mat_tL, "input/" // SysName // "/tL.dat")
CALL mat_read(mat_shL, "input/" // SysName // "/ShL.dat")
CALL mat_read(mat_stL, "input/" // SysName // "/StL.dat")
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
CALL mat_read(matlab_gs,"test/" // SysName // "/gs/gs.dat")
CALL check_pass(mat_norm(matlab_gs),mat_norm(mat_gL))
!end: test surface GR sancho

!test calc_GR
WRITE (*,*) "=========================== check calc_GR  ==========================="

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

CALL mat_read(matlab_GR,"test/" // SysName // "/GR/GR.dat")

CALL check_pass(mat_norm(matlab_GR),mat_norm(mat_GR))
!end: test_GR

!test_quad
WRITE (*,*) "============================ check quad  ============================="
CALL quad(my_negf_env, my_a, my_b, Tol, my_integral,fcnt,calc_GR)
WRITE(*,*) "MY_INTEGRAL = "
CALL mat_write(my_integral)
WRITE(*,*) "fcnt = ", fcnt
CALL mat_read(matlab_integral,"test/" // SysName // "/GR/GR.dat")
!compare Tr of the integral matrix
Tr_my_integral = mat_trace(my_integral)
open (1, file =  "test/" // SysName // "/quad/Tr_my_integral.dat" )
READ (1, FMT=*) ReTr_integral_matlab,  ImTr_integral_matlab
close(1)
open (1, file =  "test/" // SysName // "/quad/fcnt.dat" )
READ (1, FMT=*) fcnt_matlab
close(1)
Tr_integral_matlab = CMPLX(ReTr_integral_matlab,ImTr_integral_matlab )
CALL check_pass(ABS(Tr_my_integral)+ABS(fcnt),ABS(Tr_integral_matlab)+ABS(fcnt))
!end: test_quad

! finish tests
CALL test_finalise()
END PROGRAM test_quad

