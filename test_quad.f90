PROGRAM test_quad
USE test_utils, ONLY: check_pass, &
                       test_init,&
                       test_finalise
USE kinds, ONLY: dp
USE matrix_types
USE sancho_method
USE negf_env_types
USE integral
USE mathconstants, ONLY: pi

IMPLICIT NONE

REAL(KIND=dp)           :: Tol
REAL(KIND=dp)           :: eta
REAL(KIND=dp)           :: my_energy,ReTr_integral_matlab, &
                           NStates_matlab, NStates,&
                           ReE, ImE
COMPLEX(KIND=dp)        :: E,Tr_my_integral,Tr_integral

REAL(KIND=dp)           :: my_a, my_b, my_c, EF, T
COMPLEX(KIND = dp)      :: my_za, my_zb, a, b

TYPE(mat_z_obj)    :: mat_gL, mat_GR, my_integral,&
                      my_integral1, my_integral2, my_integral3,&
                      matlab_integral,matlab_gs,matlab_g11L, &
                      matlab_GR, G_less, G_f, P_direct_z,&
                      P_G_less, P_G_f
!TYPE(mat_z_obj), DIMENSION(:)    :: mat_G_f_out

TYPE(mat_d_obj)    :: mat_hL,mat_tL,mat_shL,mat_stL,mat_temp,&
                      matlab_Img11L, matlab_Reg11L, matlab_ReGR,&
                      matlab_ImGR, Peq, Pneq, P, P_direct,&
                      P_indirect,Im_G_f, mat_d_tmp1, mat_d_tmp2
TYPE(negf_env_obj) :: my_negf_env
INTEGER :: fcnt,fcnt1,fcnt2,fcnt3, sum_fcnt, iLead, ii, N1, N2, N3
CHARACTER(:), ALLOCATABLE :: SysName, FileName
CHARACTER(10), DIMENSION(3) :: FileNames
COMPLEX(KIND = dp),  DIMENSION(1000) :: EE_matrix, fF_matrix
COMPLEX (KIND = dp), DIMENSION(4) :: XX
COMPLEX(KIND = dp), DIMENSION(:), POINTER :: X1
COMPLEX(KIND = dp), DIMENSION(:), POINTER :: X2
COMPLEX(KIND = dp), DIMENSION(:), POINTER :: X3

REAL(KIND = dp), DIMENSION(1000) :: E_vector
COMPLEX(KIND = dp), DIMENSION(1000) :: Z_vector


CALL test_init(tol = 1e-4_dp)

!choose system to test. Avail: zzCNT2,4,7 and simple
SysName = "zzCNT2"
!

!initialization
open (1, file =  "test/" // SysName // "/general/Tol.dat" )
READ (1, FMT=*) Tol
open (2, file = "test/" // SysName // "/general/eta.dat" )
READ (2, FMT=*) eta
open (3, file = "test/" // SysName // "/general/E.dat" )
READ (3, FMT=*) my_energy
close(1) 
close(2) 
close(3)

Tol = 1.e-6_dp
eta = 1.e-4_dp
E   = CMPLX(my_energy,eta,KIND = dp)

EF = 0.0_dp
T  = 0.027_dp
iLead = 2

!end: initialization

!=========================test surface GR sancho====================!
WRITE (*,*) "========= check surface_GR_sancho  (LEFT)=============="
CALL mat_read(mat_hL, "input/" // SysName // "/hL.dat")
CALL mat_read(mat_tL, "input/" // SysName // "/tL.dat")
CALL mat_read(mat_shL, "input/" // SysName // "/ShL.dat")
CALL mat_read(mat_stL, "input/" // SysName // "/StL.dat")
CALL surface_GR_sancho(GR_surface = mat_gL, &
                        energy    = E, &
                        H_onsite  = mat_hL, &
                        H_hopping = mat_tL, &
                        S_onsite  = mat_shL, &
                        S_hopping = mat_stL, &
                        tolerance = Tol)

CALL mat_read(matlab_Img11L,"test/" // SysName // "/gs/Img11L.dat")
CALL mat_read(matlab_Reg11L,"test/" // SysName // "/gs/Reg11L.dat")

CALL mat_CMPLX(matlab_Reg11L, matlab_Img11L, matlab_g11L)

!WRITE (*,*) "gs="
!CALL mat_write(mat_gL)

CALL check_pass(mat_norm(matlab_g11L),mat_norm(mat_gL))
!=====================end: test surface GR sancho====================!



!========================== Create enviroment =======================!
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

!Bias u_L
   my_negf_env%obj%u_L(1) = 0.0_dp
   my_negf_env%obj%u_L(2) = 0.0_dp

   my_negf_env%obj%EF_L(1) = 0.0_dp
   my_negf_env%obj%EF_L(2) = 0.0_dp

    my_negf_env%obj%Temperature = 2.7E-2_dp  !it is 300 K in eV. Change!
    my_negf_env%obj%gamma_int   = 2._dp !it is in eV. Change!
    my_negf_env%obj%delta_int   = 2._dp !it is in eV. Change!
    my_negf_env%obj%Elow_int    = -10._dp !it is in eV. Change!
    my_negf_env%obj%Eup_int     = 10._dp !it is in eV. Change!
    my_negf_env%obj%alpha_R_int = 2._dp


!============================test calc_GR===========================!
WRITE (*,*) "=================== check calc_GR  ====================="
 
CALL calc_G(my_negf_env, CMPLX(my_energy,eta, KIND = dp), mat_GR)

!WRITE(*,*) "Tr GR = "
!CALL mat_trace(mat_GR)

CALL mat_read(matlab_ImGR,"test/" // SysName // "/GR/ImGR.dat")
CALL mat_read(matlab_ReGR,"test/" // SysName // "/GR/ReGR.dat")
CALL mat_CMPLX(matlab_ReGR, matlab_ImGR, matlab_GR)
CALL check_pass(mat_norm(matlab_GR),mat_norm(mat_GR))
!========================end: test_GR===============================!

!===================================test_quad_d======================
WRITE (*,*) "======================== check quad_d  ================="
!read integration limits
open (1, file = "test/" // SysName // "/quad/ab.dat" )
READ (1, FMT=*) my_a, my_b
close(1)

CALL quad(my_negf_env, CMPLX(my_a, eta, KIND = dp), &
                       CMPLX(my_b, eta, KIND = dp), &
                       Tol, my_integral,fcnt,calc_G)

!WRITE(*,*) "MY_INTEGRAL "
!CALL mat_write(my_integral)
WRITE(*,*) "fcnt = ", fcnt

!Find trace of the GR calculated in FORTRAN
Tr_my_integral = mat_trace(my_integral)
NStates = -1/pi*REAL(AIMAG(Tr_my_integral))
WRITE(*,*) " Nstates =  ", Nstates

!open NStates_matlab
open (1, file =  "test/" // SysName // "/quad/NStates.dat" )
READ (1, FMT=*) NStates_matlab
close(1)
WRITE(*,*) " Nstates_matlab =  ", Nstates_matlab

!compare NStates (theory (matlab)) with NStates (FORTRAN)
CALL check_pass(NStates,NStates_matlab)
!===============================end: test_quad_d====================!

!===================================test_quad_z=====================!
WRITE (*,*) "======================== check quad_z  ================="

CALL quad(my_negf_env,CMPLX(my_b, eta,     KIND = dp),&
                      CMPLX(my_b, 10.0_dp, KIND = dp),&
                      Tol, my_integral1,fcnt1,calc_G)

CALL quad(my_negf_env,CMPLX(my_b, 10._dp, KIND = dp),&
                      CMPLX(my_a, 10._dp, KIND = dp),&
                      Tol, my_integral2,fcnt2,calc_G)


CALL quad(my_negf_env,CMPLX(my_a, 10._dp, KIND = dp),&
                      CMPLX(my_a, eta,    KIND = dp), &
                      Tol, my_integral3,fcnt3,calc_G)

fcnt = fcnt1+fcnt2+fcnt3;

CALL mat_copy(my_integral1,my_integral)
CALL mat_axpy((1.0_dp, 0._dp),'N', my_integral2, my_integral);
CALL mat_axpy((1.0_dp, 0._dp), 'N', my_integral3, my_integral);
CALL mat_scale(my_integral, (-1.0_dp,0.0_dp) )

!WRITE(*,*) "MY_INTEGRAL "
!CALL mat_write(my_integral)
WRITE(*,*) "fcnt = ", fcnt

!Find trace of the GR calculated in FORTRAN
Tr_my_integral = mat_trace(my_integral)
NStates = -1/pi*REAL(AIMAG(Tr_my_integral))
WRITE(*,*) " Nstates =  ", Nstates

!open NStates_matlab
open (1, file =  "test/" // SysName // "/quad/NStates.dat" )
READ (1, FMT=*) NStates_matlab
close(1)
WRITE(*,*) " Nstates_matlab =  ", Nstates_matlab

!compare NStates (theory (matlab)) with NStates (FORTRAN)
CALL check_pass(NStates,NStates_matlab)
!===============================end: test_quad_z=====================!

!============================   calc_G_df ===========================!
WRITE (*,*) "======================== calc_G_f1  ===================="

!CALL quad(my_negf_env,CMPLX(my_a, eta,     KIND = dp),&
!                      CMPLX(my_a, 1.0_dp, KIND = dp),&
!                      Tol, my_integral1,fcnt1,calc_G_f2)

!CALL quad(my_negf_env,CMPLX(my_a, 1.0_dp, KIND = dp),&
!                      CMPLX(my_b, 1.0_dp, KIND = dp),&
!                      Tol, my_integral2,fcnt2,calc_G_f2)


!CALL quad(my_negf_env,CMPLX(my_b, 1.0_dp, KIND = dp),&
!                      CMPLX(my_b, eta,    KIND = dp), &
!                      Tol, my_integral3,fcnt3,calc_G_f2)

!fcnt = fcnt3+fcnt1+fcnt2


!CALL mat_copy(my_integral1,my_integral)
!CALL mat_axpy((1.0_dp, 0._dp),'N', my_integral2, my_integral);
!CALL mat_axpy((1.0_dp, 0._dp), 'N', my_integral3, my_integral);
!CALL mat_scale(my_integral, (-1.0_dp,0.0_dp) )

XX(1) = CMPLX(my_a, eta,     KIND = dp)
XX(2) = CMPLX(my_a, 1.0_dp, KIND = dp)
XX(3) = CMPLX(my_b, 1.0_dp, KIND = dp)
XX(4) = CMPLX(my_b, eta,    KIND = dp)
 
CALL sum_integral(my_negf_env, XX, Tol, my_integral, fcnt, calc_G_f1)


!WRITE(*,*) "MY_INTEGRAL "
!CALL mat_write(my_integral)
WRITE(*,*) "fcnt = ", fcnt

!Find trace of the GR calculated in FORTRAN
Tr_my_integral = mat_trace(my_integral)
NStates = -1/pi*REAL(AIMAG(Tr_my_integral))
WRITE(*,*) " Nstates =  ", Nstates

!open NStates_matlab
open (1, file =  "test/" // SysName // "/quad/NStates.dat" )
READ (1, FMT=*) NStates_matlab
close(1)
WRITE(*,*) " Nstates_matlab =  ", Nstates_matlab

!compare NStates (theory (matlab)) with NStates (FORTRAN)
CALL check_pass(NStates,4._dp)
!===============================end: calc_G_f =======================!

!================================deep test G_f=======================!
FileName = "out/FermiFunction.dat"
OPEN(20, file = FileName)

DO ii= 1, 1000
  ReE = 0.1_dp*ii
  ImE = 1.0_dp
  EE_matrix(ii) = CMPLX(ReE, ImE, KIND=dp)
  CALL Fermi(my_negf_env, EE_matrix(ii), iLead, fF_matrix(ii))
  WRITE (20, FMT = *)  EE_matrix(ii), fF_matrix(ii)
END DO 

CLOSE(20)
!=============================end: deep test G_f=====================!


!======================== calc_G_GammaR_G_df=========================!
WRITE (*,*) "===================== calc_G_GammaR_G_df  =============="

!Fermi function u_L
!   my_negf_env%obj%u_L(1) = -1.0_dp
!   my_negf_env%obj%u_L(2) = 1.0_dp

CALL quad(my_negf_env, CMPLX(-10, eta,  KIND = dp),&
                       CMPLX(10, eta,  KIND = dp),&
                       Tol, my_integral, fcnt, calc_G_Gamma2_G_df)


!WRITE(*,*) "MY_INTEGRAL "
!CALL mat_write(my_integral)
!WRITE(*,*) "fcnt = ", fcnt

!Find trace of the GR calculated in FORTRAN
Tr_my_integral = mat_trace(my_integral)
NStates = -1/pi*REAL(REAL(Tr_my_integral))
WRITE(*,*) " Nstates =  ", Nstates

!open NStates_matlab
open (1, file =  "test/" // SysName // "/quad/NStates.dat" )
READ (1, FMT=*) NStates_matlab
close(1)
WRITE(*,*) " Nstates_matlab =  ", Nstates_matlab

!compare NStates (theory (matlab)) with NStates (FORTRAN)
CALL check_pass(NStates, 0._dp)
!====================== end: calc_G_Gamma_R_G_df ====================!


!======================== check sum integrals =======================!
WRITE (*,*) "===================== sum integrals====================="

!my_negf_env%obj%u_L(1) = 0._dp
!my_negf_env%obj%u_L(2) = 0._dp

!my_negf_env%obj%Elow_int = -100._dp
!my_negf_env%obj%Eup_int = 100._dp


N1 = 10
N2 = 10
N3 = 10

ALLOCATE(X1(N1))
ALLOCATE(X2(N2))
ALLOCATE(X3(N3))
CALL generate_X3(my_negf_env, N1, N2, N3, X1, X2, X3)


FileNames(1) = "out/X1.dat"
FileNames(2) = "out/X2.dat"
FileNames(3) = "out/X3.dat"

!WRITE(*,*) FileNames

OPEN(10, file =   FileNames(1))
OPEN(20, file =   FileNames(2))
OPEN(30, file =   FileNames(3))
 DO ii = 1, N1
    WRITE(10, FMT="(F15.7,  F15.7)")  X1(ii)
    WRITE(20, FMT="(F15.7,  F15.7)")  X2(ii)
    WRITE(30, FMT="(F15.7,  F15.7)")  X3(ii)
 END DO
CLOSE(10)
CLOSE(20)
CLOSE(30)

!my_negf_env%obj%u_L(1) = 0._dp
!my_negf_env%obj%u_L(2) = 1._dp
!my_negf_env%obj%eps_E =1.E-8_dp
!my_negf_env%obj%eps_Sancho =1.E-8_dp
!my_negf_env%obj%Temperature =1.E-2_dp

!Tol = 1e-8


CALL sum_integral(my_negf_env, X1, Tol/N1,& 
                  my_integral1, sum_fcnt, calc_G_Gamma2_G_df)

CALL quad(my_negf_env, CMPLX(-10, my_negf_env%obj%eps_E,  KIND = dp),&
                       CMPLX(10,  my_negf_env%obj%eps_E,  KIND = dp),&
                       Tol, my_integral2, fcnt, calc_G_Gamma2_G_df)

WRITE(*,*) "calc_G_Gamma2_G_df - direct calc = ", mat_norm(my_integral2)
WRITE(*,*) "fcnt= ", fcnt
WRITE(*,*) "sum_fcnt= ", sum_fcnt

CALL check_pass(mat_norm(my_integral1),mat_norm(my_integral2))
!================================= Peq ==============================!
WRITE (*,*) "======================= check Peq======================="
CALL calc_Peq(my_negf_env, Tol, N1, N2, N3, fcnt2, fcnt3, Peq)

WRITE(*,*) "mat_trace (Peq) = ", mat_trace(Peq)
WRITE(*,*) "fcnt2 = ", fcnt2
WRITE(*,*) "fcnt3 = ", fcnt3
!========================== end: check Peq  =========================!

!================================= P ================================!
WRITE (*,*) "======================= check P ========================"

my_negf_env%obj%u_L(1) = 0._dp
my_negf_env%obj%u_L(2) = 0._dp

! contour integration
!my_negf_env%obj%u_L(1) = 0.1
CALL calc_P(my_negf_env, Tol, N1, N2, N3, fcnt1, fcnt2, fcnt3, P)
WRITE(*,*) "CONTOUR INTEGRATION"
WRITE(*,*) "mat_trace (P) = ", mat_trace(P)
WRITE(*,*) "fcnt1 = ", fcnt1
WRITE(*,*) "fcnt2 = ", fcnt2
WRITE(*,*) "fcnt3 = ", fcnt3


!Tol = 1.0E-12
!>direct integration
CALL calc_P_direct(my_negf_env, Tol, fcnt, P_direct)
WRITE(*,*) "DIRECT INTEGRATION"
WRITE(*,*) "mat_trace (P) = ", mat_trace(P_direct)
WRITE(*,*) "fcnt = ", fcnt


!============================ end: check P  =========================!

!E = E - 3._dp

CALL linspace(-10._dp,10._dp, 1000, E_vector)

OPEN(10, file = "out/TrGless.dat")
OPEN(20, file =  "out/TrImGf.dat")

DO ii = 1, 1000 

E = CMPLX(E_vector(ii),my_negf_env%obj%eps_E, KIND = dp)

CALL calc_G_less(my_negf_env, E, G_less)
!WRITE(*,*) "Tr G_less (E) = ", mat_trace(G_less)
WRITE(10, FMT="(F15.7,  F15.7)"), E_vector(ii),&
                                  aimag(mat_trace(G_less))

CALL calc_GR_f1(my_negf_env, E, G_f)
!CALL calc_P_direct(my_negf_env, Tol, fcnt, P)
a = mat_trace(G_f)
a = CMPLX(0.0_dp, -2*real(aimag(a)))
!WRITE(*,*) "Tr - 2*i* Im G_f (E) = ", a
WRITE(20, FMT="(F15.7,  F15.7)"), E_vector(ii), &
                                  aimag(a)

END DO

CLOSE(10)
CLOSE(20)

!===========================Compare P========================!


CALL calc_P_comp(my_negf_env, Tol, fcnt, fcnt1,&
                 P_direct, P_indirect)
WRITE(*,*) "P_direct",   mat_trace(P_direct)
WRITE(*,*) "P_indirect", mat_trace(P_indirect)


!========================= end: check P-2  =========================!
!CALL linspace(-10._dp,10._dp, 1000, E_vector)
!
!OPEN(10, file = "out/TrGless.dat")
!OPEN(20, file =  "out/TrImGf.dat")
!OPEN(30, file = "out/Gless-ImGf.dat")
!
!my_a = 0.0_dp
!my_b = 0.0_dp
!my_c = 0.0_dp
!
!DO ii = 1, 1000
!
!E = CMPLX(E_vector(ii),my_negf_env%obj%eps_E, KIND = dp)
!
!CALL calc_P_GR_f1(my_negf_env, E, P_G_f)
!CALL mat_imag(P_G_f, mat_d_tmp2)
!WRITE(20, FMT="(F15.7,  F15.7)"), E_vector(ii), &
!                                  mat_trace(mat_d_tmp2)
!
!my_c = my_c + mat_trace(mat_d_tmp2)*0.02_dp
!
!CALL calc_P_G_less(my_negf_env, E, P_G_less)
!CALL mat_real(P_G_less, mat_d_tmp1)
!WRITE(10, FMT="(F15.7,  F15.7)"), E_vector(ii),&
!                                  mat_trace(mat_d_tmp1)
!
!my_b = my_b + mat_trace(mat_d_tmp1)*0.02_dp
!
!CALL mat_axpy(-1.0_dp, 'N', mat_d_tmp2, mat_d_tmp1)
!WRITE (30, FMT='(F15.7, F15.7)'), E_vector(ii), &
!                                  mat_norm(mat_d_tmp1)
!
!my_a = my_a + mat_norm(mat_d_tmp1) * 0.02_dp
!
!END DO
!
!PRINT *, "error in integral expected to be = ", my_a
!PRINT *, "integral of trace Gless = ", my_b
!PRINT *, "integral of trace GR_f1 = ", my_c
!
!CLOSE(10)
!CLOSE(20)
!CLOSE(30)

WRITE(*,*) "!================= FINAL TEST P vs P ================!"

CALL &
quad(my_negf_env, CMPLX(-10, my_negf_env%obj%eps_E,  KIND = dp),&
                  CMPLX(10,  my_negf_env%obj%eps_E,  KIND = dp),&
                  Tol, my_integral2, fcnt2, calc_P_GR_f1)

CALL &
quad(my_negf_env, CMPLX(-10, my_negf_env%obj%eps_E,  KIND = dp),&
                  CMPLX(10,  my_negf_env%obj%eps_E,  KIND = dp),&
                  Tol, my_integral1, fcnt1, calc_P_G_less)

WRITE(*,*) "Tr P_G_less", mat_trace(my_integral1)
WRITE(*,*) "fcnt 1 = ", fcnt1
WRITE(*,*) "Tr P_GR_f1", mat_trace(my_integral2)
WRITE(*,*) "fcnt 2 = ", fcnt2


CALL test_finalise()
END PROGRAM test_quad

