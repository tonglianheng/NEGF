PROGRAM test_Fermi
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

REAL(KIND=dp)           :: my_a, my_b, EF, T
COMPLEX(KIND = dp)      :: my_za, my_zb

TYPE(mat_z_obj)    :: mat_gL, mat_GR, my_integral,&
                      my_integral1, my_integral2, my_integral3,&
                      matlab_integral,matlab_gs,matlab_g11L, matlab_GR
!TYPE(mat_z_obj), DIMENSION(:)    :: mat_G_f_out

TYPE(mat_d_obj)    :: mat_hL,mat_tL,mat_shL,mat_stL,mat_temp,&
                      matlab_Img11L, matlab_Reg11L, matlab_ReGR, matlab_ImGR
TYPE(negf_env_obj) :: my_negf_env
INTEGER :: fcnt,fcnt1,fcnt2,fcnt3, iLead, ii
CHARACTER(:), ALLOCATABLE :: SysName
COMPLEX(KIND = dp),  DIMENSION(:), POINTER :: EE_matrix


END  PROGRAM test_Fermi

