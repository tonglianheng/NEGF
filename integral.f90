MODULE integral

USE kinds, ONLY: dp
USE matrix_types, ONLY: mat_d_obj, &
                        mat_z_obj, &
                        mat_norm, &
                        mat_copy, &
                        mat_create, &
                        mat_release, &
                        mat_nrows, &
                        mat_axpy, &
                        mat_real_to_complex, &
                        mat_scale,&
                        mat_imag,&
                        mat_real
                        
USE negf_env_types, ONLY: negf_env_obj,&
                          calc_g_less,&
                          calc_g_f1,&
                          calc_g_gamma2_g_df,&
                          calc_g 
USE mathconstants, ONLY: pi

#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  !public methods
  PUBLIC :: quad, sum_integral,&
            generate_X3,&
            linspace,&
            calc_Peq,&
            calc_Pneq,&
            calc_P,&
            calc_P_direct,&
            calc_P_comp
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'integral'
  INTEGER, PRIVATE, SAVE :: last_negf_env_id = 0

INTERFACE linspace
MODULE PROCEDURE linspace_d 
MODULE PROCEDURE linspace_z
END INTERFACE


CONTAINS

PURE FUNCTION eps(x) RESULT(res)
  REAL(KIND=dp), INTENT(IN) :: x
  REAL(KIND=dp) :: res
  INTEGER :: ee
  REAL(KIND=dp) :: abs_x
  ee = 0
  abs_x = ABS(x)
  IF (abs_x .LT. 1._dp) THEN
    DO WHILE (2._dp**ee .GT. abs_x)
      ee = ee - 1
    END DO
  ELSE
    DO WHILE (2._dp**ee .LE. abs_x)
      ee = ee + 1
    END DO
    ee = ee - 1
  END IF
  res = 2._dp**(ee - 52)
END FUNCTION eps


SUBROUTINE quad(my_negf_env, &
                a,b, &
                Tol, &
                mat_integral,&
                fcnt,&
                func)
  
  TYPE(negf_env_obj) :: my_negf_env
  TYPE(mat_z_obj), INTENT(OUT) :: mat_integral
  COMPLEX(KIND=dp), INTENT(IN)    :: a,b
  REAL(KIND=dp), INTENT(IN)  :: Tol
  INTEGER, INTENT(OUT),OPTIONAL :: fcnt

  INTERFACE
  SUBROUTINE func(my_negf_env, EE, output)
  USE matrix_types, ONLY: mat_z_obj
  USE negf_env_types, ONLY: negf_env_obj
  USE KINDS, ONLY: dp
  COMPLEX(KIND = dp), INTENT(IN):: EE  
  TYPE(negf_env_obj),INTENT(IN) :: my_negf_env
  TYPE(mat_z_obj),INTENT(INOUT) :: output
  END SUBROUTINE func
  END INTERFACE

 !technical variables
  INTEGER                        :: ii
  COMPLEX(KIND=dp)               :: h
  REAL(KIND=dp)                  :: hmin
  COMPLEX(KIND=dp),   DIMENSION(7)  :: X
  TYPE(mat_z_obj), DIMENSION(7)  :: mat_Y
  TYPE(mat_z_obj)                :: mat_Q1, mat_Q2, mat_Q3
  INTEGER                        :: warn, warn1, warn2, warn3
!default settings
 ! IF (.NOT.(PRESENT(Tol)))  Tol = 1.E-06_dp

  !Initialize with three unequal subintervals
  h = CMPLX(0.13579_dp,0.0_dp,KIND = dp)*(b-a)
  X = (/ a,&
         a+h,&
         a+(2._dp,0._dp)*h,&
         (a+b)/(2._dp,0._dp), &
         b-(2._dp,0._dp)*h,&
         b-h, &
         b /)
  DO ii=1,7
    CALL func(my_negf_env, X(ii), mat_Y(ii))
  END DO
  fcnt = 7;
  
  !Fudge endpoints to avoid infinities
  !IF isinfinite(mat_y(1)) THEN
  !  x1 = a+eps*(b-a);
  !  CALL calc_GR(my_negf_env,energy = x1,Gretarded = mat_Y(1))
  !  fcnt = fcnt+1;
  !END IF

  !Fudge endpoints to avoid infinities
  !IF isinfinite(mat_y(7)) THEN
  !  x1 = b-eps*(b-a);
  !  CALL calc_GR(my_negf_env,energy = x1,Gretarded = mat_Y(7))
  !  fcnt = fcnt+1;
  !END IF

  !Call the recursive core integrator
 
  ! this is how to call warnings
  !CALL cp_warn(__LOCATION__, "message dfafs a")



! Call the recursive core integrator.
hmin = eps(ABS(b-a))/1024._dp;
CALL quadstep(my_negf_env,X(1), X(3), mat_Y(1), mat_Y(2), mat_Y(3),&
              Tol, fcnt, hmin,  mat_Q1, warn1, func)
CALL quadstep(my_negf_env,X(3), X(5), mat_Y(3), mat_Y(4), mat_Y(5),&
              Tol, fcnt, hmin,  mat_Q2, warn2,func)
CALL quadstep(my_negf_env,X(5), X(7), mat_Y(5), mat_Y(6), mat_Y(7),&
              Tol, fcnt, hmin,  mat_Q3, warn3,func)


CALL mat_copy(mat_Q1,mat_integral)
CALL mat_axpy((1._dp,0._dp),'N',mat_Q2,mat_integral)
CALL mat_axpy((1._dp,0._dp),'N',mat_Q3,mat_integral)

warn = max(warn1,warn2,warn3)

SELECT CASE (warn)
  CASE (1)
CALL cp_warn(__LOCATION__, "Minimum step size reached; singularity possible")
  CASE (2)
CALL cp_warn(__LOCATION__, "Maximum function count exceeded. singularity likely")
END SELECT

    !cleanup
    CALL mat_release(mat_Q1)
    CALL mat_release(mat_Q2)
    CALL mat_release(mat_Q3)

    DO ii=1, 7
    CALL mat_release(mat_Y(ii))
    END DO
   ! DEALLOCATE (mat_Y)
END SUBROUTINE quad


RECURSIVE SUBROUTINE quadstep(my_negf_env, a,b, mat_fa, mat_fc, mat_fb, Tol, fcnt, hmin, mat_Q, warn,func)
!in-out
TYPE(negf_env_obj),INTENT(IN)   :: my_negf_env
COMPLEX(KIND=dp),INTENT(IN)     :: a,b
TYPE(mat_z_obj),INTENT(IN)      :: mat_fa, mat_fb, mat_fc
REAL(KIND=dp),INTENT(IN)        :: Tol
REAl(KIND=dp), INTENT(IN)       :: hmin
INTEGER,INTENT(INOUT)           :: fcnt        
INTEGER, INTENT(OUT)            :: warn
TYPE(mat_z_obj),INTENT(OUT)     :: mat_Q

  INTERFACE
  SUBROUTINE func(my_negf_env, EE, output)
  USE matrix_types, ONLY: mat_z_obj
  USE negf_env_types, ONLY: negf_env_obj
  USE KINDS, ONLY: dp
  COMPLEX(KIND = dp), INTENT(IN):: EE
  TYPE(negf_env_obj),INTENT(IN) :: my_negf_env
  TYPE(mat_z_obj),INTENT(INOUT) :: output
  END SUBROUTINE func
  END INTERFACE

!other
INTEGER,PARAMETER               :: maxfcnt=10000
COMPLEX(KIND=dp)                :: h,c
REAL(KIND=dp)                   :: NElements
COMPLEX(KIND=dp),DIMENSION(1:2) :: X
TYPE(mat_z_obj)                 :: mat_fd, mat_fe, mat_Q1,mat_Q2,&
                                   mat_Qac, mat_Qcb, mat_difQ
INTEGER                         :: Nrows,ii
INTEGER                         :: warnac,warncb
!Evaluate integrand twice in interior of subinterval [a,b]
h = b - a
c = (a + b)/(2._dp,0._dp)
X(1) = (a + c)/(2._dp,0._dp)
X(2) = (c + b)/(2._dp,0._dp)
CALL func(my_negf_env, X(1), mat_fd)
CALL func(my_negf_env, X(2), mat_fe)
fcnt = fcnt + 2


!Three point Simpson's rule:
! mat_Q1 = (h/6._dp)*(mat_fa + 4*mat_fc + mat_fb)
CALL mat_copy(mat_fa,mat_Q1)
CALL mat_axpy((4._dp,0._dp),'N',mat_fc,mat_Q1)
CALL mat_axpy((1._dp,0._dp),'N',mat_fb,mat_Q1)
CALL mat_scale(mat_Q1,h/(6._dp,0._dp))

!Five point double Simpson's rule:
! mat_Q2 = (h/12)*(fa + 4*fd + 2*fc + 4*fe + fb)
CALL mat_copy(mat_fa,mat_Q2)
CALL mat_axpy((4._dp,0._dp),'N',mat_fd,mat_Q2)
CALL mat_axpy((2._dp,0._dp),'N',mat_fc,mat_Q2)
CALL mat_axpy((4._dp,0._dp),'N',mat_fe,mat_Q2)
CALL mat_axpy((1._dp,0._dp),'N',mat_fb,mat_Q2)
CALL mat_scale(mat_Q2,h/(12._dp,0._dp))

!One step of Romberg extrapolation
!mat_Q = may_Q2 + (mat_Q2 - mat_Q1)/15
CALL mat_copy(mat_Q2,mat_Q)
CALL mat_axpy((-1._dp,0._dp),'N',mat_Q1,mat_Q)
CALL mat_scale(mat_Q,(1._dp,0._dp)/(15._dp,0._dp))
CALL mat_axpy((1._dp,0._dp),'N',mat_Q2,mat_Q)

!Check termination criterea

IF (fcnt .GT. maxfcnt) THEN
  !Maximum function count exceeded; singularity likely.
  warn = 2
 !cleanup
    CALL mat_release(mat_fd)
    CALL mat_release(mat_fe)
    CALL mat_release(mat_Q1)
    CALL mat_release(mat_Q2)
    CALL mat_release(mat_difQ)
 !end:cleanup
  RETURN
END IF

!check mat_nrows (mat_Q) for the optimality!!! Maybe should be passed as an input?
Nrows = mat_nrows(mat_Q)
NElements = REAL((Nrows), KIND=dp)
CALL mat_copy(mat_Q,mat_difQ)
CALL mat_axpy((-1._dp,0._dp),'N',mat_Q2,mat_difQ)
IF (ABS((mat_norm(mat_difQ))) .LE. Tol) THEN
  !Accuracy over this subinterval is acceptable
  warn = 0
 !cleanup
    CALL mat_release(mat_fd)
    CALL mat_release(mat_fe)
    CALL mat_release(mat_Q1)
    CALL mat_release(mat_Q2)
    CALL mat_release(mat_difQ)
 !end:cleanup
  RETURN
END IF

IF ((ABS(ABS(h)).LT.hmin).OR.(c.EQ.a).OR.(c.EQ.b)) THEN
  !Minimum step size reached; singularity possible.
  warn = 1
 !cleanup
    CALL mat_release(mat_fd)
    CALL mat_release(mat_fe)
    CALL mat_release(mat_Q1)
    CALL mat_release(mat_Q2)
    CALL mat_release(mat_difQ)
 !end:cleanup
  RETURN
END IF

!Subdivide into two subintervals
CALL quadstep(my_negf_env, a, c, mat_fa, mat_fd, mat_fc,&
              Tol, fcnt, hmin,  mat_Qac, warnac,func)
CALL quadstep(my_negf_env, c, b, mat_fc, mat_fe, mat_fb,&
              Tol, fcnt, hmin,  mat_Qcb, warncb,func)

CALL mat_copy(mat_Qac,mat_Q)
CALL mat_axpy((1._dp,0._dp),'N',mat_Qcb,mat_Q)
warn = MAX(warnac,warncb)

    !cleanup
    CALL mat_release(mat_fd)
    CALL mat_release(mat_fe)
    CALL mat_release(mat_Q1)
    CALL mat_release(mat_Q2)
    CALL mat_release(mat_Qac)
    CALL mat_release(mat_Qcb)
    CALL mat_release(mat_difQ)

END SUBROUTINE quadstep

SUBROUTINE sum_integral(my_negf_env, X, Tol, integral_sum, sum_fcnt,func)
!calculate a set of integrals given by points in X(i)
! to be paralelized!!!
TYPE(negf_env_obj), INTENT(IN) :: my_negf_env
COMPLEX(KIND = dp), DIMENSION(:), INTENT(IN) :: X
TYPE(mat_z_obj), INTENT(OUT) :: integral_sum
INTEGER, INTENT(OUT) :: sum_fcnt
REAL(KIND = dp), INTENT(IN) :: Tol
INTEGER :: fcnt
INTEGER :: N_int, ii
TYPE(mat_z_obj) :: temp_integral
  INTERFACE
  SUBROUTINE func(my_negf_env, EE, output)
  USE matrix_types, ONLY: mat_z_obj
  USE negf_env_types, ONLY: negf_env_obj
  USE KINDS, ONLY: dp
  COMPLEX(KIND = dp), INTENT(IN):: EE
  TYPE(negf_env_obj),INTENT(IN) :: my_negf_env
  TYPE(mat_z_obj),INTENT(INOUT) :: output
  END SUBROUTINE func
  END INTERFACE


N_int = size(X) - 1
CPASSERT (N_int .GE. 1)
sum_fcnt = 0
CALL quad(my_negf_env, X(1), X(2),Tol, integral_sum,sum_fcnt,func)


IF (N_int .GE. 2) THEN
  DO ii = 2, N_int
    CALL quad(my_negf_env,X(ii), X(ii+1),Tol, temp_integral,fcnt,func)
    CALL mat_axpy((1._dp,0._dp),'N' , temp_integral, integral_sum)   
    sum_fcnt = fcnt + sum_fcnt
  END DO
END IF

CALL mat_release(temp_integral)

END SUBROUTINE sum_integral

SUBROUTINE generate_X3 (my_negf_env, N1, N2, N3, X1, X2, X3)
!> generates three sets of points for sum_integral
!> X1 (real axis)
!> X2 (L) 
!> X3 (C)
TYPE(negf_env_obj), INTENT(IN) :: my_negf_env
COMPLEX(KIND = dp), DIMENSION(:), INTENT(OUT) :: X1, X2, X3
INTEGER,  INTENT(IN) :: N1, N2, N3
REAL(KIND = dp) :: R, x_tmp, x1_tmp, y1_tmp, phi1, phi2, mu_L, mu_R
REAL(KIND = dp), DIMENSION(N3) :: PHI
COMPLEX(KIND = dp) :: Z1, Z2, Z3, ZRe1, ZRe2, xC
INTEGER :: ii

mu_L = my_negf_env%obj%u_L(1)
mu_R = my_negf_env%obj%u_L(2)

!1. Horisontal sub-contour L
R =  (mu_L - my_negf_env%obj%Elow_int)&
 / 2._dp*(my_negf_env%obj%alpha_R_int)
!center of C contour
xC = mu_L - R - my_negf_env%obj%gamma_int
x_tmp = sqrt(R**2 - (my_negf_env%obj%delta_int)**2)


![Z2, Z3] is L contour
x1_tmp = xC + x_tmp
y1_tmp = my_negf_env%obj%delta_int
Z2 = CMPLX(x1_tmp, y1_tmp, KIND = dp)

x1_tmp = mu_L + my_negf_env%obj%Eup_int
!y1_tmp = my_negf_obj%obj%delta_int
Z3 = CMPLX(x1_tmp, y1_tmp, KIND = dp)

! Z1 -circle-> Z2 is C-contour
x1_tmp = xC - R
y1_tmp = my_negf_env%obj%eps_E
Z1 = CMPLX(x1_tmp, y1_tmp, KIND = dp)

! "real" interval is ZRe1 --> ZRe2
x1_tmp = min(mu_L, mu_R) - my_negf_env%obj%Eup_int
y1_tmp =  my_negf_env%obj%eps_E
ZRe1 = CMPLX(x1_tmp, y1_tmp, KIND = dp)

x1_tmp = max(mu_L, mu_R) + my_negf_env%obj%Eup_int 
ZRe2 = CMPLX(x1_tmp, y1_tmp, KIND = dp)

!borders of intervals are ready. get X1, X2, X3
!X1
CALL linspace(ZRe1, ZRe2, N1, X1)
!X2 
CALL linspace(Z2, Z3, N2, X2)
!X3 is a bit more complicated

phi1 = ASIN(my_negf_env%obj%eps_E/R)
phi2 = ASIN(my_negf_env%obj%delta_int/R) 

CALL linspace(pi-phi1, phi2, N3, PHI)

DO ii = 1, N3
  x1_tmp = cos(PHI(ii))*R + xC
  y1_tmp = sin(PHI(ii))*R
  X3(ii) = CMPLX(x1_tmp, y1_tmp, KIND = dp)
END DO

END SUBROUTINE generate_X3

SUBROUTINE linspace_z(d1,d2,n,grid)

INTEGER, INTENT(IN) :: n
COMPLEX(KIND = dp), INTENT(IN) :: d1, d2
COMPLEX(KIND = dp), DIMENSION(n), INTENT(OUT) :: grid

COMPLEX (KIND = dp) :: tmp1, tmp2
REAL(KIND = dp) :: Re1, Re2
INTEGER :: indxi

grid(1) = d1
    Re1 = DBLE(n-1)
    tmp1 = CMPLX( Re1 ,0._dp , KIND = dp)
DO indxi= 1,n-2
    Re2 = DBLE(indxi)
    tmp2 = CMPLX( Re2 ,0._dp , KIND = dp)
   grid(indxi+1) = d1+(tmp2*(d2-d1))/tmp1
END DO
grid(n) = d2

END SUBROUTINE linspace_z


SUBROUTINE linspace_d(d1,d2,n,grid)

INTEGER, INTENT(IN) :: n
REAL(KIND = dp), INTENT(IN) :: d1, d2
REAL(KIND = dp), DIMENSION(n), INTENT(OUT) :: grid

REAL (KIND = dp) :: tmp1, tmp2
INTEGER :: indxi

grid(1) = d1
    tmp1 = CMPLX( DBLE(n-1) ,0._dp , KIND = dp)
DO indxi= 0,n-2
    tmp2 = CMPLX( DBLE(indxi) ,0._dp , KIND = dp)
   grid(indxi+1) = d1+(tmp2*(d2-d1))/tmp1
END DO
grid(n) = d2

END SUBROUTINE linspace_d

!only for test purpose
SUBROUTINE calc_Peq(my_negf_env, Tol, N1, N2, N3, fcnt2, fcnt3, Peq)
!> equilibrium part of the density matrix
TYPE(negf_env_obj), INTENT(IN) :: my_negf_env
INTEGER, INTENT(IN) :: N1, N2, N3
INTEGER, INTENT(OUT) ::  fcnt2, fcnt3
TYPE(mat_d_obj), INTENT(OUT) :: Peq
REAL(KIND = dp), INTENT(IN) :: Tol

COMPLEX(KIND = dp), DIMENSION(N1) :: X1
COMPLEX(KIND = dp), DIMENSION(N2) :: X2
COMPLEX(KIND = dp), DIMENSION(N3) :: X3

TYPE(mat_z_obj) :: Int2, Int3, sum_res_G_f1

CALL generate_X3(my_negf_env, N1, N2, N3, X1, X2, X3)

!> Int1, Int2, Int 3 --> Real axis, L, C - contours 
!CALL sum_integral(my_negf_env, X1, Tol/N1, Int1, fcnt1,&
!                  calc_G_Gamma2_G_df)

CALL sum_integral(my_negf_env, X2, Tol/N2, Int2, fcnt2,&
                  calc_G_f1)

CALL sum_integral(my_negf_env, X3, Tol/N3, Int3, fcnt3,&               
                  calc_G_f1)


CALL calc_sum_res_G_f(my_negf_env, 1, sum_res_G_f1)

!Peq = -1/pi*mat_imag(Int2 + Int3 - sum_res_G_f1 )
CALL mat_axpy((1._dp, 0._dp), 'N', Int2, Int3)
CALL mat_axpy((1._dp, 0._dp), 'N', Int3, sum_res_G_f1)
CALL mat_imag(sum_res_G_f1, Peq)
CALL mat_scale(Peq, -1._dp/pi)

CALL mat_release(Int2)
CALL mat_release(Int3)
CALL mat_release(sum_res_G_f1)

END SUBROUTINE calc_Peq

!only for test purpose
SUBROUTINE calc_Pneq(my_negf_env, Tol, N1, N2, N3, fcnt1, Pneq)
!> equilibrium part of the density matrix
TYPE(negf_env_obj), INTENT(IN) :: my_negf_env
INTEGER, INTENT(IN) :: N1, N2, N3
INTEGER, INTENT(OUT) ::  fcnt1
TYPE(mat_d_obj), INTENT(OUT) :: Pneq
REAL(KIND = dp), INTENT(IN) :: Tol


COMPLEX(KIND = dp), DIMENSION(N1) :: X1
COMPLEX(KIND = dp), DIMENSION(N2) :: X2
COMPLEX(KIND = dp), DIMENSION(N3) :: X3
TYPE(mat_z_obj) :: Int1

CALL generate_X3(my_negf_env, N1, N2, N3, X1, X2, X3)

!> Int1, Int2, Int 3 --> Real axis, L, C - contours
CALL sum_integral(my_negf_env, X1, Tol/N1, Int1, fcnt1,&
                  calc_G_Gamma2_G_df)

!CALL sum_integral(my_negf_env, X2, Tol/N2, Int2, fcnt2,&
!                  calc_G_f1)

!CALL sum_integral(my_negf_env, X3, Tol/N3, Int3, fcnt3,&
!                  calc_G_f1)

!Pneq = 1/(2*pi)*mat_real(Int1)
CALL mat_real(Int1, Pneq)
CALL mat_scale(Pneq, 1._dp/(2._dp*pi))

CALL mat_release(Int1)

END SUBROUTINE calc_Pneq


SUBROUTINE calc_P(my_negf_env, Tol, N1, N2, N3, fcnt1, fcnt2, fcnt3, P)
!> density matrix full (Peq + Pneq)
TYPE(negf_env_obj), INTENT(IN) :: my_negf_env
INTEGER, INTENT(IN) :: N1, N2, N3
INTEGER, INTENT(OUT) :: fcnt1,  fcnt2, fcnt3
TYPE(mat_d_obj), INTENT(OUT) :: P
REAL(KIND = dp), INTENT(IN) :: Tol

COMPLEX(KIND = dp), DIMENSION(N1) :: X1
COMPLEX(KIND = dp), DIMENSION(N2) :: X2
COMPLEX(KIND = dp), DIMENSION(N3) :: X3

TYPE(mat_d_obj) :: Pneq
TYPE(mat_z_obj) :: Int1, Int2, Int3, sum_res_G_f1

CALL generate_X3(my_negf_env, N1, N2, N3, X1, X2, X3)

!> Int1, Int2, Int 3 --> Real axis, L, C - contours 
CALL sum_integral(my_negf_env, X1, Tol/N1, Int1, fcnt1,&
                  calc_G_Gamma2_G_df)

CALL sum_integral(my_negf_env, X2, Tol/N2, Int2, fcnt2,&
                  calc_G_f1)

CALL sum_integral(my_negf_env, X3, Tol/N3, Int3, fcnt3,&               
                  calc_G_f1)


CALL calc_sum_res_G_f(my_negf_env, 1, sum_res_G_f1)

!>Peq = -1/pi*mat_imag(Int2 + Int3 - sum_res_G_f1 )
CALL mat_axpy((1._dp, 0._dp), 'N', Int2, Int3)
CALL mat_axpy((1._dp, 0._dp), 'N', Int3, sum_res_G_f1)
CALL mat_imag(sum_res_G_f1, P)
CALL mat_scale(P, -1._dp/pi)

!>Pneq = 1/(2*pi) * Int1
CALL mat_real(Int1, Pneq)
CALL mat_scale(Pneq, 1._dp/(2._dp*pi))
CALL mat_axpy(1._dp, 'N', Pneq, P)

!>clean up
CALL mat_release(Int1)
CALL mat_release(Int2)
CALL mat_release(Int3)
CALL mat_release(sum_res_G_f1)
CALL mat_release(Pneq)

END SUBROUTINE calc_P

!used for test purposes (direct calc of density matrix)
SUBROUTINE calc_P_direct(my_negf_env, Tol, fcnt, P_direct)
TYPE(negf_env_obj), INTENT(IN) :: my_negf_env
TYPE(mat_d_obj), INTENT(OUT) :: P_direct
REAL(KIND = dp), INTENT(IN) :: Tol
INTEGER, INTENT(OUT) :: fcnt

REAL(KIND = dp) :: Re_a, Re_b, mu_L, mu_R, eta
COMPLEX(KIND = dp) :: a, b, tmp
TYPE(mat_z_obj) :: P_cmplx

mu_L = my_negf_env%obj%u_L(1)
mu_R = my_negf_env%obj%u_L(2)
eta =  my_negf_env%obj%eps_E

Re_a = my_negf_env%obj%Elow_int * my_negf_env%obj%alpha_R_int + min(mu_L,mu_R)
Re_b = my_negf_env%obj%Eup_int +  max(mu_L, mu_R)
a = CMPLX(Re_a, eta, KIND = dp)
b = CMPLX(Re_b, eta, KIND = dp)
CALL quad(my_negf_env, a, b, Tol, P_cmplx, fcnt, calc_G_less)

tmp =  CMPLX(0.0_dp, -0.5_dp/pi, KIND = dp)
CALL mat_scale(P_cmplx, tmp)
CALL mat_real(P_cmplx, P_direct)

CALL mat_release(P_cmplx)

END SUBROUTINE calc_P_direct

!used for test purposes (direct calc of density matrix)
SUBROUTINE calc_P_comp(my_negf_env, Tol, fcnt, fcnt1, P_direct, P_indirect)
TYPE(negf_env_obj), INTENT(IN) :: my_negf_env
TYPE(mat_d_obj), INTENT(OUT) :: P_direct, P_indirect
REAL(KIND = dp), INTENT(IN) :: Tol
INTEGER, INTENT(OUT) :: fcnt, fcnt1

REAL(KIND = dp) :: Re_a, Re_b, mu_L, mu_R, eta
COMPLEX(KIND = dp) :: a, b, tmp, tmp1
TYPE(mat_z_obj) :: P_cmplx, P_cmplx1

mu_L = my_negf_env%obj%u_L(1)
mu_R = my_negf_env%obj%u_L(2)
eta =  my_negf_env%obj%eps_E

Re_a = my_negf_env%obj%Elow_int * my_negf_env%obj%alpha_R_int + min(mu_L,mu_R)

Re_b = my_negf_env%obj%Eup_int +  max(mu_L, mu_R)
a = CMPLX(Re_a, eta, KIND = dp)
b = CMPLX(Re_b, eta, KIND = dp)
CALL quad(my_negf_env, a, b, Tol, P_cmplx, fcnt, calc_G_less)

CALL quad(my_negf_env, a, b, Tol, P_cmplx1, fcnt1, calc_G_f1)


tmp =  CMPLX(0.0_dp, -0.5_dp/pi, KIND = dp)
CALL mat_scale(P_cmplx, tmp)
CALL mat_real(P_cmplx, P_direct)

tmp1 = CMPLX(-1._dp/pi, 0.0_dp, KIND = dp)
CALL mat_scale(P_cmplx1, tmp1)
CALL mat_imag(P_cmplx1, P_indirect)

CALL mat_release(P_cmplx)
CALL mat_release(P_cmplx1)

END SUBROUTINE calc_P_comp

SUBROUTINE calc_sum_res_G_f(my_negf_env, iLead, sum_res_G_f)
!>sum residues of the f(iLead) (or f(iLead))
TYPE(negf_env_obj), INTENT(IN) :: my_negf_env
INTEGER, INTENT(IN) :: iLead
TYPE(mat_z_obj), INTENT(OUT) :: sum_res_G_f

TYPE(mat_z_obj) :: G
COMPLEX(KIND = dp) ::  tmp_cmplx, EE_tmp
INTEGER :: nrows, niu_max, ii
REAL(KIND = dp) :: Temp

Temp = my_negf_env%obj%Temperature

niu_max = floor(my_negf_env%obj%delta_int / (2*pi*Temp) - 0.5_dp) + 1

!zero matrix
nrows = mat_nrows(my_negf_env%obj%H_S)
CALL mat_create(sum_res_G_f, nrows, nrows)


EE_tmp = CMPLX(my_negf_env%obj%u_L(iLead), 0._dp, KIND = dp)

ii=1

DO WHILE (ii .LE.  niu_max)
  tmp_cmplx = &
  CMPLX(0._dp, pi*Temp*(2._dp*ii-1._dp), KIND = dp)
  CALL calc_G(my_negf_env, EE_tmp + tmp_cmplx, G)
  !sum_res_G_f = sum_res_G_f + G
  CALL mat_axpy((1._dp,0._dp), 'N', G, sum_res_G_f)
  ii = ii + 1
END DO

CALL mat_scale(sum_res_G_f, CMPLX(0._dp, -2*pi*Temp, KIND = dp))

CALL mat_release(G)

END SUBROUTINE calc_sum_res_G_f


END MODULE integral
