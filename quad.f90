MODULE integral

USE kinds, ONLY: dp
USE matrix_types!, ONLY: mat_d_obj, &
                !       mat_z_obj, &
                !        mat_inv_lu, &
                !        mat_norm, &
                !        mat_copy, &
                !        mat_create, &
                !        mat_release, &
                !        mat_nrows, &
                !        mat_ncols, &
                !        mat_axpy, &
                !        mat_mult, &
                !        mat_real_to_complex
USE negf_env_types

#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  !public methods
  PUBLIC :: quad
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'integral'
  INTEGER, PRIVATE, SAVE :: last_negf_env_id = 0

 
CONTAINS

SUBROUTINE quad(my_negf_env, &
                a,b, &
                Tol, &
                mat_integral)
  
  TYPE(negf_env_obj) :: my_negf_env
  TYPE(mat_z_obj), INTENT(OUT) :: mat_integral
  REAL(KIND=dp), INTENT(IN)    :: a,b
  REAL(KIND=dp), INTENT(IN), OPTIONAL    :: Tol
  !technical variables
  INTEGER                        :: fcnt,ii
  REAL(KIND=dp)                  :: h,hmin
  REAL(KIND=dp),   DIMENSION(7)  :: X
  TYPE(mat_z_obj), DIMENSION(7)  :: mat_Y
  REAL(KIND=dp)                  :: x1
  TYPE(mat_z_obj)                :: mat_Q1, mat_Q2, mat_Q3
  INTEGER                        :: warn1, warn2, warn3
  !default settings
  IF (.NOT.(PRESENT(Tol)))  Tol = 1.E-06_dp

  !Initialize with three unequal subintervals
  h = 0.13579_dp*(b-a)
  X = (/ a,&
         a+h,&
         a+2._dp*h,&
        (a+b)/2._dp, &
         b-2._dp*h,&
         b-h, &
         b /)
  DO ii=1,7
    CALL calc_GR(my_negf_env,energy = X(ii),Gretarded = mat_Y(ii))
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
hmin = eps(b-a)/1024._dp;
CALL quadstep(my_negf_env,X(1), X(3), mat_Y(1), mat_Y(2), mat_Y(3),&
              Tol, fcnt, hmin,  mat_Q1, warn1)
CALL quadstep(my_negf_env,X(3), X(5), mat_Y(3), mat_Y(4), mat_Y(5),&
              Tol, fcnt, hmin,  mat_Q2, warn2)
CALL quadstep(my_negf_env,X(5), X(7), mat_Y(5), mat_Y(6), mat_Y(7),&
              Tol, fcnt, hmin,  mat_Q3, warn3)

CALL mat_copy(mat_Q1,mat_integral)
CALL mat_axpy((1._dp,0._dp),'N',mat_Q2,mat_integral)
CALL mat_axpy((1._dp,0._dp),'N',mat_Q3,mat_integral)

SELECT CASE (warn)
  CASE 1
CALL cp_warn(__LOCATION__, "Minimum step size reached; singularity possible")
  CASE 2
CALL cp_warn(__LOCATION__, "Maximum function count exceeded. singularity likely")
END SELECT

END SUBROUTINE quad


RECURSIVE SUBROUTINE quadstep(my_negf_env, a,b, mat_fa, mat_fc, mat_fb, Tol, fcnt, hmin, mat_Q, warn)
!in-out
TYPE(negf_env_obj),INTENT(IN)   :: my_negf_env
REAL(KIND=dp),INTENT(IN)        :: a,b
TYPE(mat_z_obj),INTENT(IN)      :: mat_fa, mat_fb, mat_fc
REAL(kind=dp),INTENT(IN)        :: Tol
INTEGER,INTENT(INOUT)           :: fcnt        
INTEGER, INTENT(OUT)            :: warn
TYPE(kind=dp),INTENT(OUT)       :: mat_Q
!others
INTEGER,PARAMETER               :: maxfcnt=10000
REAL(TYPE=dp)                   :: h,c
REAL(TYPE=dp)                   :: NElements
REAL(TYPE=dp),DIMENSION(1:2)    :: X
TYPE(mat_z_obj),DIMENSION(1:2)  :: mat_Y
TYPE(mat_z_obj)                 :: mat_fd, mat_fe, mat_Q1,mat_Q2,&
                                   mat_Qac, matQcb
INTEGER                         :: Nrows
INTEGER                         :: warnac,warncb
!Evaluate integrand twice in interior of subinterval [a,b]
h = b - a
c = (a + b)/2._dp
X(1) = (a + c)/2._dp
X(2) = (c + b)/2._dp
CALL calc_GR(my_negf_env,energy = X(1),Gretarded = mat_Y(1))
CALL calc_GR(my_negf_env,energy = X(2),Gretarded = mat_Y(2))
fcnt = fcnt + 2
mat_fd = mat_Y(1)
mat_fe = mat_Y(2)

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

IF (fcnt .GT. maxfcnt)
  !Maximum function count exceeded; singularity likely.
  warn = 2
  RETURN
END IF

!check mat_nrows (mat_Q) for the optimality!!! Maybe should be passed as an input?
Nrows = mat_nrows(mat_Q)
NElements = REAL((Nrows*Nrows), KIND=dp)
IF ((mat_norm(mat_axpy((-1._dp,0._dp),'N',mat_Q2,mat_Q))/NElements) &
                                                       .LE. Tol) THEN
  !Accuracy over this subinterval is acceptable
  warn = 0
  RETURN
END IF

IF ((ABS(h).LT.hmin).OR.(c.EQ.a).OR.(c.EQ.b))
  !Minimum step size reached; singularity possible.
  warn = 1
  RETURN
END IF

!Subdivide into two subintervals
CALL quadstep(my_negf_env, a, c, mat_fa, mat_fd, mat_fc,&
              Tol, fcnt, hmin,  mat_Qac, warnac)
CALL quadstep(my_negf_env, c, b, mat_fc, mat_fe, mat_fb,&
              Tol, fcnt, hmin,  mat_Qcb, warncb)
CALL mat_copy(mat_Qac,mat_Q)
CALL mat_axpy((1._dp,0._dp),'N',mat_Qcb,mat_Q)
warn = MAX(warnac,warncb)
END SUBROUTINE quadstep

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
END FUNCTION 

END MODULE integral

