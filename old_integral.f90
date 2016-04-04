MODULE integral

USE kinds, ONLY: dp
USE matrix_types!, ONLY: mat_d_obj, &
                !        mat_z_obj, &
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


  INTERFACE quad
     MODULE PROCEDURE quad_d
     MODULE PROCEDURE quad_z
  END INTERFACE quad
 
CONTAINS


SUBROUTINE quad_d(my_negf_env, &
                a,b, &
                Tol, &
                mat_integral,&
                fcnt,&
                func)
  
  TYPE(negf_env_obj) :: my_negf_env
  TYPE(mat_z_obj), INTENT(OUT) :: mat_integral
  REAL(KIND=dp), INTENT(IN)    :: a,b
  REAL(KIND=dp), INTENT(IN)  :: Tol
  INTEGER, INTENT(OUT),OPTIONAL :: fcnt


  INTERFACE
  SUBROUTINE func(my_negf_env, energy, output)
  USE matrix_types, ONLY: mat_z_obj
  USE negf_env_types, ONLY: negf_env_obj
  USE KINDS, ONLY: dp
  REAL(KIND = dp),INTENT(IN)    :: energy
  TYPE(negf_env_obj),INTENT(IN) :: my_negf_env
  TYPE(mat_z_obj),INTENT(INOUT) :: output
  END SUBROUTINE func
  END INTERFACE

 !technical variables
  INTEGER                        :: ii
  REAL(KIND=dp)                  :: h,hmin
  REAL(KIND=dp),   DIMENSION(7)  :: X
  TYPE(mat_z_obj), DIMENSION(7)  :: mat_Y
  REAL(KIND=dp)                  :: x1
  TYPE(mat_z_obj)                :: mat_Q1, mat_Q2, mat_Q3
  INTEGER                        :: warn, warn1, warn2, warn3
!default settings
 ! IF (.NOT.(PRESENT(Tol)))  Tol = 1.E-06_dp

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
hmin = eps(b-a)/1024._dp;
CALL quadstep_d(my_negf_env,X(1), X(3), mat_Y(1), mat_Y(2), mat_Y(3),&
              Tol, fcnt, hmin,  mat_Q1, warn1,func)
CALL quadstep_d(my_negf_env,X(3), X(5), mat_Y(3), mat_Y(4), mat_Y(5),&
              Tol, fcnt, hmin,  mat_Q2, warn2,func)
CALL quadstep_d(my_negf_env,X(5), X(7), mat_Y(5), mat_Y(6), mat_Y(7),&
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
END SUBROUTINE quad_d


RECURSIVE SUBROUTINE quadstep_d(my_negf_env, a,b, mat_fa, mat_fc, mat_fb, Tol, fcnt, hmin, mat_Q, warn,func)
!in-out
TYPE(negf_env_obj),INTENT(IN)   :: my_negf_env
REAL(KIND=dp),INTENT(IN)        :: a,b
TYPE(mat_z_obj),INTENT(IN)      :: mat_fa, mat_fb, mat_fc
REAL(KIND=dp),INTENT(IN)        :: Tol
INTEGER,INTENT(INOUT)           :: fcnt        
INTEGER, INTENT(OUT)            :: warn
TYPE(mat_z_obj),INTENT(OUT)     :: mat_Q

  INTERFACE
  SUBROUTINE func(my_negf_env, energy, output)
  USE matrix_types, ONLY: mat_z_obj
  USE negf_env_types, ONLY: negf_env_obj
  USE KINDS, ONLY: dp
  REAL(KIND = dp),INTENT(IN) :: energy
  TYPE(negf_env_obj),INTENT(IN) :: my_negf_env
  TYPE(mat_z_obj),INTENT(INOUT) :: output
  END SUBROUTINE func
  END INTERFACE

!others
INTEGER,PARAMETER               :: maxfcnt=100000
REAL(KIND=dp)                   :: hmin,h,c
REAL(KIND=dp)                   :: NElements
REAL(KIND=dp),DIMENSION(1:2)    :: X
!TYPE(mat_z_obj),DIMENSION(1:2)    :: mat_Y
TYPE(mat_z_obj)                 :: mat_fd, mat_fe, mat_Q1,mat_Q2,&
                                   mat_Qac, mat_Qcb, mat_difQ
INTEGER                         :: Nrows,ii
INTEGER                         :: warnac,warncb
!Evaluate integrand twice in interior of subinterval [a,b]
h = b - a
c = (a + b)/2._dp
X(1) = (a + c)/2._dp
X(2) = (c + b)/2._dp
CALL func(my_negf_env, X(1),mat_fd)
CALL func(my_negf_env, X(2),mat_fe)
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
  RETURN
END IF

!check mat_nrows (mat_Q) for the optimality!!! Maybe should be passed as an input?
Nrows = mat_nrows(mat_Q)
NElements = REAL((Nrows), KIND=dp)
 CALL mat_copy(mat_Q,mat_difQ)
CALL mat_axpy((-1._dp,0._dp),'N',mat_Q2,mat_difQ)

IF ((mat_norm(mat_difQ)) .LE. Tol) THEN
  !Accuracy over this subinterval is acceptable
  
    warn = 0
    
    !cleanup
    CALL mat_release(mat_fd)
    CALL mat_release(mat_fe)
    CALL mat_release(mat_Q1)
    CALL mat_release(mat_Q2)
    CALL mat_release(mat_difQ)
    RETURN

END IF

IF ((ABS(h).LT.hmin).OR.(c.EQ.a).OR.(c.EQ.b)) THEN
  !Minimum step size reached; singularity possible.
  warn = 1
  RETURN
END IF

!Subdivide into two subintervals
CALL quadstep_d(my_negf_env, a, c, mat_fa, mat_fd, mat_fc,&
              Tol, fcnt, hmin,  mat_Qac, warnac,func)
CALL quadstep_d(my_negf_env, c, b, mat_fc, mat_fe, mat_fb,&
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


END SUBROUTINE quadstep_d

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


SUBROUTINE quad_z(my_negf_env, &
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
  COMPLEX(KIND = dp),INTENT(IN)    :: EE
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
CALL quadstep_z(my_negf_env,X(1), X(3), mat_Y(1), mat_Y(2), mat_Y(3),&
              Tol, fcnt, hmin,  mat_Q1, warn1,func)
CALL quadstep_z(my_negf_env,X(3), X(5), mat_Y(3), mat_Y(4), mat_Y(5),&
              Tol, fcnt, hmin,  mat_Q2, warn2,func)
CALL quadstep_z(my_negf_env,X(5), X(7), mat_Y(5), mat_Y(6), mat_Y(7),&
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
END SUBROUTINE quad_z


RECURSIVE SUBROUTINE quadstep_z(my_negf_env, a,b, mat_fa, mat_fc, mat_fb, Tol, fcnt, hmin, mat_Q, warn,func)
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
  COMPLEX(KIND = dp),INTENT(IN) :: EE
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
CALL func(my_negf_env, X(1),mat_fd)
CALL func(my_negf_env, X(2),mat_fe)
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
IF ((mat_norm(mat_difQ)) .LE. Tol) THEN
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
CALL quadstep_z(my_negf_env, a, c, mat_fa, mat_fd, mat_fc,&
              Tol, fcnt, hmin,  mat_Qac, warnac,func)
CALL quadstep_z(my_negf_env, c, b, mat_fc, mat_fe, mat_fb,&
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

END SUBROUTINE quadstep_z

END MODULE integral

