MODULE negf_env_types

  USE kinds, ONLY: dp
  USE matrix_types, ONLY: mat_d_obj, &
                          mat_z_obj, &
                          mat_create, &
                          mat_delete

#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! public types
  PUBLIC :: negf_env_obj

  ! public methods
  PUBLIC :: negf_env_associate, &
            negf_env_create, &
            negf_env_has_data, &
            negf_env_nullify, &
            ! negf_env_get, &
            negf_env_retain, &
            negf_env_release

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_env_types'
  INTEGER, PRIVATE, SAVE :: last_negf_env_id = 0

! ************************************************************************
!> \brief data type for NEGF calculations
!> \param nterminals : number of leads
!> \param H_S, S_S, GR_S : the hamiltonian, overlap and retarded Green's
!>                         function blocks for the scattering (central)
!>                         region
!> \param H_LS, S_LS     : the hamiltonian and overlap blocks of the
!>                         interface region between the terminals (leads)
!>                         and the scattering region
!> \param H_L, S_L       : The hamiltonian and overlap blocks of the lead
!>                         regions. H_L(i) corresponds to i-th lead
!> \param Sigma_L        : Self energies of the leads
!> \param eps_E          : infinitesmal for retarded Green's function
!>                         GR = (E + i * eps_E)**(-1)
!> \param id_nr, ref_count : Object book-keeping counters
! ************************************************************************
  TYPE negf_env_data
     INTEGER :: id_nr, ref_count
     INTEGER :: nterminals
     REAL(KIND=dp) :: eps_E, E_min, E_max, eps_Sancho
     TYPE(mat_d_obj) :: H_S, S_S
     TYPE(mat_z_obj) :: GR_S
     TYPE(mat_d_obj), DIMENSION(:), POINTER :: H_LS, S_LS, H_L, S_L
     TYPE(mat_z_obj), DIMENSION(:), POINTER :: Sigma_L, Gamma_L
     TYPE(mat_z_obj), DIMENSION(:), POINTER :: GR0_L
  END TYPE negf_env_data

  TYPE negf_env_obj
     TYPE(negf_env_data), POINTER, PRIVATE :: obj
  END TYPE negf_env_obj

CONTAINS


  SUBROUTINE negf_env_create(negf_env, nterminals, E_min, E_max, eps_E, eps_Sancho)
    TYPE(negf_env_obj), INTENT(INOUT) :: negf_env
    INTEGER, INTENT(IN), OPTIONAL :: nterminals
    REAL(KIND=dp), INTENT(IN), OPTIONAL :: eps_E, eps_Sancho, E_min, E_max

    CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_create', &
                                   routineP = moduleN//':'//routineN

    CPASSERT(.NOT.negf_env_has_data(negf_env))
    IF (PRESENT(nterminals)) THEN
       negf_env%obj%nterminals = nterminals
    ELSE
       negf_env%obj%nterminals = 2
    END IF
    IF (PRESENT(E_min)) THEN
       negf_env%obj%E_min = E_min
    ELSE
       negf_env%obj%E_min = -1.0
    END IF
    IF (PRESENT(E_max)) THEN
       negf_env%obj%E_max = E_max
    ELSE
       negf_env%obj%E_max = 1.0
    END IF
    IF (PRESENT(eps_E)) THEN
       negf_env%obj%eps_E = eps_E
    ELSE
       negf_env%obj%eps_E = 1.E-05_dp
    END IF
    IF (PRESENT(eps_Sancho)) THEN
       negf_env%obj%eps_Sancho = eps_Sancho
    ELSE
       negf_env%obj%eps_Sancho = 1.E-05_dp
    END IF
    ! allocate arrays
    ALLOCATE(negf_env%obj)
    ALLOCATE(negf_env%obj%H_LS(nterminals))
    ALLOCATE(negf_env%obj%H_L(nterminals))
    ALLOCATE(negf_env%obj%S_LS(nterminals))
    ALLOCATE(negf_env%obj%S_L(nterminals))
    ALLOCATE(negf_env%obj%S_L(nterminals))
    ALLOCATE(negf_env%obj%Sigma_L(nterminals))
    ALLOCATE(negf_env%obj%Gamma_L(nterminals))
    ALLOCATE(negf_env%obj%GR0_L(nterminals))
    ! book-keeping stuff
    negf_env%obj%ref_count = 1
    negf_env%obj%id_nr = last_negf_env_id + 1
    last_negf_env_id = negf_env%obj%id_nr
  END SUBROUTINE negf_env_create

  PURE FUNCTION negf_env_has_data(negf_env) RESULT(res)
    TYPE(negf_env_obj), INTENT(IN) :: negf_env
    LOGICAL :: res
    res = ASSOCIATED(negf_env%obj)
  END FUNCTION negf_env_has_data

  SUBROUTINE negf_env_associate(a, b)
    TYPE(negf_env_obj), INTENT(OUT) :: a
    TYPE(negf_env_obj), INTENT(IN) :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_associate', &
      routineP = moduleN//':'//routineN
    
    a%obj => b%obj
    CALL negf_env_retain(a)
  END SUBROUTINE negf_env_associate
  

  SUBROUTINE negf_env_nullify(negf_env)
    TYPE(negf_env_obj), INTENT(OUT) :: negf_env
    
    CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_nullify', &
      routineP = moduleN//':'//routineN

    NULLIFY(negf_env%obj)
  END SUBROUTINE negf_env_nullify


  ! SUBROUTINE negf_env_get(negf_env, &
  !                         nterminals, &
  !                         eps_E, &
  !                         H_S, &
  !                         S_S, &
  !                         GR_S, &
  !                         H_LS, &
  !                         S_LS, &
  !                         H_L, &
  !                         S_L, &
  !                         Sigma_L, &
  !                         Gamma_L, &
  !                         GR0_L)
  !   TYPE(negf_env_obj), INTENT(INOUT) :: negf_env
  !   INTEGER, INTENT(OUT), OPTIONAL :: nterminals
  !   REAL(KIND=dp), INTENT(OUT), OPTIONAL :: eps_E
  !   TYPE(mat_d_obj), INTENT(OUT), OPTIONAL :: H_S, S_S
  !   TYPE(mat_z_obj), INTENT(OUT), OPTIONAL :: GR_S
  !   TYPE(mat_d_obj), DIMENSION(:), POINTER, OPTIONAL :: H_LS, S_LS, H_L, S_L
  !   TYPE(mat_z_obj), DIMENSION(:), POINTER, OPTIONAL :: Sigma_L, Gamma_L
  !   TYPE(mat_z_obj), DIMENSION(:), POINTER, OPTIONAL :: GR0_L

  !   CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_get', &
  !                                  routineP = moduleN//':'//routineN
    
  !   IF (PRESENT(nterminals)) nterminals = negf_env%obj%nterminals
  !   IF (PRESENT(eps_E)) eps_E = negf_env%obj%eps_E

  !   IF (PRESENT(H_S)) H_S%p => negf_env%obj%H_S%p
  !   IF (PRESENT(S_S)) S_S%p => negf_env%obj%S_S%p
  !   IF (PRESENT(GR_S)) GR_S%p => negf_env%obj%GR_S%p

  !   IF (PRESENT(H_LS)) H_LS => negf_env%obj%H_LS
  !   IF (PRESENT(S_LS)) S_LS => negf_env%obj%S_LS 
  !   IF (PRESENT(H_L)) H_L => negf_env%obj%H_L
  !   IF (PRESENT(S_L)) S_L => negf_env%obj%S_L
  !   IF (PRESENT(Sigma_L)) Sigma_L => negf_env%obj%Sigma_L
  !   IF (PRESENT(Gamma_L)) Gamma_L => negf_env%obj%Gamma_L
  !   IF (PRESENT(GR0_L)) GR0_L => negf_env%obj%GR0_L
  ! END SUBROUTINE negf_env_get


  SUBROUTINE negf_env_retain(negf_env)
    TYPE(negf_env_obj), INTENT(INOUT) :: negf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_retain', &
                                   routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(negf_env%obj))
    negf_env%obj%ref_count = negf_env%obj%ref_count + 1
  END SUBROUTINE negf_env_retain

  SUBROUTINE negf_env_release(negf_env)
    TYPE(negf_env_obj), INTENT(INOUT) :: negf_env

    CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_release', &
                                   routineP = moduleN//':'//routineN
    
    INTEGER :: ii

    IF (ASSOCIATED(negf_env%obj)) THEN
       CPASSERT(negf_env%obj%ref_count .GT. 0)
       negf_env%obj%ref_count = negf_env%obj%ref_count - 1
       IF (negf_env%obj%ref_count .EQ. 0) THEN
          negf_env%obj%ref_count = 1
          CALL mat_delete(negf_env%obj%H_S)
          CALL mat_delete(negf_env%obj%S_S)
          CALL mat_delete(negf_env%obj%GR_S)
          IF (ASSOCIATED(negf_env%obj%H_LS)) THEN
             DO ii = 1, SIZE(negf_env%obj%H_LS)
                CALL mat_delete(negf_env%obj%H_LS(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%H_L)) THEN
             DO ii = 1, SIZE(negf_env%obj%H_L)
                CALL mat_delete(negf_env%obj%H_L(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%S_LS)) THEN
             DO ii = 1, SIZE(negf_env%obj%S_LS)
                CALL mat_delete(negf_env%obj%S_LS(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%S_L)) THEN
             DO ii = 1, SIZE(negf_env%obj%S_L)
                CALL mat_delete(negf_env%obj%S_L(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%Sigma_L)) THEN
             DO ii = 1, SIZE(negf_env%obj%Sigma_L)
                CALL mat_delete(negf_env%obj%Sigma_L(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%Gamma_L)) THEN
             DO ii = 1, SIZE(negf_env%obj%Gamma_L)
                CALL mat_delete(negf_env%obj%Gamma_L(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%GR0_L)) THEN
             DO ii = 1, SIZE(negf_env%obj%GR0_L)
                CALL mat_delete(negf_env%obj%GR0_L(ii))
             END DO
          END IF
          negf_env%obj%ref_count = 0
          DEALLOCATE(negf_env%obj)
       END IF
    ELSE
       NULLIFY(negf_env%obj)
    END IF
  END SUBROUTINE negf_env_release

END MODULE negf_env_types
