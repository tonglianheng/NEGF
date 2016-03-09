MODULE negf_env_types

  USE kinds, ONLY: dp
  USE matrix_types, ONLY: mat_d_obj, &
                          mat_z_obj, &
                          mat_create, &
                          mat_release

#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! public types
  PUBLIC :: negf_env_obj

  ! public methods
  PUBLIC :: negf_env_associate, &
            negf_env_create, &
            negf_env_get, &
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
     INTEGER :: nE
     REAL(KIND=dp) :: eps_E, E_min, E_max, eps_Sancho, dE
     TYPE(mat_d_obj), DIMENSION(:), POINTER :: H_S, S_S, H_LS, S_LS, &
                                               H_L_onsite, S_L_onsite, &
                                               H_L_hopping, S_L_hopping
     TYPE(mat_z_obj), DIMENSION(:), POINTER :: GR_S, Sigma_L, Gamma_L, GR0_L
  END TYPE negf_env_data

  TYPE negf_env_obj
     TYPE(negf_env_data), POINTER, PRIVATE :: obj => NULL()
  END TYPE negf_env_obj

CONTAINS

  SUBROUTINE negf_env_create(negf_env, nterminals, E_min, E_max, nE, eps_E, eps_Sancho)
    TYPE(negf_env_obj), INTENT(INOUT) :: negf_env
    INTEGER, INTENT(IN), OPTIONAL :: nterminals
    INTEGER, INTENT(IN), OPTIONAL :: nE
    REAL(KIND=dp), INTENT(IN), OPTIONAL :: eps_E, eps_Sancho, E_min, E_max

    CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_create', &
                                   routineP = moduleN//':'//routineN

    CPASSERT(.NOT.ASSOCIATED(negf_env%obj))
    ! default settings
    negf_env%obj%nterminals = 2
    negf_env%obj%nE = 1000
    negf_env%obj%E_min = -1.0
    negf_env%obj%E_max = 1.0
    negf_env%obj%eps_E = 1.E-05_dp
    negf_env%obj%eps_Sancho = 1.E-05_dp
    ! modified settings
    IF (PRESENT(nterminals)) negf_env%obj%nterminals = nterminals
    IF (PRESENT(nE))         negf_env%obj%nE = nE
    IF (PRESENT(E_min))      negf_env%obj%E_min = E_min
    IF (PRESENT(E_max))      negf_env%obj%E_max = E_max
    IF (PRESENT(eps_E))      negf_env%obj%eps_E = eps_E
    IF (PRESENT(eps_Sancho)) negf_env%obj%eps_Sancho = eps_Sancho
    ! calculate energy step
    dE = (negf_env%obj%E_max - negf_env%obj%E_min) / real(nE,KIND=dp)
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


  SUBROUTINE negf_env_associate(a, b)
    TYPE(negf_env_obj), INTENT(IN) :: a
    TYPE(negf_env_obj), INTENT(INOUT) :: b

    CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_associate', &
      routineP = moduleN//':'//routineN

    CALL negf_env_release(b)
    b%obj => a%obj
    CALL negf_env_retain(b)
  END SUBROUTINE negf_env_associate


  SUBROUTINE negf_env_get(negf_env, &
                          nterminals, &
                          eps_E, &
                          eps_Sancho, &
                          H_S, &
                          S_S, &
                          GR_S, &
                          H_LS, &
                          S_LS, &
                          H_L, &
                          S_L, &
                          Sigma_L, &
                          Gamma_L, &
                          GR0_L)
    TYPE(negf_env_obj), INTENT(INOUT) :: negf_env
    INTEGER, INTENT(OUT), OPTIONAL :: nterminals
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: eps_E
    TYPE(mat_d_obj), INTENT(OUT), OPTIONAL :: H_S, S_S
    TYPE(mat_z_obj), INTENT(OUT), OPTIONAL :: GR_S
    TYPE(mat_d_obj), DIMENSION(:), POINTER, OPTIONAL :: H_LS, S_LS, H_L, S_L
    TYPE(mat_z_obj), DIMENSION(:), POINTER, OPTIONAL :: Sigma_L, Gamma_L
    TYPE(mat_z_obj), DIMENSION(:), POINTER, OPTIONAL :: GR0_L

    CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_get', &
                                   routineP = moduleN//':'//routineN

    CPASSERT(ASSOCIATED(negf_env%obj))
    ! scalars
    IF (PRESENT(nterminals)) nterminals = negf_env%obj%nterminals
    IF (PRESENT(eps_E)) eps_E = negf_env%obj%eps_E
    IF (PRESENT(eps_Sancho)) eps_Sancho = negf_env%obj%eps_Sancho
    ! scattering region matrices
    IF (PRESENT(H_S)) CALL mat_associate(negf_env%obj%H_S, H_S)
    IF (PRESENT(S_S)) CALL mat_associate(negf_env%obj%S_S, S_S)
    IF (PRESENT(GR_S)) CALL mat_associate(negf_env%obj%GR_S, GR_S)
    ! interface region matrices
    IF (PRESENT(H_LS)) CALL mat_associate(negf_env%H_LS, H_LS)
    IF (PRESENT(S_LS)) CALL mat_associate(negf_env%S_LS, S_LS)
    ! leads matrices
    IF (PRESENT(H_L_onsite)) H_L_onsite => negf_env%obj%H_L_onsite
    IF (PRESENT(H_L_hopping)) H_L_hopping => negf_env%obj%H_L_hopping
    IF (PRESENT(S_L_onsite)) S_L_onsite => negf_env%obj%S_L_onsite
    IF (PRESENT(S_L_onsite)) S_L_hopping => negf_env%obj%S_L_hopping
    IF (PRESENT(Sigma_L)) Sigma_L => negf_env%obj%Sigma_L
    IF (PRESENT(Gamma_L)) Gamma_L => negf_env%obj%Gamma_L
    IF (PRESENT(GR0_L)) GR0_L => negf_env%obj%GR0_L
  END SUBROUTINE negf_env_get


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
          CALL mat_release(negf_env%obj%H_S)
          CALL mat_release(negf_env%obj%S_S)
          CALL mat_release(negf_env%obj%GR_S)
          IF (ASSOCIATED(negf_env%obj%H_LS)) THEN
             DO ii = 1, SIZE(negf_env%obj%H_LS)
                CALL mat_release(negf_env%obj%H_LS(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%H_L_onsite)) THEN
             DO ii = 1, SIZE(negf_env%obj%H_L_onsite)
                CALL mat_release(negf_env%obj%H_L_onsite(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%H_L_hopping)) THEN
             DO ii = 1, SIZE(negf_env%obj%H_L_hopping)
                CALL mat_release(negf_env%obj%H_L_onsite(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%S_LS)) THEN
             DO ii = 1, SIZE(negf_env%obj%S_LS)
                CALL mat_release(negf_env%obj%S_LS(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%S_L_onsite)) THEN
             DO ii = 1, SIZE(negf_env%obj%S_L_onsite)
                CALL mat_release(negf_env%obj%S_L_onsite(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%S_L_hopping)) THEN
             DO ii = 1, SIZE(negf_env%obj%S_L_hopping)
                CALL mat_release(negf_env%obj%S_L_hopping(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%Sigma_L)) THEN
             DO ii = 1, SIZE(negf_env%obj%Sigma_L)
                CALL mat_release(negf_env%obj%Sigma_L(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%Gamma_L)) THEN
             DO ii = 1, SIZE(negf_env%obj%Gamma_L)
                CALL mat_release(negf_env%obj%Gamma_L(ii))
             END DO
          END IF
          IF (ASSOCIATED(negf_env%obj%GR0_L)) THEN
             DO ii = 1, SIZE(negf_env%obj%GR0_L)
                CALL mat_release(negf_env%obj%GR0_L(ii))
             END DO
          END IF
          negf_env%obj%ref_count = 0
          DEALLOCATE(negf_env%obj)
       END IF
    END IF
  END SUBROUTINE negf_env_release

  SUBROUTINE calc_lead_SE(negf_env, energy, iLead, Sigma)
    TYPE(negf_env_obj), INTENT(IN) :: negf_env
    REAL(KIND=dp), INTENT(IN) :: energy
    INTEGER, INTENT(IN) :: iLead
    TYPE(mat_z_obj), INTENT(INOUT) :: Sigma

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_lead_SE', &
                                   routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp) :: EE
    TYPE(mat_z_obj) :: GR_lead_surface, &
                       work1, work2

    CALL mat_release(Sigma)
    CALL mat_create(Sigma, &
                    mat_ncols(negf_env%obj%H_LS(iLead)), &
                    mat_ncols(negf_env%obj%H_LS(iLead)))
    EE = CMPLX(energy, negf_env%obj%eps_E, KIND=dp)
    CALL surface_GR_sancho(GR_lead_surface, &
                           EE, &
                           negf_env%obj%H_L_onsite(iLead), &
                           negf_env%obj%H_L_hopping(iLead), &
                           negf_env%obj%S_L_onsite(iLead), &
                           negf_env%obj%S_L_hopping(iLead), &
                           negf_env%obj%eps_Sancho)
    CALL mat_copy(negf_env%obj%H_LS(iLead), work1)
    CALL mat_scale(work, (-1.0_dp,0.0_dp))
    CALL mat_axpy(EE, 'N', &
                  negf_env%obj%S_LS(iLead), &
                  work1)
    CALL mat_create(work2, &
                    mat_nrows(GR_lead_surface), &
                    mat_ncols(work1))
    CALL mat_mult('N', 'N', (1.0_dp,0.0_dp), &
                  GR_lead_surface, work1, work2)
    CALL mat_mult('H', 'N', (1.0_dp,0.0_dp), &
                  work1, work2, Sigma)
    ! cleanup
    CALL mat_release(work1)
    CALL mat_release(work2)
    CALL mat_release(GR_lead_surface)
  END SUBROUTINE calc_lead_SE

  SUBROUTINE calc_lead_gamma(negf_env, Sigma, Gamma)
    TYPE(negf_env_obj), INTENT(IN) :: negf_en
    TYPE(mat_z_obj), INTENT(IN) :: Sigma
    TYPE(mat_z_obj), INTENT(INOUT) :: Gamma

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_lead_gamma', &
                                   routineP = moduleN//':'//routineN

    TYPE(mat_z_obj) :: Sigma

    CALL mat_release(Gamma)
    CALL mat_copy(Sigma, Gamma)
    CALL mat_scale(Gamma, (0.0_dp,1.0_dp))
    CALL mat_axpy((0.0_dp,-1.0_dp), 'H', Sigma, Gamma)
    ! cleanup
    CALL mat_release(Sigma)
  END SUBROUTINE calc_lead_gamma

  SUBROUTINE calc_transmission(negf_env, energy, transmission, GRetarded)
    TYPE(negf_env_obj), INTENT(IN) :: negf_env
    REAL(KIND=dp), INTENT(IN) :: energy
    REAL(KIND=dp), INTENT(OUT) :: transmission
    TYPE(mat_z_obj), INTENT(INOUT) :: GRetarded

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_scatter_GR', &
                                   routineP = moduleN//':'//routineN

    TYPE(mat_z_obj), DIMENSION(:), ALLOCATABLE :: Sigma
    TYPE(mat_z_obj) :: work1, work2, Gamma
    INTEGER :: ii
    COMPLEX(KIND=dp) :: EE

    ! calculate retarded Green function
    CALL mat_release(GRetarded)
    EE = CMPLX(energy, negf_env%obj%eps_E, KIND=dp)
    ALLOCATE(Sigma(negf_env%obj%nterminals))
    CALL mat_real_to_complex(negf_env%H_S, work1)
    CALL mat_scale(work1, (-1.0_dp,0.0_dp))
    DO ii = 1, negf_env%obj%nterminals
       CALL calc_lead_SE(negf_env, energy, ii, Sigma(ii))
       CALL mat_axpy((-1.0_dp,0.0_dp), 'N', Sigma(ii), work1)
    END DO
    CALL mat_axpy(EE, 'N', negf_env%obj%S_S, work1)
    CALL mat_inv_LU(work, GRetarded)

    ! calculate transmission coefficient
    CALL mat_create(work2, mat_nrows(work1), mat_ncols(work1))
    CALL calc_lead_gamma(negf_env, Sigma(1), Gamma)
    CALL mat_mult('N', 'N', (1.0_dp,0.0_dp), Gamma, GRetarded, &
                  (0.0_dp,0.0_dp), work1)
    CALL calc_lead_gamma(negf_env, Sigma(2), Gamma)
    CALL mat_mult('N', 'N', (1.0_dp,0.0_dp), work1, Gamma, &
                  (0.0_dp,0.0_dp), work2)
    CALL mat_mult('N', 'H', (1.0_dp,0.0_dp), work2, GRetarded, &
                  (0.0_dp,0.0_dp), work1)
    transmission = mat_trace(work1)

    ! cleanup
    CALL mat_release(work1)
    CALL mat_release(work2)
    DO ii = 1, SIZE(Sigma)
       CALL mat_release(Sigma(ii))
    END DO
    DEALLOCATE(Sigma)
  END SUBROUTINE calc_transmission

END MODULE negf_env_types
