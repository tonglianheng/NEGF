MODULE negf_env_types

  USE kinds, ONLY: dp
  USE matrix_types, ONLY: mat_d_obj,   &
                          mat_z_obj,   &
                          mat_create,  &
                          mat_release, &
                          mat_nrows,   &
                          mat_ncols,   &
                          mat_trace,   &
                          mat_release, &
                          mat_real_to_complex,&
                          mat_scale,   &
                          mat_axpy,    &
                          mat_inv_lu,  &
                          mat_mult,    &
                          mat_copy,    &
                          mat_associate
  USE sancho_method, ONLY: surface_GR_sancho 
  USE physcon, ONLY: boltzmann
  USE mathconstants, ONLY: pi
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
            negf_env_release, &
            Fermi, &
            calc_GR, &
            calc_G, &
            calc_GR_fL,&
            calc_G_GammaR_G_fRmfL
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
!>                         and the scattering region. H_LS(i) corresponds
!>                         to i-th lead
!> \param H_L, S_L       : The hamiltonian and overlap blocks of the lead
!>                         regions. H_L(i) corresponds to i-th lead
!> \param Sigma_L        : Self energies of the leads. Sigma_L(i) 
!>                         corresponds to i-th lead
!> \param eps_E          : infinitesmal for retarded Green's function
!>                         GR = (E + i * eps_E)**(-1)
!> \param id_nr, ref_count : Object book-keeping counters
!> \param u_L            : are -e*mu_L(i) where mu_L(i) are chemical 
!>                         potentials in the i-th  lead
! ************************************************************************
  TYPE negf_env_data
     INTEGER :: id_nr, ref_count
     INTEGER :: nterminals
     REAL(KIND=dp) :: eps_E, eps_Sancho
     TYPE(mat_d_obj) :: H_S, S_S  
     TYPE(mat_d_obj), DIMENSION(:), POINTER :: H_LS, S_LS, &
                                               H_L_onsite, S_L_onsite, &
                                               H_L_hopping, S_L_hopping
     REAL(KIND=dp), DIMENSION(:), POINTER ::   u_L
  END TYPE negf_env_data

  TYPE negf_env_obj
     TYPE(negf_env_data), POINTER, PUBLIC :: obj => NULL()
  END TYPE negf_env_obj

  INTERFACE Fermi
    MODULE PROCEDURE Fermi_d
    MODULE PROCEDURE Fermi_z
  END INTERFACE Fermi

CONTAINS

  SUBROUTINE negf_env_create(negf_env, nterminals, eps_E, eps_Sancho)
    TYPE(negf_env_obj), INTENT(INOUT) :: negf_env
    INTEGER, INTENT(IN), OPTIONAL :: nterminals
    REAL(KIND=dp), INTENT(IN), OPTIONAL :: eps_E, eps_Sancho
    REAL(KIN=dp), INTENT(IN), DIMENSION(:), POINTER: u_L 
    CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_create', &
                                   routineP = moduleN//':'//routineN

    CPASSERT(.NOT.ASSOCIATED(negf_env%obj))
    ALLOCATE(negf_env%obj)
    ! default settings
    negf_env%obj%nterminals = 2
    negf_env%obj%eps_E = 1.E-05_dp
    negf_env%obj%eps_Sancho = 1.E-05_dp
    ! modified settings
    IF (PRESENT(nterminals)) negf_env%obj%nterminals = nterminals
                          nterminals, &
                          eps_E, &
                          eps_Sancho, &                          
                          H_S, &
                          S_S, &
                          H_LS, &
                          S_LS, &
                          H_L, &
                          S_L, &
                          )
    TYPE(negf_env_obj), INTENT(INOUT) :: negf_env
    INTEGER, INTENT(OUT), OPTIONAL :: nterminals
    REAL(KIND=dp), INTENT(OUT), OPTIONAL :: eps_E, eps_Sancho
    TYPE(mat_d_obj), INTENT(OUT), OPTIONAL :: H_S, S_S
    TYPE(mat_d_obj), DIMENSION(:), POINTER, OPTIONAL :: H_LS, S_LS, H_L, S_L

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
          negf_env%obj%ref_count = 0
          DEALLOCATE(negf_env%obj)
       END IF
    END IF
  END SUBROUTINE negf_env_release

  SUBROUTINE calc_lead_SE_R(negf_env, energy, Sigma_R)
    TYPE(negf_env_obj), INTENT(IN) :: negf_env
    REAL(KIND=dp), INTENT(IN) :: energy
    TYPE(mat_z_obj), INTENT(INOUT) :: Sigma

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_lead_SE_R', &
                                   routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp) :: EE
    TYPE(mat_z_obj) :: GR_lead_surface_R, &
                       work1, work2

    CALL mat_release(Sigma_R)
    CALL mat_create(Sigma_R, &
                    mat_ncols(negf_env%obj%H_LS), &
                    mat_ncols(negf_env%obj%H_LS)
    !include bias
    EE = CMPLX(energy-negf_env%obj%u_L(iLead), negf_env%obj%eps_E, KIND=dp)
    CALL surface_GR_sancho_R(GR_lead_surface_R, &
                           EE, &
                           negf_env%obj%H_R_onsite, &
                           negf_env%obj%H_R_hopping, &
                           negf_env%obj%S_R_onsite, &
                           negf_env%obj%S_R_hopping, &
                           negf_env%obj%eps_Sancho)
    CALL mat_copy(negf_env%obj%H_LS(work1)
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
  END SUBROUTINE calc_lead_SE_R

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
    
    CALL mat_real_to_complex(negf_env%obj%S_LS(iLead),S_S)
    CALL mat_axpy(EE, 'N', &
                  S_S, &
                  work1)
    
    CALL mat_create(work2, &
                    mat_nrows(GR_lead_surface), &
                    mat_ncols(work1))
    CALL mat_mult('N', 'C', (1.0_dp,0.0_dp), &
                  GR_lead_surface, work1,(1.0_dp,0.0_dp), work2)
    CALL mat_mult('N', 'N', (1.0_dp,0.0_dp), &
                  work1, work2,(1.0_dp,0.0_dp), Sigma)
    ! cleanup
    CALL mat_release(work1)
    CALL mat_release(work2)
    CALL mat_release(S_S)
    CALL mat_release(GR_lead_surface)
  END SUBROUTINE calc_lead_SE

  SUBROUTINE calc_lead_SE_z(negf_env, EE, iLead, Sigma)
    TYPE(negf_env_obj), INTENT(IN) :: negf_env
    COMPLEX(KIND=dp), INTENT(IN) :: EE
    INTEGER, INTENT(IN) :: iLead
