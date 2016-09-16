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
                          mat_associate,&
                          mat_real, &
                          mat_imag, &
                          mat_read
  USE sancho_method, ONLY: surface_GR_sancho 
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
            !calc_GR, &
            calc_G, &
            calc_G_f2,&
            calc_G_Gamma2_G_df,&
            calc_G_f1,&
            calc_G_less,&
            calc_GR_f1,&
            calc_P_G_less,&
            calc_P_GR_f1,&
            negf_env_read_matrices

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
!> \param eps_Sancho     : tolerance of convergence for Sancho's method
!> \prarm nsteps_Sanncho : maximum allowed number of steps for Sancho's method
!> \param id_nr, ref_count : Object book-keeping counters
!> \param u_L            : u_L(i) corresponds to -e*mu_L(i), mu_L(i) - 
!>                           chemical potential of the i-th lead
!> \param EF_L             : EF_L(i) Fermi level in the i-th lead
!> \param EF_S             : Fermi level in the scattering region
!> \param Temperature      : of electron gas
!> \param gamma_int        : length of the L-contour
!> \param delta_int        : height of the L-contour
!> \param Elow_int         : max lowest energy level of leads/sys
!>                           in respect to the correspondent EF
!> \param Eup_int          : upper integration limit
!> \param alpha_R_int      : coef., increases the radius of C-contour (>1)
! ************************************************************************
  TYPE negf_env_data
     INTEGER :: id_nr, ref_count
     INTEGER :: nterminals
     REAL(KIND=dp) :: eps_E, eps_Sancho
     INTEGER :: nsteps_Sancho
     TYPE(mat_d_obj) :: H_S, S_S  
     TYPE(mat_d_obj), DIMENSION(:), POINTER :: H_LS, S_LS, &
                                               H_L_onsite, S_L_onsite, &
                                               H_L_hopping, S_L_hopping
     REAL(KIND=dp), DIMENSION(:), POINTER ::   EF_L, u_L
     REAL(KIND=dp)  :: EF_S, Temperature
     !_int means these are parameters of integrals
     REAL(KIND=dp) :: gamma_int, delta_int, Elow_int, Eup_int, alpha_R_int
     ! TYPE(mat_z_obj), DIMENSION(:), POINTER :: GR_S, Sigma_L, Gamma_L, GR0_L
  END TYPE negf_env_data

  TYPE negf_env_obj
     TYPE(negf_env_data), POINTER, PUBLIC :: obj => NULL()
  END TYPE negf_env_obj

CONTAINS

  SUBROUTINE negf_env_create(negf_env, nterminals, eps_E, eps_Sancho,&
             EF_L, EF_S, u_L, Temperature,&
             gamma_int, delta_int, Elow_int, Eup_int, alpha_R_int)
    TYPE(negf_env_obj), INTENT(INOUT)   :: negf_env
    INTEGER, INTENT(IN), OPTIONAL       :: nterminals
    REAL(KIND=dp), INTENT(IN), OPTIONAL :: eps_E, eps_Sancho
    REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL :: EF_L, u_L
    REAL(KIND=dp), INTENT(IN), OPTIONAL :: EF_S, Temperature,&
             gamma_int, delta_int, Elow_int, Eup_int, alpha_R_int
    CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_create', &
                                   routineP = moduleN//':'//routineN

    CPASSERT(.NOT.ASSOCIATED(negf_env%obj))
    ALLOCATE(negf_env%obj)
    ! default settings
    negf_env%obj%nterminals = 2
    negf_env%obj%eps_E = 1.E-05_dp
    negf_env%obj%eps_Sancho = 1.E-05_dp
    negf_env%obj%nsteps_Sancho = 100
    negf_env%obj%EF_S = 0._dp
    negf_env%obj%Temperature = 2.6E-2_dp  !it is 300 K in eV.
    negf_env%obj%gamma_int   = 1._dp      !it is in eV.
    negf_env%obj%delta_int   = 1._dp      !it is in eV.
    negf_env%obj%Elow_int    = -10._dp    !it is in eV.
    negf_env%obj%Eup_int     = 10._dp     !it is in eV.
    negf_env%obj%alpha_R_int = 2._dp 
    ALLOCATE(negf_env%obj%u_L(nterminals))    
    negf_env%obj%u_L(1) = 0.0_dp
    negf_env%obj%u_L(2) = 0.0_dp
    AllOCATE(negf_env%obj%EF_L(nterminals))
    negf_env%obj%EF_L(1) = 0.0_dp
    negf_env%obj%EF_L(2) = 0.0_dp
    ! modified settings
    IF (PRESENT(nterminals)) negf_env%obj%nterminals = nterminals
    IF (PRESENT(eps_E))      negf_env%obj%eps_E = eps_E
    IF (PRESENT(eps_Sancho)) negf_env%obj%eps_Sancho = eps_Sancho
    IF (PRESENT(EF_S))       negf_env%obj%EF_S = EF_S
    IF (PRESENT(Temperature)) negf_env%obj%Temperature = Temperature
    IF (PRESENT(Temperature)) negf_env%obj%gamma_int = gamma_int
    IF (PRESENT(Temperature)) negf_env%obj%delta_int = delta_int
    IF (PRESENT(Temperature)) negf_env%obj%Elow_int = Elow_int
    IF (PRESENT(Temperature)) negf_env%obj%Eup_int = Eup_int
    IF (PRESENT(Temperature)) negf_env%obj%alpha_R_int = alpha_R_int
    IF (PRESENT(nsteps_Sancho)) negf_env%obj%nsteps_Sancho = nsteps_Sancho
    ! allocate arrays
    ALLOCATE(negf_env%obj%H_LS(nterminals))
    ALLOCATE(negf_env%obj%S_LS(nterminals))
    ALLOCATE(negf_env%obj%H_L_onsite(nterminals))
    ALLOCATE(negf_env%obj%H_L_hopping(nterminals))
    ALLOCATE(negf_env%obj%S_L_onsite(nterminals))
    ALLOCATE(negf_env%obj%S_L_hopping(nterminals))
    ! set bias
    !negf_env%obj%u_L = 0._dp
    IF (PRESENT(u_L)) THEN
       CPASSERT(SIZE(u_L) .EQ. negf_env%obj%nterminals)
       negf_env%obj%u_L = u_L
    END IF
    IF (PRESENT(EF_L)) THEN
       CPASSERT(SIZE(EF_L) .EQ. negf_env%obj%nterminals)
       negf_env%obj%EF_L = EF_L
    END IF
    ! ALLOCATE(negf_env%obj%Sigma_L(nterminals))
    ! ALLOCATE(negf_env%obj%Gamma_L(nterminals))
    ! ALLOCATE(negf_env%obj%GR0_L(nterminals))
    ! ALLOCATE(negf_env%obj%GR0_S(nterminals))

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
                          H_LS, &
                          S_LS, &
                          H_L_onsite, &
                          S_L_onsite,&
                          H_L_hopping,&
                          S_L_hopping,&
                          EF_L,&
                          EF_S,&
                          u_L,&
                          Temperature,&
                          gamma_int,&
                          delta_int,&
                          Elow_int,&
                          Eup_int,&
                          alpha_R_int)
     
     TYPE(negf_env_obj), INTENT(IN) :: negf_env
     INTEGER,         INTENT(OUT), OPTIONAL :: nterminals
     INTEGER :: ii
     REAL(KIND=dp),   INTENT(OUT), OPTIONAL :: eps_E, eps_Sancho
     TYPE(mat_d_obj), INTENT(OUT), OPTIONAL :: H_S, S_S
     TYPE(mat_d_obj), DIMENSION(:), POINTER, OPTIONAL :: H_LS, S_LS, &
                                                         H_L_onsite,  S_L_onsite, &
                                                         H_L_hopping, S_L_hopping
     REAL(KIND=dp), INTENT(OUT), DIMENSION(:), POINTER, OPTIONAL :: &
          EF_L, u_L, &
          gamma_int, delta_int, Elow_int, Eup_int, alpha_R_int
     
     REAL(KIND=dp), INTENT(OUT), OPTIONAL :: EF_S, Temperature
     
     CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_get', &
                                    routineP = moduleN//':'//routineN

     CPASSERT(ASSOCIATED(negf_env%obj))
     ! scalars
     IF (PRESENT(nterminals))    nterminals = negf_env%obj%nterminals
     IF (PRESENT(eps_E))         eps_E      = negf_env%obj%eps_E
     IF (PRESENT(eps_Sancho))    eps_Sancho = negf_env%obj%eps_Sancho
     IF (PRESENT(EF_S))          EF_S       = negf_env%obj%EF_S
     IF (PRESENT(Temperature))   Temperature= negf_env%obj%Temperature
     IF (PRESENT(gamma_int))   Temperature= negf_env%obj%gamma_int
     IF (PRESENT(delta_int))   Temperature= negf_env%obj%delta_int
     IF (PRESENT(Elow_int))   Temperature= negf_env%obj%Elow_int
     IF (PRESENT(Eup_int))   Temperature= negf_env%obj%Eup_int
     IF (PRESENT(alpha_R_int))   Temperature= negf_env%obj%alpha_R_int
     
     ! scattering region matrices
     IF (PRESENT(H_S)) CALL mat_associate(negf_env%obj%H_S, H_S)
     IF (PRESENT(S_S)) CALL mat_associate(negf_env%obj%S_S, S_S)
    ! IF (PRESENT(GR_S)) CALL mat_associate(negf_env%obj%GR_S, GR_S)

    ! interface region matrices
    IF (PRESENT(H_LS)) H_LS => negf_env%obj%H_LS
    IF (PRESENT(S_LS)) S_LS => negf_env%obj%S_LS
    ! leads matrices
    IF (PRESENT(H_L_onsite)) H_L_onsite => negf_env%obj%H_L_onsite
    IF (PRESENT(H_L_hopping)) H_L_hopping => negf_env%obj%H_L_hopping
    IF (PRESENT(S_L_onsite)) S_L_onsite => negf_env%obj%S_L_onsite
    IF (PRESENT(S_L_hopping)) S_L_hopping => negf_env%obj%S_L_hopping
    IF (PRESENT(EF_L)) EF_L => negf_env%obj%EF_L
    IF (PRESENT(u_L)) u_L => negf_env%obj%u_L
    ! IF (PRESENT(Sigma_L)) Sigma_L => negf_env%obj%Sigma_L
    ! IF (PRESENT(Gamma_L)) Gamma_L => negf_env%obj%Gamma_L
    ! IF (PRESENT(GR0_L)) GR0_L => negf_env%obj%GR0_L
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
          ! CALL mat_release(negf_env%obj%GR_S)
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
          !IF (ASSOCIATED(negf_env%obj%u_L)) THEN
          !   DO ii = 1, SIZE(negf_env%obj%u_L)
          !      release(negf_env%obj%u_L(ii))
          !   END DO
          !END IF
          ! IF (ASSOCIATED(negf_env%obj%Sigma_L)) THEN
          !    DO ii = 1, SIZE(negf_env%obj%Sigma_L)
          !       CALL mat_release(negf_env%obj%Sigma_L(ii))
          !    END DO
          ! END IF
          ! IF (ASSOCIATED(negf_env%obj%Gamma_L)) THEN
          !    DO ii = 1, SIZE(negf_env%obj%Gamma_L)
          !       CALL mat_release(negf_env%obj%Gamma_L(ii))
          !    END DO
          ! END IF
          ! IF (ASSOCIATED(negf_env%obj%GR0_L)) THEN
          !    DO ii = 1, SIZE(negf_env%obj%GR0_L)
          !       CALL mat_release(negf_env%obj%GR0_L(ii))
          !    END DO
          ! END IF
          negf_env%obj%ref_count = 0
          DEALLOCATE(negf_env%obj)
       END IF
    END IF
  END SUBROUTINE negf_env_release

  SUBROUTINE calc_lead_SE(negf_env, EE, iLead, Sigma)
     ! Self-energy of the lead iLead
     ! Sigma(EE)= (EE*S_L(i) - H_L(i)) *g_s(i)* (EE*S_L(i) - H_L(i))' 
    TYPE(negf_env_obj), INTENT(IN) :: negf_env
    COMPLEX(KIND=dp), INTENT(IN) :: EE
    INTEGER, INTENT(IN) :: iLead
    TYPE(mat_z_obj), INTENT(INOUT) :: Sigma
    COMPLEX(KIND=dp)  :: EE_temp
    CHARACTER(len=*), PARAMETER :: routineN = 'calc_lead_SE', &
                                   routineP = moduleN//':'//routineN
    TYPE(mat_z_obj) :: GR_lead_surface, &
                       work1, work2, S_S

    CALL mat_release(Sigma)
    CALL mat_create(Sigma, &
                    mat_ncols(negf_env%obj%H_LS(iLead)), &
                    mat_ncols(negf_env%obj%H_LS(iLead)))
    ! take into account u_L(iLead) and EF_L(iLead)
    EE_temp = EE + &
              CMPLX(negf_env%obj%u_L(iLead),  0.0_dp, KIND = dp) + &
              CMPLX(negf_env%obj%EF_L(iLead), 0.0_dp, KIND = dp)

    CALL surface_GR_sancho(GR_lead_surface, &
                           EE_temp, &
                           negf_env%obj%H_L_onsite(iLead), &
                           negf_env%obj%H_L_hopping(iLead), &
                           negf_env%obj%S_L_onsite(iLead), &
                           negf_env%obj%S_L_hopping(iLead), &
                           negf_env%obj%eps_Sancho, &
                           negf_env%obj%nsteps_Sancho)

    CALL mat_real_to_complex(negf_env%obj%H_LS(iLead),work1)

    CALL mat_scale(work1, (-1.0_dp,0.0_dp))

    CALL mat_real_to_complex(negf_env%obj%S_LS(iLead),S_S)
    CALL mat_axpy(EE_temp, 'N', &
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

  SUBROUTINE calc_lead_gamma(negf_env, Sigma, Gamma_lead)
    TYPE(negf_env_obj), INTENT(IN) :: negf_env
    TYPE(mat_z_obj), INTENT(IN)    :: Sigma
    TYPE(mat_z_obj), INTENT(INOUT) :: Gamma_lead

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_lead_gamma', &
                                   routineP = moduleN//':'//routineN

    CALL mat_release(Gamma_lead)
    CALL mat_copy(Sigma, Gamma_lead)
    CALL mat_scale(Gamma_lead, (0.0_dp,1.0_dp))
    CALL mat_axpy((0.0_dp,-1.0_dp), 'H', Sigma, Gamma_lead)
    ! cleanup
  END SUBROUTINE calc_lead_gamma

!  SUBROUTINE calc_transmission(negf_env, EE, transmission, GRetarded)
!    TYPE(negf_env_obj), INTENT(IN) :: negf_env
!    COMPLEX(KIND=dp), INTENT(IN) :: EE
!    REAL(KIND=dp), INTENT(OUT) :: transmission
!    TYPE(mat_z_obj), INTENT(INOUT) :: GRetarded
!
!    CHARACTER(len=*), PARAMETER :: routineN = 'calc_transmission', &
!                                   routineP = moduleN//':'//routineN
!
!    TYPE(mat_z_obj), DIMENSION(:), ALLOCATABLE :: Sigma
!    TYPE(mat_z_obj) :: work1, work2, Gamma, S_S
!    INTEGER :: ii
!    COMPLEX(KIND=dp) :: EE_temp
!
!
!    ! calculate retarded Green function
!    CALL mat_release(GRetarded)
!
!    ALLOCATE(Sigma(negf_env%obj%nterminals))
!    CALL mat_real_to_complex(negf_env%obj%H_S, work1)
!    CALL mat_scale(work1, (-1.0_dp,0.0_dp))
!    DO ii = 1, negf_env%obj%nterminals
!       EE_temp = EE - CMPLX(negf_env%obj%u_L(ii),0.0_dp, KIND = dp)
!       CALL calc_lead_SE(negf_env, EE_temp, ii, Sigma(ii))
!       CALL mat_axpy((-1.0_dp,0.0_dp), 'N', Sigma(ii), work1)
!    END DO
!    CALL mat_real_to_complex(negf_env%obj%S_S, S_S)
!    CALL mat_axpy(EE, 'N', S_S, work1)
!    CALL mat_inv_LU(work1, GRetarded)
!
!
!
!    ! calculate transmission coefficient
!    CALL mat_create(work2, mat_nrows(work1), mat_ncols(work1))
!    CALL calc_lead_gamma(negf_env, Sigma(1), Gamma)
!    CALL mat_mult('N', 'N', (1.0_dp,0.0_dp), Gamma, GRetarded, &
!                  (0.0_dp,0.0_dp), work1)
!    CALL calc_lead_gamma(negf_env, Sigma(2), Gamma)
!    CALL mat_mult('N', 'N', (1.0_dp,0.0_dp), work1, Gamma, &
!                  (0.0_dp,0.0_dp), work2)
!    CALL mat_mult('N', 'C', (1.0_dp,0.0_dp), work2, GRetarded, &
!                  (0.0_dp,0.0_dp), work1)
!    transmission = mat_trace(work1)
!
!    ! cleanup
!    CALL mat_release(work1)
!    CALL mat_release(work2)
!    CALL mat_release(S_S)
!    DO ii = 1, SIZE(Sigma)
!       CALL mat_release(Sigma(ii))
!    END DO
!    DEALLOCATE(Sigma)
!  END SUBROUTINE calc_transmission

    SUBROUTINE calc_G(negf_env, EE, G)
    !> Green function of a complex energy EE
    !> (def.) G(EE) = inv(EE*S - H_S - Sigma_L - Sigma_R)
    TYPE(negf_env_obj), INTENT(IN):: negf_env
    COMPLEX(KIND=dp), INTENT(IN)  :: EE
    TYPE(mat_z_obj), INTENT(INOUT):: G
    CHARACTER(len=*), PARAMETER   :: routineN = 'calc_G', &
                                     routineP = moduleN//':'//routineN

    TYPE(mat_z_obj), DIMENSION(:), ALLOCATABLE :: Sigma
    TYPE(mat_z_obj) :: work, S_S
    INTEGER :: ii
    COMPLEX(KIND=dp) :: EE_temp
   
    ! calculate retarded Green function
    CALL mat_release(G)
    ALLOCATE(Sigma(negf_env%obj%nterminals))
    CALL mat_real_to_complex(negf_env%obj%H_S, work)
    CALL mat_scale(work, (-1.0_dp,0.0_dp))
    DO ii = 1, negf_env%obj%nterminals
       CALL calc_lead_SE(negf_env, EE, ii, Sigma(ii))

    ! energy is already calculated with respect to E-EF_L(i)-U_L(i)
       CALL mat_axpy((-1.0_dp,0.0_dp), 'N', Sigma(ii), work)
    END DO
    EE_temp = EE + CMPLX(negf_env%obj%EF_S, 0.0_dp, KIND = dp)

    CALL mat_real_to_complex(negf_env%obj%S_S, S_S)
    CALL mat_axpy(EE_temp, 'N', S_S, work)
    CALL mat_inv_LU(work, G)

    ! cleanup
    CALL mat_release(work)
    CALL mat_release(S_S)
    DO ii = 1, SIZE(Sigma)
       CALL mat_release(Sigma(ii))
    END DO
    DEALLOCATE(Sigma)
  END SUBROUTINE calc_G

  SUBROUTINE Fermi(negf_env, EE, iLead, fF)
     ! calculate FERMI FUNCTION
     ! (def.) fF = 1/(1+exp( (EE - u_L(iLead))/kT) )
     ! for complex Energy
     ! EF at equilibrium is set to zero, only bias shifts the bands
     TYPE(negf_env_obj), INTENT(IN) :: negf_env
     COMPLEX(KIND=dp), INTENT(OUT) :: fF
     COMPLEX(KIND=dp), INTENT(IN) :: EE
     INTEGER, INTENT(IN) :: iLead
     COMPLEX(KIND = dp) :: tmp, EE_temp, I_cmplx, T_cmplx

     !avoid Inf by reperesenting Fermi function differently:
     !(a) 1/(1+exp(+))
     !(b) exp(-)/(1+exp(-))

     EE_temp = EE &
          - CMPLX(negf_env%obj%u_L(iLead),   0._dp, KIND = dp)

     T_cmplx =  CMPLX(negf_env%obj%Temperature, 0._dp, KIND = dp)

     I_cmplx = (1._dp, 0._dp)

     IF (real(EE_temp, KIND=dp) > 0.0_dp) THEN
        tmp = EXP(-(EE_temp)/T_cmplx)
        fF  = tmp/(I_cmplx + tmp)
     ELSE
        tmp = EXP(EE_temp/T_cmplx)
        fF  =  I_cmplx/(I_cmplx + tmp)
     ENDIF
  END SUBROUTINE Fermi

  ! unused; forseen for nLead >2
  ! SUBROUTINE calc_G_f(negf_env, EE, G_f, iLead)
  !    TYPE(negf_env_obj), INTENT(IN) :: negf_env
  !    TYPE(mat_z_obj), INTENT(INOUT) :: G_f
  !    COMPLEX(KIND=dp), INTENT(IN)   :: EE
  !    INTEGER, INTENT(IN)            :: iLead
     
  !    COMPLEX(KIND=dp)               :: fF
     
  !    CALL Fermi (negf_env, EE, iLead, fF)
  !    CALL calc_G(negf_env, EE, G_f)
  !    CALL mat_scale(G_f, fF)
  ! END SUBROUTINE calc_G_f
  
  SUBROUTINE calc_G_f2(negf_env, EE, G_f)
     ! (def.)  G_f(z) = G(z) * f_R(z)
     ! Properties:  P_eq = 1/(pi) * Im ( \int_{upper plane} G_f(z) dz  )
     ! if P_neq = ...*(f_L - f_R)dE  then P = P_qe + Pneq
     TYPE(negf_env_obj), INTENT(IN) :: negf_env
     TYPE(mat_z_obj), INTENT(INOUT) :: G_f
     COMPLEX(KIND=dp), INTENT(IN)   :: EE
     
     COMPLEX(KIND=dp)               :: fF
     
     ! fF at EE+EF_L(2)+u_L(2)
     CALL Fermi (negf_env, EE, 2, fF)
     ! fF at EE+EF_S
     CALL calc_G(negf_env, EE, G_f)
     CALL mat_scale(G_f, fF)
  END SUBROUTINE calc_G_f2

  SUBROUTINE calc_G_f1(negf_env, EE, G_f)
     ! (def.)  G_f(z) = G(z) * f_L(z)
     ! Properties:  P_eq = 1/(pi) * Im ( \int_{upper plane} G_f(z) dz  )
     ! if P_neq = ...*(f_R - f_L)dE  then P = P_qe + Pneq
     TYPE(negf_env_obj), INTENT(IN) :: negf_env
     TYPE(mat_z_obj), INTENT(INOUT) :: G_f
     COMPLEX(KIND=dp), INTENT(IN)   :: EE
     
     COMPLEX(KIND=dp)               :: fF
     
     ! fF at EE+EF_L(1)+u_L(1)
     CALL Fermi (negf_env, EE, 1, fF)
     ! G at EE+EF_S
     CALL calc_G(negf_env, EE, G_f)
     CALL mat_scale(G_f, fF)
  END SUBROUTINE calc_G_f1

  ! test function. remove in work version!
  SUBROUTINE calc_GR_f1(negf_env, EE, G_f)
     ! (def.) G_f =  G * fL
     ! use it only for EE = real(EE) + i*eta with eta --> 0.
     ! do not use it for G(z)f(z) with imag(z) considerably < or > 0!
     TYPE(negf_env_obj), INTENT(IN) :: negf_env
     TYPE(mat_z_obj), INTENT(INOUT) :: G_f
     COMPLEX(KIND=dp), INTENT(IN)   :: EE
     
     COMPLEX(KIND=dp)               :: fL, EE_temp
     
     ! fF at EE+EF_L(1)+u_L(1)
     EE_temp = CMPLX(real(EE), 0._dp, KIND = dp)
     CALL Fermi (negf_env, EE, 1, fL)
     ! G at EE+EF_S
     CALL calc_G(negf_env, EE, G_f)
     CALL mat_scale(G_f, fL)
  END SUBROUTINE calc_GR_f1

  ! test function. remove in work version!
  SUBROUTINE calc_P_GR_f1(negf_env, EE, P_G_f)
     ! (def.) P_G_f = -1/pi G * fL
     TYPE(negf_env_obj), INTENT(IN) :: negf_env
     TYPE(mat_z_obj), INTENT(INOUT) :: P_G_f
     COMPLEX(KIND=dp), INTENT(IN)   :: EE
     
     COMPLEX(KIND=dp)               :: fL, EE_temp
     
     ! mat_release(P_G_f)
     EE_temp = CMPLX(real(EE), 0._dp, KIND = dp)
     CALL Fermi (negf_env, EE, 1, fF)
     ! G at EE+EF_S
     CALL calc_G(negf_env, EE, P_G_f)
     CALL mat_scale(P_G_f, CMPLX(-1.0_dp/pi,0.0_dp,KIND=dp)*fL)
  END SUBROUTINE calc_P_GR_f1

  !for futher implementations with iLeads > 2
  ! SUBROUTINE calc_G_GammaR_G_df(negf_env, EE, G_GammaR_G_df, iLead)
  !    TYPE(negf_env_obj), INTENT(IN) :: negf_env
  !    TYPE(mat_z_obj), INTENT(INOUT) :: G_GammaR_G_df
  !    COMPLEX(KIND=dp), INTENT(IN)   :: EE
  !    INTEGER, INTENT(IN) :: iLead
  !    COMPLEX(KIND=dp) :: I_cmplx, O_cmplx
  !    TYPE(mat_z_obj)  :: mat, mat1
     
  !    COMPLEX(KIND=dp) :: fFL, fFR
     
  !    CALL Fermi(negf_env, EE,  1, fFL)
  !    CALL Fermi(negf_env, EE,  2, fFR)
  !    fFR = fFR - fFL
  !    CALL mat_release(G_GammaR_G_df)
  !    CALL calc_lead_SE(negf_env, EE, 2, mat)
  !    CALL calc_lead_gamma(negf_env, mat, G_GammaR_G_df)
  !    CALL calc_G(negf_env, EE, mat)
  !    I_cmplx = (1._dp, 0._dp)
  !    O_cmplx = (0._dp, 0._dp)
  !    CALL mat_mult('N', 'N', I_cmplx, mat, G_GammaR_G_df,&
  !         O_cmplx, mat1)
  !    CALL mat_mult('N', 'C', I_cmplx, mat1, mat, &
  !         O_cmplx, G_GammaR_G_df)
     
  !    CALL mat_scale(G_GammaR_G_df, fFR)
     
  !    CALL mat_release(mat)
  !    CALL mat_release(mat1)
  ! END SUBROUTINE calc_G_GammaR_G_df

  SUBROUTINE calc_G_Gamma2_G_df(negf_env, EE, G_Gamma2_G_df)
     ! (def.) G_Gamma2_G_df = G*Gamma_R*G'*(f_R - f_L)
     ! Properties: P_neq = 1/(2*pi*i) * /int G_Gamma2_G_df
     ! P = P_neq + P_eq, where P is a full density matrix
     
     TYPE(negf_env_obj), INTENT(IN) :: negf_env
     TYPE(mat_z_obj), INTENT(INOUT) :: G_Gamma2_G_df
     COMPLEX(KIND=dp), INTENT(IN)   :: EE
     COMPLEX(KIND=dp) :: I_cmplx, O_cmplx, E4Fermi
     REAL(KIND = dp) :: ReE
     TYPE(mat_z_obj)  :: mat, mat1
     
     COMPLEX(KIND=dp) :: fF1, fF2
     ReE = real(EE, KIND = dp)
     !> arg of Fermi here should be a real energy to decrease a mistake
     E4Fermi = CMPLX (ReE, 0._dp, KIND = dp)
     CALL Fermi(negf_env, E4Fermi,  1, fF1)
     CALL Fermi(negf_env, E4Fermi,  2, fF2)
     fF2 = fF2 - fF1
     CALL mat_release(G_Gamma2_G_df)
     !> SE is calculated for EE + EF_L(2) + u_L(2)
     CALL calc_lead_SE(negf_env, EE, 2, mat)
     CALL calc_lead_gamma(negf_env, mat, G_Gamma2_G_df)
     !> G is calculated for EE +E_S
     CALL calc_G(negf_env, EE, mat)
     I_cmplx = (1._dp, 0._dp)
     O_cmplx = (0._dp, 0._dp)
     CALL mat_mult('N', 'N', I_cmplx, mat, G_Gamma2_G_df,&
          O_cmplx, mat1)
     CALL mat_mult('N', 'C', I_cmplx, mat1, mat, &
          O_cmplx, G_Gamma2_G_df)

     CALL mat_scale(G_Gamma2_G_df, fF2)
     
     CALL mat_release(mat)
     CALL mat_release(mat1)
  END SUBROUTINE calc_G_Gamma2_G_df

  !> technical; isn't used in the implementation, used for tests
  SUBROUTINE calc_P_G_less(negf_env, EE, P_G_less)
     ! integrand of the density matrix (the whole)
     ! P == \int P_G_less
     !  (def.) P_G_less = i* G * [ Gamma_L*f_L + Gamma_R*f_R] G' /equiv
     !  (def.) i* [A_L*f_L + A_R*f_R]
     
     TYPE(negf_env_obj), INTENT(IN) :: negf_env
     COMPLEX(KIND = dp), INTENT(IN) :: EE
     TYPE(mat_z_obj), INTENT(INOUT) :: P_G_less
     TYPE(mat_z_obj) :: mat_tmp, mat_tmp1, mat_tmp2, G
     COMPLEX(KIND = dp) :: tmp, EE_temp
     
     CALL mat_release(P_G_less)
     EE_temp = CMPLX(real(EE), 0.0_dp, KIND = dp)
     
     ! Gamma_L*f_L
     CALL calc_lead_SE(negf_env, EE, 1, mat_tmp)
     CALL calc_lead_gamma(negf_env, mat_tmp, mat_tmp1)
     CALL Fermi(negf_env, EE_temp, 1, tmp)
     CALL mat_scale(mat_tmp1, tmp)
     
     ! Gamma_R*f_R
     CALL calc_lead_SE(negf_env, EE, 2, mat_tmp)
     CALL calc_lead_gamma(negf_env, mat_tmp, mat_tmp2)
     CALL Fermi(negf_env, EE_temp, 2, tmp)
     CALL mat_scale(mat_tmp2, tmp)
     
     ! GammaL*fL + GammaR*fR
     CALL mat_axpy((1._dp,0._dp), 'N', mat_tmp1, mat_tmp2)
     
     ! G
     CALL calc_G(negf_env, EE, G)
     
     ! G*[Gamma_L*fL + Gamma_R*fR]
     CALL mat_mult('N', 'N', (1._dp,0._dp), G, mat_tmp2, &
                   (0._dp,0._dp), mat_tmp)
     CALL mat_mult('N', 'C', (1._dp,0._dp), mat_tmp, G, &
                   (0._dp,0._dp), P_G_less)
     
     ! G*[Gamma_L*fL + Gamma_R*fR)*G'
     CALL mat_scale(P_G_less, (0._dp,1._dp))
     
     ! i* ...
     tmp =  CMPLX(0.0_dp, -0.5_dp/pi, KIND = dp)
     CALL mat_scale(P_G_less, tmp)
     
     !clean up
     CALL mat_release(mat_tmp1)
     CALL mat_release(mat_tmp2)
     CALL mat_release(mat_tmp)
     CALL mat_release(G)
     !CALL mat_release(P_G_less_CMPLX)
  END SUBROUTINE calc_P_G_less

  SUBROUTINE calc_G_less(negf_env, EE, G_less)
     ! i* G * [ Gamma_L*f_L + Gamma_R*f_R] G'
     ! technical; not used in the implementation, used for tests
     TYPE(negf_env_obj), INTENT(IN) :: negf_env
     COMPLEX(KIND = dp), INTENT(IN) :: EE
     TYPE(mat_z_obj), INTENT(INOUT) :: G_less
     
     TYPE(mat_z_obj) :: mat_tmp, mat_tmp1, mat_tmp2, G
     COMPLEX(KIND = dp) :: tmp, EE_temp
     
     EE_temp = CMPLX(real(EE), 0.0_dp, KIND = dp)
     
     ! Gamma_L*f_L
     CALL calc_lead_SE(negf_env, EE, 1, mat_tmp)
     CALL calc_lead_gamma(negf_env, mat_tmp, mat_tmp1)
     CALL Fermi(negf_env, EE_temp, 1, tmp)
     CALL mat_scale(mat_tmp1, tmp)
     ! Gamma_R*f_R
     CALL calc_lead_SE(negf_env, EE, 2, mat_tmp)
     CALL calc_lead_gamma(negf_env, mat_tmp, mat_tmp2)
     CALL Fermi(negf_env, EE_temp, 2, tmp)
     CALL mat_scale(mat_tmp2, tmp)
     ! GammaL*fL + GammaR*fR
     CALL mat_axpy((1._dp,0._dp), 'N', mat_tmp1, mat_tmp2)
     ! G
     CALL calc_G(negf_env, EE, G)
     ! G*[Gamma_L*fL + Gamma_R*fR]
     CALL mat_mult('N', 'N', (1._dp,0._dp), G, mat_tmp2, &
                   (0._dp,0._dp), mat_tmp)
     ! G*[Gamma_L*fL + Gamma_R*fR)*G'
     CALL mat_mult('N', 'C', (1._dp,0._dp), mat_tmp, G, &
                   (0._dp,0._dp), G_less)
     ! i* ...
     CALL mat_scale(G_less, (0._dp,1._dp))

     CALL mat_release(mat_tmp1)
     CALL mat_release(mat_tmp2)
     CALL mat_release(mat_tmp)
     CALL mat_release(G)
  END SUBROUTINE calc_G_less

  SUBROUTINE negf_env_read_matrices(negf_env, &
                                    H_S, &
                                    S_S, &
                                    H_L_onsite, &
                                    H_L_hopping, &
                                    S_L_onsite, &
                                    S_L_hopping, &
                                    H_LS, &
                                    S_LS)
     ! read matrices, for tests only
     TYPE(negf_env_obj), INTENT(INOUT) :: negf_env
     CHARACTER(LEN=default_string_length), INTENT(IN) :: H_S, S_S
     CHARACTER(LEN=default_string_length), DIMENSION(:), INTENT(IN) :: H_L_onsite, &
                                                                       H_L_hopping, &
                                                                       S_L_onsite, &
                                                                       S_L_hopping, &
                                                                       H_LS, &
                                                                       S_LS
     CHARACTER(len=*), PARAMETER :: routineN = 'negf_env_read_matrices', &
                                    routineP = moduleN//':'//routineN
     INTEGER :: ii
     DO ii = negf_env%obj%nterminals
        CALL mat_read(negf_env%obj%H_L_onsite(ii), H_L_onsite(ii))
        CALL mat_read(negf_env%obj%H_L_hopping(ii), H_L_hopping(ii))
        CALL mat_read(negf_env%obj%S_L_onsite(ii), S_L_onsite(ii))
        CALL mat_read(negf_env%obj%S_L_hopping(ii), S_L_hopping(ii))
        CALL mat_read(negf_env%obj%H_LS(ii), H_LS(ii))
        CALL mat_read(negf_env%obj%S_LS(ii), S_LS(ii))
     END DO
     CALL mat_read(negf_env%obj%H_S(ii), H_S(ii))
     CALL mat_read(negf_env%obj%S_S(ii), S_S(ii))
  END SUBROUTINE negf_env_read_matrices
 
END MODULE negf_env_types
