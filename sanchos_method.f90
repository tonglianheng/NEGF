MODULE sancho_method

  USE kinds, ONLY: dp
  USE matrix_types, ONLY: mat_d_obj, &
                          mat_z_obj, &
                          mat_inv_lu, &
                          mat_norm, &
                          mat_copy, &
                          mat_create, &
                          mat_delete

#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sancho_method'

CONTAINS

  SUBROUTINE surface_GR0_sancho(GR0, &
                                energy, &
                                onsite, &
                                hopping, &
                                overlap_on, &
                                overlap_hop, &
                                tolerance, &
                                GR_bulk)
    TYPE(mat_z_obj), INTENT(INOUT) :: GR0
    TYPE(mat_z_obj), INTENT(INOUT), OPTIONAL :: GR_bulk
    COMPLEX(KIND=dp), INTENT(IN) :: energy
    TYPE(mat_d_obj), INTENT(IN) :: onsite, hopping, &
                                   overlap_on, overlap_hop
    REAL(KIND=dp), INTENT(IN) :: tolerance

    TYPE(mat_z_obj) :: So, Ho, St, Ht
    TYPE(mat_z_obj) :: A, B, E, ES

    CALL mat_real_to_complex(overlap_on, So)
    CALL mat_real_to_complex(onsite, Ho)
    CALL mat_real_to_complex(overlap_hop, St)
    CALL mat_real_to_complex(hopping, Ht)

    CALL mat_create(ES, mat_nrows(So), mat_ncols(So))
    CALL mat_axpy((-1.0_dp,0.0_dp), "N", Ho, ES)
    CALL mat_axpy(energy, "N", So, ES)

    CALL mat_copy(ES, E)

    CALL mat_create(A, mat_nrows(St), mat_ncols(St))
    CALL mat_axpy((-1.0_dp,0.0_dp), "N", Ht, A)
    CALL mat_axpy(energy, "N", St, ESt_m_Ht)

    CALL mat_create(B, mat_ncols(St), mat_nrows(St))
    CALL mat_axpy((-1.0_dp,0.0_dp), "H", Ht, B)
    CALL mat_axpy(CONJG(energy), "H", St, B)

    ! clear GR0 and recreate to correct size
    CALL mat_delete(GR0)
    CALL mat_create(GR0, mat_nrows(Ho), mat_ncols(Ho))

    ! create some further temp matrices, they may be non-square
    CALL mat_create(GR0_A, mat_nrows(GR0), mat_ncols(A))
    CALL mat_create(GR0_B, mat_nrows(GR0), mat_ncols(B))
    CALL mat_create(B_copy, mat_nrows(B), mat_ncols(B))

    DO WHILE (mat_norm(A) + mat_norm(B) .GT. tolerance)
       CALL mat_inv_lu(E, GR0)
       CALL mat_mult('N', 'N', (1.0_dp,0.0_dp), GR0, A, (0.0_dp,0.0_dp), GR0_A)
       CALL mat_mult('N', 'N', (1.0_dp,0.0_dp), GR0, B, (0.0_dp,0.0_dp), GR0_B)
       CALL mat_mult('N', 'N', (-1.0_dp,0.0_dp), B, GR0_A, (1.0_dp,0.0_dp), E)
       CALL mat_mult('N', 'N', (-1.0_dp,0.0_dp), A, GR0_B, (1.0_dp,0.0_dp), E)
       CALL mat_mult('N', 'N', (-1.0_dp,0.0_dp), A, GR0_B, (1.0_dp,0.0_dp), ES)
       ! use St as temperary storage for new A
       CALL mat_mult('N', 'N', (-1.0_dp,0.0_dp), A, GR0_A, (0.0_dp,0.0_dp), St)
       CALL mat_copy(St, A)
       ! use B_copy as temprary storage for new B
       CALL mat_mult('N', 'N', (-1.0_dp,0.0_dp), B, GR0_B, (0.0_dp,0.0_dp), B_copy)
       CALL mat_copy(B_copy, B)
    END DO

    CALL mat_inv_lu(ES, GR0)

    IF (PRESENT(GR_bulk)) THEN
       CALL mat_inv_lu(E, GR_bulk)
    END IF

    ! cleanup
    CALL mat_delete(So)
    CALL mat_delete(Ho)
    CALL mat_delete(St)
    CALL mat_delete(Ht)
    CALL mat_delete(A)
    CALL mat_delete(B)
    CALL mat_delete(E)
    CALL mat_delete(ES)
    CALL mat_delete(GR0_A)
    CALL mat_delete(GR0_B)
    CALL mat_delete(B_copy)
  END SUBROUTINE surface_G0_sancho

END MODULE sancho_method
