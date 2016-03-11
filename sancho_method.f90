MODULE sancho_method

  USE kinds, ONLY: dp
  USE matrix_types, ONLY: mat_d_obj, &
                          mat_z_obj, &
                          mat_inv_lu, &
                          mat_norm, &
                          mat_copy, &
                          mat_create, &
                          mat_release, &
                          mat_nrows, &
                          mat_ncols, &
                          mat_axpy, &
                          mat_mult, &
                          mat_real_to_complex

#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: surface_GR_sancho

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'sancho_method'

CONTAINS

  SUBROUTINE surface_GR_sancho(GR_surface, &
                               energy, &
                               H_onsite, &
                               H_hopping, &
                               S_onsite, &
                               S_hopping, &
                               tolerance, &
                               GR_bulk)
    TYPE(mat_z_obj), INTENT(INOUT) :: GR_surface
    TYPE(mat_z_obj), INTENT(INOUT), OPTIONAL :: GR_bulk
    COMPLEX(KIND=dp), INTENT(IN) :: energy
    TYPE(mat_d_obj), INTENT(IN) :: H_onsite, H_hopping, &
                                   S_onsite, S_hopping
    REAL(KIND=dp), INTENT(IN) :: tolerance

    TYPE(mat_z_obj) :: So, Ho, St, Ht
    TYPE(mat_z_obj) :: A, B, E, ES, B_copy
    TYPE(mat_z_obj) :: GR_surface_A, GR_surface_B

    CALL mat_real_to_complex(S_onsite, So)
    CALL mat_real_to_complex(H_onsite, Ho)
    CALL mat_real_to_complex(S_hopping, St)
    CALL mat_real_to_complex(H_hopping, Ht)

    CALL mat_create(ES, mat_nrows(So), mat_ncols(So))
    CALL mat_axpy((-1.0_dp,0.0_dp), "N", Ho, ES)
    CALL mat_axpy(energy, "N", So, ES)

    CALL mat_copy(ES, E)

    CALL mat_create(A, mat_nrows(St), mat_ncols(St))
    CALL mat_axpy((-1.0_dp,0.0_dp), "N", Ht, A)
    CALL mat_axpy(energy, "N", St, A)

    CALL mat_create(B, mat_ncols(St), mat_nrows(St))
    CALL mat_axpy((-1.0_dp,0.0_dp), "H", Ht, B)
    CALL mat_axpy(CONJG(energy), "H", St, B)

    ! clear GR_surface and recreate to correct size
    CALL mat_release(GR_surface)
    CALL mat_create(GR_surface, mat_nrows(Ho), mat_ncols(Ho))

    ! create some further temp matrices, they may be non-square
    CALL mat_create(GR_surface_A, mat_nrows(GR_surface), mat_ncols(A))
    CALL mat_create(GR_surface_B, mat_nrows(GR_surface), mat_ncols(B))
    CALL mat_create(B_copy, mat_nrows(B), mat_ncols(B))

    DO WHILE (mat_norm(A) + mat_norm(B) .GT. tolerance)
       CALL mat_inv_lu(E, GR_surface)
       CALL mat_mult('N', 'N', (1.0_dp,0.0_dp), GR_surface, A, (0.0_dp,0.0_dp), GR_surface_A)
       CALL mat_mult('N', 'N', (1.0_dp,0.0_dp), GR_surface, B, (0.0_dp,0.0_dp), GR_surface_B)
       CALL mat_mult('N', 'N', (-1.0_dp,0.0_dp), B, GR_surface_A, (1.0_dp,0.0_dp), E)
       CALL mat_mult('N', 'N', (-1.0_dp,0.0_dp), A, GR_surface_B, (1.0_dp,0.0_dp), E)
       CALL mat_mult('N', 'N', (-1.0_dp,0.0_dp), A, GR_surface_B, (1.0_dp,0.0_dp), ES)
       ! use St as temperary storage for new A
       CALL mat_mult('N', 'N', (-1.0_dp,0.0_dp), A, GR_surface_A, (0.0_dp,0.0_dp), St)
       CALL mat_copy(St, A)
       ! use B_copy as temprary storage for new B
       CALL mat_mult('N', 'N', (-1.0_dp,0.0_dp), B, GR_surface_B, (0.0_dp,0.0_dp), B_copy)
       CALL mat_copy(B_copy, B)
    END DO

    CALL mat_inv_lu(ES, GR_surface)

    IF (PRESENT(GR_bulk)) THEN
       CALL mat_inv_lu(E, GR_bulk)
    END IF

    ! cleanup
    CALL mat_release(So)
    CALL mat_release(Ho)
    CALL mat_release(St)
    CALL mat_release(Ht)
    CALL mat_release(A)
    CALL mat_release(B)
    CALL mat_release(E)
    CALL mat_release(ES)
    CALL mat_release(GR_surface_A)
    CALL mat_release(GR_surface_B)
    CALL mat_release(B_copy)
  END SUBROUTINE surface_GR_sancho

END MODULE sancho_method
