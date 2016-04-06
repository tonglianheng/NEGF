MODULE cp2k_negf_interface



CONTAINS

  SUBROUTINE construct_para_env_sub(negf_env, para_env, para_env_sub, blacs_env_sub)


    CALL negf_env_get(negf_env=negf_env, nprocs_per_E=nprocs_per_E)
    color_sub = para_env%mepos/negf_env%nprocs_per_E
    CALL mp_comm_split_direct(para_env%group, comm_sub, colour_sub)
    CALL cp_para_env_create(para_env_sub, comm_sub)
    ! create the corresponding BLACS environments
    blacs_grid_layout = BLACS_GRID_SQUARE
    blacs_repeatable = .TRUE.
    CALL cp_blacs_env_create(blacs_env_sub, &
                             para_env_sub, &
                             blacs_grid_layout, &
                             )



  END SUBROUTINE construct_negf_mp_env


  SUBROUTINE get_onsite_matrix(negf_env, dbcsr_mat, onsite_mat, name)
    TYPE(negf_env_obj), INTENT(IN)  :: negf_env
    TYPE(cp_dbcsr_type), INTENT(IN) :: dbcsr_mat
    TYPE(cp_fm_type), POINTER :: onsite_mat
    CHARACTER(len=*), INTENT(in) :: name

    TYPE(cp_fm_struct_type), POINTER :: fmstruct
    TYPE(cp_fm_type), POINTER :: fm_mat
    INTEGER, DIMENSION(:), POINTER :: scatter_atoms, row_blk_size, &
                                      col_blk_size, row_blk_offset, &
                                      col_blk_offset
    INTEGER :: blacs_context_all, iatom_global, nrows, ncols, &
               nrows_onsite, ncols_onsite 

    ! convert DBCSR_mat to fm format
    CALL copy_dbcsr_to_fm(dbcsr_mat, fm_mat)

    ! get information on scattering atoms and BLACS envs
    CALL negf_env_get(negf_env=negf_env, &
                      scatter_atoms=scatter_atoms, &
                      blacs_env=blacs_env, &
                      blacs_context_all=blacs_context_all)

    ! get the correct sub matrix dimensions
    n_scatter_atoms = SIZE(scatter_atoms)
    CPASSERT(n_scatter_atoms .GE. 1)
    CALL cp_dbcsr_get_info(dbcsr_mat, &
                           row_blk_size=row_blk_size, &
                           col_blk_size=col_blk_size, &
                           row_blk_offset=row_blk_offset, &
                           col_blk_offset=col_blk_offset)
    nrows_onsite = 0
    ncols_onsite = 0
    DO ii = 1, n_scatter_atoms
       iatom_global = scatter_atoms(ii)
       nrows_onsite = nrows_onsite + row_blk_size(iatom_global)
       ncols_onsite = ncols_onsite + col_blk_size(iatom_global)
    END DO

    ! create onsite_mat
    CALL cp_fm_struct_create(fmstruct=fmstruct, &
                             para_env=para_env, &
                             context=blacs_env, &
                             nrows_global=nrows_onsite, &
                             ncols_global=ncols_global)
    CALL cp_fm_create(onsite_mat, fmstruct, name=name)
    CALL cp_fm_struct_release(fm_struct)

    ! copy matrix data
    DO ii = 1, n_scatter_atoms
       iatom_global = scatter_atoms(ii)
       irow_start_src = row_blk_offset(iatom_global)
       icol_start_src = col_blk_offset(iatom_global)
       nrows = row_blk_size(iatom_global)
       ncols = col_blk_size(iatom_global)
       irow_start_dest = irow_start_dest + nrows
       icol_start_dest = icol_start_dest + ncols
       ! get the relevant elements from fm and fill onsite_mat
       CALL cp_fm_to_fm_submat_general(fm_mat, &
                                       onsite_mat, &
                                       nrows, &
                                       ncols, &
                                       irow_start_src, &
                                       icol_start_src, &
                                       irow_start_dest, &
                                       icol_start_dest, &
                                       blacs_context_all)
    END DO
  END SUBROUTINE get_onsite_matrix


END MODULE cp2k_negf_interface
