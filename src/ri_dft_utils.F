! Routines for the calculation of RIDFT hartree and core contributions
! Based in in parts on RT-BSE code from Stepan Marek
! Currently RIDFT is mostly implemented with QM calculations with static electrostatic point charge embedding
! To be tested in detail (Development version!)

MODULE ri_dft_utils
    USE rt_bse_utils,                    ONLY: get_hartree_noenv, & 
                                              create_hartree_ri_3c, & 
                                              create_2c_coulomb_rep_v
    USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
    USE qs_energy_types,                 ONLY: qs_energy_type
    USE cp_dbcsr_api,                    ONLY: dbcsr_create, dbcsr_clear, dbcsr_desymmetrize, &
                                               dbcsr_p_type, dbcsr_release, dbcsr_type, dbcsr_type_antisymmetric, &
                                               dbcsr_type_symmetric, dbcsr_transposed, dbcsr_type_no_symmetry, &
                                               dbcsr_add
    USE cp_fm_types,                     ONLY: cp_fm_type, cp_fm_create, &
                                               cp_fm_set_all, cp_fm_release
    USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release, &
                                              cp_fm_struct_type
    USE cp_blacs_env,                    ONLY: cp_blacs_env_type
    USE libint_2c_3c,                    ONLY: libint_potential_type
    USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
    USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
    USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
    USE particle_types,                  ONLY: particle_type
    USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
    USE kpoint_types,                    ONLY: kpoint_type
    USE kinds,                           ONLY: dp
    USE cp_cfm_types,                    ONLY: cp_cfm_type, &
                                              cp_cfm_create,&
                                              cp_fm_to_cfm, & 
                                              cp_cfm_release
    USE scf_control_types,               ONLY: scf_control_type
    USE qs_scf_types,                    ONLY: qs_scf_env_type
    USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
    USE cp_control_types,                ONLY: dft_control_type
    USE input_constants,                 ONLY:do_potential_truncated,&
                                              do_potential_coulomb
    USE message_passing,                 ONLY: mp_para_env_type
    
    USE qs_density_mixing_types,         ONLY: gspace_mixing_nr
    USE qs_integral_utils,               ONLY: basis_set_list_setup
    USE libint_wrapper, ONLY: cp_libint_static_init, cp_libint_static_cleanup
    USE cp_dbcsr_operations,             ONLY:copy_dbcsr_to_fm, &
                                              copy_fm_to_dbcsr
    USE cp_dbcsr_contrib,                ONLY: dbcsr_dot
    USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type, &
                                              release_neighbor_list_sets
    USE kpoint_types,                    ONLY: get_kpoint_info
    USE core_ae,                         ONLY: build_core_erf_ridft
    USE pw_poisson_types,                ONLY: pw_poisson_type, &
                                               pw_poisson_none
    USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
    USE pw_methods,                     ONLY: pw_zero

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "ri_dft_utils"

   PUBLIC :: init_ridft,&
             get_hartree_energy, & 
             get_hartree_and_core_contributions_energy
             

CONTAINS
! **************************************************************************************************
!> \brief Checks the inputs are compatible with RIDFT and initializes necessary parameters
!>        (Sets poisson_env%parameters%solver = pw_poisson_none!)
!> \param qs_env qs-env storing most of the needed input information (scf_control,particle_set,...)
!> \param poisson_env Use to set the poisson_env%parameters%solver = pw_poisson_none
!> \date 08.2026
!> \author Jelena Schmitz
! **************************************************************************************************
    SUBROUTINE init_ridft(qs_env, poisson_env)
        CHARACTER(len=*), PARAMETER                        :: routineN = "init_ridft"
        INTEGER                                            :: handle
        TYPE(qs_environment_type), POINTER, INTENT(IN)     :: qs_env
        TYPE(dft_control_type), POINTER                    :: dft_control
        TYPE(pw_poisson_type), POINTER                     :: poisson_env
        TYPE(section_vals_type), POINTER                   :: poisson_section
        LOGICAL                                            :: is_given
        TYPE(qs_scf_env_type), POINTER                     :: scf_env_ridft

        CALL timeset(routineN, handle)
        CALL get_qs_env(qs_env, dft_control = dft_control,scf_env=scf_env_ridft)
        CPASSERT(ASSOCIATED(qs_env))
        poisson_section => section_vals_get_subs_vals(qs_env%input, "DFT%POISSON")
        CALL section_vals_val_get(poisson_section, "POISSON_SOLVER", explicit=is_given)

        !Check QM/MM not being active 
        IF (qs_env%qmmm) THEN
               CPABORT("QM/MM is not implemented with RIDFT, use point charges for embedding instead.")
        END IF
        !Check for do_ddapc == True
        IF (dft_control%qs_control%ddapc_restraint) THEN
            CPABORT("'do_ddapc/ddapc_restraint' is not compatible with RIDFT.")
        END IF
        IF (ASSOCIATED(qs_env%cp_ddapc_ewald)) THEN
            IF (qs_env%cp_ddapc_ewald%do_decoupling .OR. &
                 qs_env%cp_ddapc_ewald%do_qmmm_periodic_decpl .OR. &
                 qs_env%cp_ddapc_ewald%do_solvation) THEN
                    CPABORT("'do_ddapc/cp_ddapc_ewald' is not compatible with RIDFT.")
            END IF
        END IF
        !RIDFT needs no Poisson solver and sets it to none. If one was given raise an error!
        IF (is_given) THEN
               CPABORT("Do not specify a Poisson solver with RIDFT.")
        END IF

        !Remove the default poisson solver
        poisson_env%parameters%solver = pw_poisson_none
        !Check it was sucessfully removed
        IF (.NOT. (poisson_env%parameters%solver .EQ. pw_poisson_none)) THEN
               CPABORT("Do not specify a Poisson solver with RIDFT.")
        END IF

        !Use DIRECT_P_MIXING with ridft
        IF (ASSOCIATED(scf_env_ridft)) THEN
            IF (scf_env_ridft%mixing_method >= gspace_mixing_nr) THEN
                CPABORT("RIDFT does not suppport g-space mixing methods.")
            END IF
        END IF

        CALL timestop(handle)

    END SUBROUTINE init_ridft
! **************************************************************************************************
!> \brief Calculation of the Hartree energy contribution to the total energy in Resolution of Identity DFT
!>        (Does not contain the core parts as the GPW version would (see addition in qs_ks_methods))
!> \param qs_env qs-env storing most of the needed input information (scf_control,particle_set,...)
!> \param rho_ao Input density in ao basis at the current SCF step 
!> \param energy Stores energy%hartree, the result which will be returned 
!> \param v_ao_dbcsr Returns the v_hartree element  directly 
!> \date 08.2026
!> \author Jelena Schmitz
!> \author Stepan Marek [Parts taken from rt_bse.F]
!> \author Jan Wilhelm [Parts taken from gw_utils.F]
! **************************************************************************************************
    SUBROUTINE get_hartree_energy(qs_env, rho_ao, energy,v_ao_dbcsr)
        CHARACTER(len=*), PARAMETER                        :: routineN = "get_hartree_energy"
        INTEGER                                            :: handle
        !Input/Output:
        TYPE(qs_environment_type), POINTER, INTENT(IN)     :: qs_env
        TYPE(qs_energy_type), POINTER, INTENT(INOUT)       :: energy
        TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER,&
                                           INTENT(IN)      :: rho_ao
        TYPE(dbcsr_type), INTENT(INOUT)                    :: v_ao_dbcsr                                   
        !Routine specific:
        TYPE(cp_fm_type), DIMENSION(:), ALLOCATABLE        :: v_fm
        TYPE(cp_fm_struct_type), POINTER                   :: rho_fm_struct
        TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
        TYPE(cp_fm_type), POINTER                          :: mo_coeff_ri
        TYPE(dbcsr_type)                                   :: v_hartree_T, v_hartree_T_unsym
        TYPE(dbcsr_type), ALLOCATABLE                      :: v_dbcsr
        INTEGER, DIMENSION(:),ALLOCATABLE                  :: sizes_RI
        INTEGER                                            :: i, iatom, n_ao, n_atom, &
                                                              n_kind, n_ri, nsgf, unit_nr, &
                                                              ikind
        INTEGER, DIMENSION(:), ALLOCATABLE                 :: i_RI_start_from_atom, &
                                                              i_RI_end_from_atom
        TYPE(libint_potential_type)                        :: ri_metric = libint_potential_type()
        TYPE(gto_basis_set_p_type),DIMENSION(:),&
          ALLOCATABLE                                      :: basis_set_RI_p_type, basis_set_AO_p_type
        TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
        TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
        TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
        TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
        INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
        REAL(KIND=dp),DIMENSION(:,:,:), POINTER            :: int_3c_ptr => NULL()
        REAL(dp)                                           :: e_hartree_ri_local, spin_degeneracy
        TYPE(cp_cfm_type), DIMENSION(:), POINTER           :: rho_ao_ridft, rho_ao_ridft_copy
        TYPE(scf_control_type), POINTER                    :: scf_control
        TYPE(cp_fm_type)                                   :: temp_ridft_fm
        TYPE(dbcsr_type), TARGET                           :: rho_ao_calc, rho_ao_ridft_copy_dbcsr
        TYPE(section_vals_type), POINTER                   :: gw_sec, post_scf_bands_section
        TYPE(dft_control_type), POINTER                    :: dft_control
        TYPE(mp_para_env_type), POINTER                    :: para_env





        CALL timeset(routineN, handle)
        CALL get_qs_env(qs_env, scf_control=scf_control,qs_kind_set=qs_kind_set, particle_set=particle_set,&
                         mos=mos, atomic_kind_set=atomic_kind_set, dft_control=dft_control,&
                         para_env=para_env)

        n_kind = SIZE(qs_kind_set)
        n_atom = SIZE(particle_set)
        CALL get_mo_set(mo_set=mos(1), mo_coeff=mo_coeff_ri, nao=n_ao)
        blacs_env => mo_coeff_ri%matrix_struct%context

        CALL dbcsr_create(rho_ao_calc,'rho_ao_calc', template=rho_ao(1,1)%matrix,&
                              matrix_type=dbcsr_type_no_symmetry)

        CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

        ALLOCATE(basis_set_AO_p_type(n_kind))
        ALLOCATE(basis_set_RI_p_type(n_kind))
        CALL basis_set_list_setup(basis_set_RI_p_type, "RI_AUX", qs_kind_set)
        CALL basis_set_list_setup(basis_set_AO_p_type, "ORB", qs_kind_set)

        !From gw_utils: (get_RI_basis_and_basis_function_indices)
        ALLOCATE (i_RI_start_from_atom(n_atom))
        ALLOCATE (i_RI_end_from_atom(n_atom))
        ALLOCATE(sizes_RI(n_atom))
        sizes_RI = 0
        n_RI = 0
        DO iatom = 1, n_atom
            i_RI_start_from_atom(iatom) = n_RI + 1
            ikind = kind_of(iatom)
            CALL get_qs_kind(qs_kind=qs_kind_set(ikind), nsgf=nsgf, basis_type="RI_AUX")
            sizes_RI(iatom) = nsgf
            n_RI = n_RI + nsgf
            i_RI_end_from_atom(iatom) = n_RI
        END DO

        post_scf_bands_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%BANDSTRUCTURE")
        gw_sec => section_vals_get_subs_vals(post_scf_bands_section, "GW")

        !CALL section_vals_val_get(gw_sec, "CUTOFF_RADIUS_RI", r_val=ri_metric%cutoff_radius)
        ri_metric%cutoff_radius = 10.0_dp
        ri_metric%omega = 0.0_dp
        ri_metric%filename = "t_c_g.dat" !and use do_potential_truncated
        !ri_metric%potential_type = do_potential_truncated
        ri_metric%potential_type = do_potential_coulomb 

        unit_nr = 1

        CALL cp_libint_static_init()
        CALL cp_fm_struct_create(fmstruct=rho_fm_struct, context=blacs_env,&
                                  nrow_global = n_ao, ncol_global = n_ao, para_env= para_env)

        ALLOCATE (rho_ao_ridft(dft_control%nspins))
        ALLOCATE (rho_ao_ridft_copy(dft_control%nspins))

        CALL cp_fm_create(temp_ridft_fm, rho_fm_struct)
        DO i = 1, dft_control%nspins
            CALL cp_cfm_create(rho_ao_ridft(i), matrix_struct=rho_fm_struct)
            CALL cp_cfm_create(rho_ao_ridft_copy(i), matrix_struct=rho_fm_struct)
            CALL copy_dbcsr_to_fm(rho_ao(i, 1)%matrix, temp_ridft_fm)
            CALL cp_fm_to_cfm(msourcer=temp_ridft_fm, mtarget=rho_ao_ridft(i))
            CALL cp_fm_to_cfm(msourcer=temp_ridft_fm, mtarget=rho_ao_ridft_copy(i))
        END DO

        CALL create_hartree_ri_3c(rho_dbcsr = rho_ao_calc, int_3c = int_3c_ptr, &
                                   n_ao=n_ao, n_RI=n_RI,&
                                   basis_set_AO = basis_set_AO_p_type,&
                                   basis_set_RI = basis_set_RI_p_type, &
                                   i_RI_start_from_atom = i_RI_start_from_atom, ri_metric = ri_metric,&
                                   qs_env = qs_env, unit_nr = unit_nr)                      
          
        ALLOCATE (v_fm(dft_control%nspins))

        DO i = 1, dft_control%nspins
            CALL cp_fm_create(v_fm(i), rho_fm_struct) 
            CALL cp_fm_set_all(v_fm(i), 0.0_dp)
        END DO
   
        ALLOCATE (v_dbcsr)

        CALL dbcsr_create(v_ao_dbcsr,'v_ao_dbcsr', template = rho_ao_calc, matrix_type=dbcsr_type_symmetric)

        CALL create_2c_coulomb_rep_v(v_dbcsr, basis_set_RI_p_type, qs_env,&
                                      n_RI, blacs_env, para_env, ri_metric)
       
        !rho_fm Density matrix corresponding to single spin species, in atomic orbital basis
        IF (dft_control%nspins == 1) spin_degeneracy = 2.0_dp
        IF (dft_control%nspins == 2) spin_degeneracy = 1.0_dp
        energy%hartree = 0.0_dp
        DO i = 1, dft_control%nspins
            CALL get_hartree_noenv(v_fm = v_fm(i),&
                                 rho_fm = rho_ao_ridft(i),&
                                 int_3c = int_3c_ptr,&
                                 v_dbcsr = v_dbcsr,&
                                 n_RI = n_RI, sizes_RI = sizes_RI,&
                                 para_env = para_env,&
                                 rho_dbcsr = rho_ao_calc,&
                                 v_ao_dbcsr = v_ao_dbcsr)

            CALL dbcsr_clear(v_ao_dbcsr)

            !Calculate Hartree energy contribution from potential 
            e_hartree_ri_local = 0.0_dp
            
            !Transfer calculated v_hartree from fm to dbcsr
            CALL copy_fm_to_dbcsr(v_fm(i), v_ao_dbcsr)

            !Tranpose the Hartree Pot 
            CALL dbcsr_create(v_hartree_T, 'v_hartree_T', template=v_ao_dbcsr)
            CALL dbcsr_transposed(v_hartree_T,v_ao_dbcsr)
            CALL dbcsr_desymmetrize(v_hartree_T, v_hartree_T_unsym)

            CALL copy_dbcsr_to_fm(rho_ao(1, 1)%matrix, temp_ridft_fm)
            CALL dbcsr_create(rho_ao_ridft_copy_dbcsr, 'rho_ao_ridft_copy_dbcsr', template=rho_ao_calc)
            CALL copy_fm_to_dbcsr(temp_ridft_fm, rho_ao_ridft_copy_dbcsr)

            !Calculate Trace using dot(rho, v_hartree_transpose)
            CALL dbcsr_dot(rho_ao_ridft_copy_dbcsr, &
                  v_hartree_T_unsym, &
                  e_hartree_ri_local)
            e_hartree_ri_local = e_hartree_ri_local *0.5_dp

            energy%hartree = energy%hartree + e_hartree_ri_local
        END DO
        
        CALL dbcsr_release(v_hartree_T)
        CALL dbcsr_release(v_hartree_T_unsym)
        CALL cp_fm_release(temp_ridft_fm)
        CALL cp_fm_release(v_fm)
        CALL cp_fm_struct_release(rho_fm_struct)
        DO i = 1, dft_control%nspins
            CALL cp_cfm_release(rho_ao_ridft(i))
            CALL cp_cfm_release(rho_ao_ridft_copy(i))
        END DO
        DEALLOCATE(basis_set_AO_p_type)
        DEALLOCATE(basis_set_RI_p_type)
        DEALLOCATE(i_RI_start_from_atom)
        DEALLOCATE(i_RI_end_from_atom)
        DEALLOCATE(sizes_RI)
        DEALLOCATE(v_dbcsr)
        DEALLOCATE(int_3c_ptr)
        CALL cp_libint_static_cleanup()
        CALL timestop(handle)
    END SUBROUTINE get_hartree_energy
! **************************************************************************************************
!> \brief Calculation the equivalent of the energy contributions stored in the hartree energy for GPW
!>        in Resolution of Identity DFT framework 
!> \param qs_env qs-env storing most of the needed input information (scf_control,particle_set,...)
!> \param rho_ao Input density in ao basis at the current SCF step 
!> \param energy Stores energy%hartree, the result which will be returned 
!> \param v_ao_dbcsr Returns the v_hartree element  directly 
!> \param ks_matrix Kohn Sham matrix updated with the hartree potential contributions
!> \date 08.2026
!> \author Jelena Schmitz
! **************************************************************************************************
    SUBROUTINE get_hartree_and_core_contributions_energy(qs_env, rho_ao, energy, ks_matrix,&
                                                         rho_tot_gspace,v_hartree_rspace)
        CHARACTER(len=*), PARAMETER                        :: routineN = "get_hartree_and_core_contributions_energy"
        INTEGER                                            :: handle
        !Input/Output:
        TYPE(qs_environment_type), POINTER, INTENT(IN)     :: qs_env
        TYPE(qs_energy_type), POINTER, INTENT(INOUT)       :: energy
        TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER,&
                                           INTENT(IN)      :: rho_ao
        TYPE(dbcsr_type)                                   :: v_ao_dbcsr 
        TYPE(dft_control_type), POINTER                    :: dft_control
        TYPE(neighbor_list_set_p_type), DIMENSION(:), &
            POINTER                                        :: sab_orb_ridft, sac_ae_ridft, sac_ppl_ridft, &
                                                               sac_erf_ridft
        TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
        REAL(KIND=dp)                                      :: e_nuc_erf, e_nn_erf
        INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_ridft
        TYPE(kpoint_type), POINTER                         :: kpoints_ridft
        TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
        TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
        INTEGER                                            :: nimages
        TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
        TYPE(pw_c1d_gs_type)                               :: rho_tot_gspace
        TYPE(pw_r3d_rs_type), POINTER                      :: v_hartree_rspace
        

        CALL timeset(routineN, handle)

        !Zero the created rho_tot_gspace (should not be filled if ridft = True, but can create bugs if not zeroed),
        CALL pw_zero(rho_tot_gspace)
        CALL pw_zero(v_hartree_rspace)

        CALL dbcsr_create(v_ao_dbcsr,'v_ao_dbcsr', template = rho_ao(1,1)%matrix, matrix_type=dbcsr_type_symmetric)

        CALL get_qs_env(qs_env, sab_orb=sab_orb_ridft, sac_ae=sac_ae_ridft, &
                            sac_ppl=sac_ppl_ridft, sac_erf=sac_erf_ridft, kpoints=kpoints_ridft,&
                            qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set, particle_set=particle_set,&
                            dft_control = dft_control)

        !Calculate the hartree energy contribution:
        CALL get_hartree_energy(qs_env, rho_ao, energy, v_ao_dbcsr)

        !Update the ks matrix with the hartree potential:
        CALL dbcsr_add(ks_matrix(1,1)%matrix,v_ao_dbcsr,alpha_scalar=1.0_dp, beta_scalar=1.0_dp)

        !energy%hartree (Poisson) += e_nuc_erf + e_nn_erf - energy%core_self
        NULLIFY (sab_orb_ridft, sac_ae_ridft, sac_ppl_ridft, sac_erf_ridft, cell_to_index_ridft, kpoints_ridft)
        nimages = dft_control%nimages
        !JS no assc check for sab_orb_ridft?
        IF (nimages > 1 .AND. ASSOCIATED(kpoints_ridft)) THEN
            CALL get_kpoint_info(kpoint=kpoints_ridft, cell_to_index=cell_to_index_ridft)
        END IF
        IF (ASSOCIATED(sac_ae_ridft) .OR. ASSOCIATED(sac_ppl_ridft) .OR. ASSOCIATED(sac_erf_ridft)) THEN
            CALL build_core_erf_ridft(ks_matrix, rho_ao, e_nuc_erf, e_nn_erf, &
                                          qs_kind_set, atomic_kind_set, particle_set, &
                                          sab_orb_ridft, sac_ae_ridft, sac_ppl_ridft, sac_erf_ridft, &
                                          nimages, cell_to_index_ridft, "ORB")
                                          
            energy%hartree = energy%hartree + e_nuc_erf + e_nn_erf - energy%core_self
        END IF
        CALL dbcsr_release(v_ao_dbcsr)
        CALL timestop(handle)

    END SUBROUTINE get_hartree_and_core_contributions_energy
END MODULE ri_dft_utils