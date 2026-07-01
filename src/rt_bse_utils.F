! Subroutine from rt_bse.F moved to utils to be useabel in ri-dft without creating circular dependencies

MODULE rt_bse_utils
   USE post_scf_bandstructure_types, ONLY: post_scf_bandstructure_type
   USE cp_fm_types, ONLY: cp_fm_type, cp_fm_create, &
                          cp_fm_set_all, cp_fm_release
   USE cp_cfm_types, ONLY: cp_cfm_type, &
                           cp_cfm_to_fm
   USE kinds, ONLY: dp
   USE cp_dbcsr_contrib, ONLY: dbcsr_print
   USE cp_dbcsr_api, ONLY: dbcsr_type, &
                           dbcsr_clear, &
                           dbcsr_iterator_type, &
                           dbcsr_iterator_start, &
                           dbcsr_iterator_stop, &
                           dbcsr_iterator_blocks_left, &
                           dbcsr_iterator_next_block, &
                           dbcsr_reserve_blocks, &
                           dbcsr_get_num_blocks, &
                           dbcsr_get_block_p, &
                           dbcsr_get_info, & 
                           dbcsr_p_type, &
                           dbcsr_create, & 
                           dbcsr_release, & 
                           dbcsr_copy
   USE cp_dbcsr_operations, ONLY: copy_dbcsr_to_fm, &
                                  copy_fm_to_dbcsr,&
                                 cp_dbcsr_sm_fm_multiply
   USE message_passing, ONLY: mp_para_env_type
   USE basis_set_types, ONLY: gto_basis_set_p_type
   USE libint_2c_3c, ONLY: libint_potential_type
   USE qs_environment_types, ONLY: qs_environment_type, &
                                   get_qs_env
   USE gw_integrals, ONLY: build_3c_integral_block

   !JS for 2c int creation func:
   USE qs_neighbor_list_types, ONLY: neighbor_list_set_p_type, &
                                     release_neighbor_list_sets
   USE input_section_types, ONLY: section_vals_get,&
                                  section_vals_get_subs_vals, &
                                  section_vals_type, & 
                                  section_vals_val_get
   USE cp_blacs_env,                    ONLY:  cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type 
   USE parallel_gemm_api, ONLY: parallel_gemm     
   USE qs_tensors, ONLY: build_2c_integrals, &
                         build_2c_neighbor_lists    
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_type  
   USE cp_fm_basic_linalg, ONLY: cp_fm_invert    
   
   !JS for mat_munu part:
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types,                  ONLY: molecule_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                              distribution_2d_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_distribution_type, &
                                              dbcsr_type_no_symmetry
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE distribution_methods,            ONLY: distribute_molecules_2d
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_dist2d_to_dist
   USE cp_blacs_env,                    ONLY: BLACS_GRID_SQUARE,&
                                              cp_blacs_env_create
   USE particle_methods,                ONLY: get_particle_set
   USE qs_integral_utils,               ONLY: basis_set_list_setup

   USE cell_types,                      ONLY: cell_type,&
                                              get_cell

   

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = "rt_bse_utils"

   !#:include "rt_bse_macros.fypp"

   PUBLIC :: create_hartree_ri_3c, &
             get_hartree_noenv, &
             create_2c_coulomb_rep_v

CONTAINS

! **************************************************************************************************
!> \brief Calculates the Hartree matrix in the atomic orbital basis, given a density matrix, in local arrays
!>        Calculates the values for single spin species present in given rho
!> \param v_fm Hartree potential in atomic orbital basis - is overwritten by the updated potential
!> \param rho_fm Density matrix corresponding to single spin species, in atomic orbital basis
!> \param int_3c Previously allocated array (best to use create_hartree_ri_3c) for 3c integrals
!> \param v_dbcsr Previously calculated 2c Coulomb repulsion between RI orbitals
!> \param n_RI Number of RI basis orbitals
!> \param sizes_RI Number of RI basis orbitals per atom
!> \param para_env MPI Parallel environment (used for summation across ranks)
!> \param rho_dbcsr Previously created dbcsr matrix, used as workspace
!> \param v_ao_dbcsr Previously created dbcsr matrix, used as workspace
!> \author Stepan Marek
!> \date 01.2025
! **************************************************************************************************
   SUBROUTINE get_hartree_noenv(v_fm, rho_fm, int_3c, v_dbcsr, n_RI, sizes_RI, para_env, rho_dbcsr, v_ao_dbcsr)
      TYPE(cp_fm_type)                                   :: v_fm
      TYPE(cp_cfm_type), INTENT(IN)                      :: rho_fm
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: int_3c
      TYPE(dbcsr_type)                                   :: v_dbcsr
      INTEGER                                            :: n_RI
      INTEGER, DIMENSION(:)                              :: sizes_RI
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type)                                   :: rho_dbcsr, v_ao_dbcsr
      CHARACTER(len=*), PARAMETER                        :: routineN = "get_hartree"
      TYPE(dbcsr_iterator_type)                          :: iterator_matrix
      INTEGER                                            :: i, j, k, n, nblocks, ind_1, ind_2, row_offset, col_offset, &
                                                            row_size, col_size, j_n_AO, k_n_AO, i_n_RI, &
                                                            ri_offset, ind_i, handle
      REAL(kind=dp), DIMENSION(:), ALLOCATABLE           :: Pvector, Qvector
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: block_matrix
      INTEGER                                            :: nblkrows_local, nblkcols_local, j_blk, k_blk, j_offset, k_offset
      INTEGER, DIMENSION(:), POINTER                     :: local_blk_rows, local_blk_cols
      LOGICAL                                            :: found
      !LOGICAL                                            :: transposed

      MARK_USED(i_n_RI)
      MARK_USED(ri_offset)
      MARK_USED(ind_i)

      ! No memory optimisation so far - calculate all 3cs an all ranks
      ! Importantly - dbcsr blocks are ordered by atoms - i.e. ethene with 6 atoms will have 6x6 block structure
      ! Number of basis states on each basis set is known is post_scf_bandstructure env

      CALL timeset(routineN, handle)

      ! Allocate the Q and Pvector on each rank
      ALLOCATE (Qvector(n_RI), source=0.0_dp)
      ALLOCATE (Pvector(n_RI), source=0.0_dp)

      ! First step - analyze the structure of copied dbcsr matrix on all ranks
      CALL dbcsr_clear(rho_dbcsr)
      ! Only the real part of the density matrix contributes
      ! Use v_fm as workspace
      CALL cp_cfm_to_fm(msource=rho_fm, mtargetr=v_fm)
      CALL copy_fm_to_dbcsr(v_fm, rho_dbcsr)
      j_offset = 0
      CALL dbcsr_get_info(rho_dbcsr, nblkrows_local=nblkrows_local, nblkcols_local=nblkcols_local, &
                          local_rows=local_blk_rows, local_cols=local_blk_cols)
      DO j_blk = 1, nblkrows_local
         k_offset = 0
         DO k_blk = 1, nblkcols_local
            ! Check whether we can retrieve the rho block
            ! TODO : Handle transposed case?
            CALL dbcsr_get_block_p(rho_dbcsr, local_blk_rows(j_blk), local_blk_cols(k_blk), &
                                   block=block_matrix, found=found, row_size=row_size, col_size=col_size)
            ! If the block is not found, then the density matrix here has below threshold values
            IF (.NOT. found) CYCLE
            ! With the block retrieved, add its contributions to the Q-vector
            !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j,k) &
            !$OMP SHARED(n_RI, row_size, col_size, Qvector, int_3c, j_offset, k_offset, block_matrix)
            DO i = 1, n_RI
               DO j = 1, row_size
                  DO k = 1, col_size
                     Qvector(i) = Qvector(i) + int_3c(j_offset + j, k_offset + k, i)*block_matrix(j, k)
                  END DO
               END DO
            END DO
            !$OMP END PARALLEL DO
            ! Increment k-offset - setup for the next block
            k_offset = k_offset + col_size
         END DO
         ! Increments the j-offset - row_size is carried over from the last iteration
         j_offset = j_offset + row_size
      END DO
      ! Now, each rank has contributions from D_jk within its scope
      ! Need to sum over different ranks to get the total vector on all ranks
      CALL para_env%sum(Qvector)
      ! Once this is done, Pvector is current on all ranks
      ! Continue with V_PQ summation
      nblocks = dbcsr_get_num_blocks(v_dbcsr)
      CALL dbcsr_iterator_start(iterator_matrix, v_dbcsr)
      DO n = 1, nblocks
         ! TODO : Try OMP parallelisation over different blocks - expect many more available speedup for large systems
         CALL dbcsr_iterator_next_block(iterator_matrix, ind_1, ind_2, block_matrix, &
                                        row_offset=row_offset, col_offset=col_offset, row_size=row_size, col_size=col_size)
         ! TODO : Better names for RI
         j_n_AO = sizes_RI(ind_1)
         k_n_AO = sizes_RI(ind_2)
         ! The allocations are as follows
         !$OMP PARALLEL DO DEFAULT(none) PRIVATE(j,k) &
         !$OMP SHARED(block_matrix, Pvector, Qvector,j_n_AO,k_n_AO,row_offset,col_offset)
         DO j = 1, j_n_AO
            DO k = 1, k_n_AO
               Pvector(j + row_offset - 1) = Pvector(j + row_offset - 1) + block_matrix(j, k)*Qvector(k + col_offset - 1)
            END DO
         END DO
         !$OMP END PARALLEL DO
      END DO
      CALL dbcsr_iterator_stop(iterator_matrix)
      ! Again, make sure that the P vector is present on all ranks
      CALL para_env%sum(Pvector)
      ! Now, for the final trick, iterate over local blocks of v_ao_dbcsr to get the Hartree as dbcsr, then convert to fm
      ! TODO : Clear or set blocks to zero
      ! CALL dbcsr_clear(v_ao_dbcsr)
      j_offset = 0
      DO j_blk = 1, nblkrows_local
         k_offset = 0
         DO k_blk = 1, nblkcols_local
            ! Check whether we can retrieve the rho block
            ! TODO : Handle transposed case?
            CALL dbcsr_get_block_p(v_ao_dbcsr, local_blk_rows(j_blk), local_blk_cols(k_blk), &
                                   block=block_matrix, found=found, row_size=row_size, col_size=col_size)
            ! If the block is not found, reserve it
            IF (.NOT. found) THEN
               ! Reservations
               CALL dbcsr_reserve_blocks(v_ao_dbcsr, local_blk_rows(j_blk:j_blk), local_blk_cols(k_blk:k_blk))
               ! Rerun the getter to get the new block
               CALL dbcsr_get_block_p(v_ao_dbcsr, local_blk_rows(j_blk), local_blk_cols(k_blk), &
                                      block=block_matrix, found=found, row_size=row_size, col_size=col_size)
            END IF
            ! With the block retrieved, contract with the P vector
            !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j,k) &
            !$OMP SHARED(row_size, col_size, n_RI, block_matrix, Pvector, int_3c, j_offset, k_offset)
            DO j = 1, row_size
               DO k = 1, col_size
                  block_matrix(j, k) = 0.0_dp
                  DO i = 1, n_RI
                     block_matrix(j, k) = block_matrix(j, k) + Pvector(i)*int_3c(j_offset + j, k_offset + k, i)
                  END DO
               END DO
            END DO
            !$OMP END PARALLEL DO
            ! Increment k-offset - setup for the next block
            k_offset = k_offset + col_size
         END DO
         ! Increments the j-offset - row_size is carried over from the last iteration
         j_offset = j_offset + row_size
      END DO
      ! Since P vector was present on all the ranks, v_dbcsr_ao has the complete Hartree result
      ! copy_dbcsr_to_fm should set all values in v_fm to zero
      CALL copy_dbcsr_to_fm(v_ao_dbcsr, v_fm)
      DEALLOCATE (Qvector)
      DEALLOCATE (Pvector)

      CALL timestop(handle)
   END SUBROUTINE get_hartree_noenv
! **************************************************************************************************
!> \brief Separated method for allocating the 3c integrals for RI Hartree
!> \note RI method calculates the Hartree contraction without the use of DBT, as it cannot emulate vectors
!> \param rho_dbcsr matrix used for the description of shape of 3c array
!> \param int_3c 3-center integral array to be allocated and filled
!> \param n_ao Number of atomic orbitals
!> \param n_RI Number of auxiliary RI orbitals
!> \param basis_set_AO AO basis set
!> \param basis_set_RI RI auxiliary basis set
!> \param i_RI_start_from_atom Array of indices where functions of a given atom in RI basis start
!> \param unit_nr Unit number used for printing information about the size of int_3c
!> \author Stepan Marek
!> \date 01.2025
! **************************************************************************************************
   SUBROUTINE create_hartree_ri_3c(rho_dbcsr, int_3c, n_ao, n_RI, basis_set_AO, basis_set_RI, &
                                   i_RI_start_from_atom, ri_metric, qs_env, unit_nr)
      TYPE(dbcsr_type)                                  :: rho_dbcsr
      REAL(kind=dp), DIMENSION(:, :, :), POINTER          :: int_3c
      INTEGER                                           :: n_ao, n_RI
      TYPE(gto_basis_set_p_type), DIMENSION(:)          :: basis_set_AO, &
                                                           basis_set_RI
      INTEGER, DIMENSION(:)                             :: i_RI_start_from_atom
      TYPE(libint_potential_type)                       :: ri_metric
      TYPE(qs_environment_type), POINTER                :: qs_env
      INTEGER                                           :: unit_nr
      REAL(kind=dp)                                     :: size_mb
      INTEGER                                           :: nblkrows_local, &
                                                           nblkcols_local, &
                                                           i_blk_local, &
                                                           j_blk_local, &
                                                           nrows_local, &
                                                           ncols_local, &
                                                           col_local_offset, &
                                                           row_local_offset, &
                                                           start_col_index, &
                                                           end_col_index, &
                                                           start_row_index, &
                                                           end_row_index
      INTEGER, DIMENSION(:), POINTER                    :: local_blk_rows, &
                                                           local_blk_cols, &
                                                           row_blk_size, &
                                                           col_blk_size
      PRINT *, 'create_hartree_ri_3c started'                                                     
      ! TODO : Implement option/decision to not precompute all the 3c integrals
      size_mb = REAL(n_ao, kind=dp)*REAL(n_ao, kind=dp)*REAL(n_RI, kind=dp)* &
                REAL(STORAGE_SIZE(size_mb), kind=dp)/8.0_dp/1024.0_dp/1024.0_dp
      IF (unit_nr > 0) WRITE (unit_nr, '(A44,E32.2E3,A4)') &
         " RTBSE| Approximate size of the 3c integrals", size_mb, " MiB"  
      ! Get the number of block rows and columns
      CALL dbcsr_get_info(rho_dbcsr, nblkrows_local=nblkrows_local, nblkcols_local=nblkcols_local)
      ! Get the global indices of local rows and columns
      CALL dbcsr_get_info(rho_dbcsr, local_rows=local_blk_rows, local_cols=local_blk_cols)
      ! Get the sizes of all blocks
      CALL dbcsr_get_info(rho_dbcsr, row_blk_size=row_blk_size, col_blk_size=col_blk_size)  
      ! Get the total required local rows and cols
      nrows_local = 0
      DO i_blk_local = 1, nblkrows_local
         nrows_local = nrows_local + row_blk_size(local_blk_rows(i_blk_local))
      END DO
      ncols_local = 0
      DO j_blk_local = 1, nblkcols_local
         ncols_local = ncols_local + col_blk_size(local_blk_cols(j_blk_local))
      END DO
      ! Allocate the appropriate storage
      ALLOCATE (int_3c(nrows_local, ncols_local, n_RI))
      ! Fill the storage with appropriate values, block by block
      row_local_offset = 1
      DO i_blk_local = 1, nblkrows_local
         col_local_offset = 1
         DO j_blk_local = 1, nblkcols_local
            start_row_index = row_local_offset
            end_row_index = start_row_index + row_blk_size(local_blk_rows(i_blk_local)) - 1
            start_col_index = col_local_offset
            end_col_index = start_col_index + col_blk_size(local_blk_cols(j_blk_local)) - 1
            CALL build_3c_integral_block(int_3c(start_row_index:end_row_index, &
                                                start_col_index:end_col_index, &
                                                1:n_RI), &
                                         qs_env, potential_parameter=ri_metric, &
                                         basis_j=basis_set_AO, basis_k=basis_set_AO, &
                                         basis_i=basis_set_RI, &
                                         atom_j=local_blk_rows(i_blk_local), &
                                         atom_k=local_blk_cols(j_blk_local), &
                                         i_bf_start_from_atom=i_RI_start_from_atom)
            col_local_offset = col_local_offset + col_blk_size(local_blk_cols(j_blk_local))
         END DO  
         row_local_offset = row_local_offset + row_blk_size(local_blk_rows(i_blk_local))
      END DO  
   END SUBROUTINE create_hartree_ri_3c
! **************************************************************************************************
!> \brief 
!> \note 
!> \param 
!JS First try to outsource creation of ' v_dbcsr Previously calculated 2c Coulomb repulsion between RI orbitals'
!JS needed for get_hartree_no_env as its own subroutine for rt_bse and qs_ks_method ridft
!JS parctically init_hartree from rt_bse without bs_env dependency
! **************************************************************************************************
   SUBROUTINE create_2c_coulomb_rep_v(v_dbcsr, basis_set_RI_p_type, qs_env, n_RI,&
                                     blacs_env, para_env, ri_metric, bs_env)

      TYPE(dbcsr_type)                                   :: v_dbcsr
      TYPE(dbcsr_type), ALLOCATABLE                      ::matrix
      TYPE(post_scf_bandstructure_type), POINTER, OPTIONAL         :: bs_env
      TYPE(libint_potential_type)                        :: coulomb_op
      TYPE(cp_fm_type)                                   :: v_fm, metric_fm, &
                                                            metric_inv_fm, work_fm, fm_RI_RI
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: v_dbcsr_a, &
                                                            metric_dbcsr
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(qs_environment_type), POINTER , INTENT(IN)                :: qs_env
      TYPE(section_vals_type), POINTER                   :: gw_sec, post_scf_bands_section
      REAL(KIND=dp)                                      :: eps_filter
      TYPE(libint_potential_type) , INTENT(IN)           :: ri_metric 
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: basis_set_RI_p_type(:)
      REAL(KIND=dp)                                      ::  input_regularization_RI, regularization_RI
      INTEGER, INTENT(IN)                                            :: n_RI
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_RI_global
      TYPE(mp_para_env_type), POINTER                    :: para_env

      !For dbcsr mit RI, RI dim creation without create_mat_munu Call as it is in mp2_gpw 
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist_sub
      TYPE(distribution_2d_type), POINTER                :: distribution_2d_sub
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes
      INTEGER :: blacs_grid_layout
      LOGICAL :: blacs_repeatable
      TYPE(mp_para_env_type), POINTER                    :: para_env_tensor
      INTEGER                                            :: color_sub, num_pe, natom, nkind
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_tensor
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_ri_aux
      TYPE(dbcsr_p_type)                                 :: mat_RI_RI
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      INTEGER, DIMENSION(3)                              :: periodic
      
      PRINT *, 'create_2c_coulomb_rep_v started' 
      ALLOCATE (v_dbcsr_a(1))
      ALLOCATE (metric_dbcsr(1))

      IF (.NOT. PRESENT(bs_env)) THEN
         NULLIFY (gw_sec)
         !JS Init post_scf_bands_section from qs_energy_utils.F 
         NULLIFY (post_scf_bands_section)
         post_scf_bands_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%BANDSTRUCTURE")
         gw_sec => section_vals_get_subs_vals(post_scf_bands_section, "GW")
         CALL get_qs_env(qs_env=qs_env, cell=cell)
         CALL get_cell(cell=cell, periodic=periodic, h=hmat)
         CALL section_vals_val_get(gw_sec, "REGULARIZATION_RI", r_val=input_regularization_RI)
         !Copied from gw_utils:
         IF (input_regularization_RI > -1.0E-12_dp) THEN
            regularization_RI = input_regularization_RI
         ELSE
            ! default case:
            ! 1. for periodic systems, we use the regularized resolution of the identity per default
            regularization_RI = 1.0E-2_dp
            ! 2. for molecules, no regularization is necessary
            IF (SUM(periodic) == 0) regularization_RI = 0.0_dp
         END IF

         CALL section_vals_val_get(gw_sec, "EPS_FILTER", r_val=eps_filter)

         ! ! !JS bs_env%fm_RI_RI%matrix_struct in Stepans Code needs to be created first
         NULLIFY (fm_struct_RI_global)
         CALL cp_fm_struct_create(fm_struct_RI_global, context=blacs_env, nrow_global=n_RI, &
                                    ncol_global=n_RI, para_env=para_env)
         CALL cp_fm_create(fm_RI_RI, fm_struct_RI_global)
         !CALL copy_fm_to_dbcsr(fm_RI_RI, dbcsr_ri_temp)

         !JS from mp2_gpw create_mat_munu and gw_utils:
         num_pe = para_env%num_pe 
         color_sub = para_env%mepos/num_pe

         ALLOCATE (para_env_tensor)
         CALL para_env_tensor%from_split(para_env, color_sub)

         NULLIFY (blacs_env_tensor)
         CALL cp_blacs_env_create(blacs_env=blacs_env_tensor, para_env=para_env_tensor)
         
         CALL get_qs_env(qs_env, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      molecule_set=molecule_set, &
                      molecule_kind_set=molecule_kind_set)
         NULLIFY (distribution_2d_sub)
         CALL distribute_molecules_2d(cell=cell, &
                                    atomic_kind_set=atomic_kind_set, &
                                    qs_kind_set=qs_kind_set, &
                                    particle_set=particle_set, &
                                    molecule_kind_set=molecule_kind_set, &
                                    molecule_set=molecule_set, &
                                    distribution_2d=distribution_2d_sub, &
                                    blacs_env=blacs_env_tensor, &
                                    force_env_section=qs_env%input)
         ALLOCATE (dbcsr_dist_sub)
         CALL cp_dbcsr_dist2d_to_dist(distribution_2d_sub, dbcsr_dist_sub)
         natom = SIZE(particle_set)
         ALLOCATE (row_blk_sizes(natom))
         ALLOCATE (mat_RI_RI%matrix)
         nkind = SIZE(atomic_kind_set)
         ALLOCATE (basis_set_ri_aux(nkind))
         CALL basis_set_list_setup(basis_set_ri_aux, "RI_AUX", qs_kind_set)
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=row_blk_sizes, basis=basis_set_ri_aux)
         DEALLOCATE (basis_set_ri_aux)
         CALL dbcsr_create(matrix=mat_RI_RI%matrix, &
                           name="(ai|munu)", &
                           dist=dbcsr_dist_sub, matrix_type=dbcsr_type_no_symmetry, &
                           row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes)
      ELSE 
         !JS To do. Implement params extraction from bs_env or interface
      END IF

      CALL dbcsr_create(v_dbcsr_a(1),'v_dbcsr',&
                        template = mat_RI_RI%matrix)
      CALL dbcsr_create(metric_dbcsr(1), name="RI_metric_dbcsr",&
                        template = mat_RI_RI%matrix)

      !JS From rt_bse.F hartree_init: 
      ! Calculate full coulomb RI basis elements - V _ (PQ) matrix
      NULLIFY (nl_2c)
      CALL build_2c_neighbor_lists(nl_2c, basis_set_RI_p_type, basis_set_RI_p_type, &
                                   coulomb_op, "Coulomb_neighbor_2c_list", qs_env, &
                                   sym_ij=.FALSE., molecular=.TRUE.)
      CALL build_2c_integrals(v_dbcsr_a, eps_filter, qs_env, nl_2c, &
                              basis_set_RI_p_type, basis_set_RI_p_type, coulomb_op, &
                              do_kpoints=.FALSE., regularization_RI=regularization_RI)
      ! Calculate the RI metric elements
      ! nl_2c is automatically rewritten (even reallocated) in this routine
      CALL build_2c_neighbor_lists(nl_2c, basis_set_RI_p_type, basis_set_RI_p_type, &
                                   ri_metric, "Metric_neighbor_2c_list", qs_env, &
                                   sym_ij=.FALSE., molecular=.TRUE.)
      CALL build_2c_integrals(metric_dbcsr, eps_filter, qs_env, nl_2c, &
                              basis_set_RI_p_type, basis_set_RI_p_type, ri_metric, &
                              do_kpoints=.FALSE., regularization_RI=regularization_RI)
      !JS remove 
      PRINT *, '--------------------'
      PRINT *, 'Coulomb v_dbcsr_a'
      ! !PRINT *, regularization_RI
      CALL dbcsr_print(v_dbcsr_a(1))
      ! PRINT *, '--------------------'
      ! nl_2c no longer needed
      CALL release_neighbor_list_sets(nl_2c)

      !JS create matrix, inverse matrix 
      CALL cp_fm_create(metric_fm, fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(metric_fm, 0.0_dp)
      CALL cp_fm_create(metric_inv_fm, fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(metric_inv_fm, 0.0_dp)
      CALL cp_fm_create(work_fm, fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(work_fm, 0.0_dp)

      CALL copy_dbcsr_to_fm(metric_dbcsr(1), metric_fm)
      CALL cp_fm_invert(metric_fm, metric_inv_fm)

      CALL cp_fm_create(v_fm, fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(v_fm, 0.0_dp)

      ! Multiply by the inverse from each side (M^-1 is symmetric)
      CALL cp_dbcsr_sm_fm_multiply(v_dbcsr_a(1), metric_inv_fm, &
                                    work_fm, n_RI)
      CALL parallel_gemm("N", "N", n_RI, n_RI, n_RI, &
                           1.0_dp, metric_inv_fm, work_fm, 0.0_dp, v_fm)
      ! Now, create the tensor from the matrix
      ! First, convert full matrix to dbcsr
      CALL dbcsr_clear(v_dbcsr_a(1))
      CALL copy_fm_to_dbcsr(v_fm, v_dbcsr_a(1))
      CALL dbcsr_create(v_dbcsr, "Hartree ri", v_dbcsr_a(1))
      CALL dbcsr_copy(v_dbcsr, v_dbcsr_a(1))
      ! Create and copy distinctly, so that unnecessary objects can be destroyed
      ! Destroy all unnecessary matrices
      CALL dbcsr_release(v_dbcsr_a(1))
      CALL dbcsr_release(metric_dbcsr(1))
      DEALLOCATE (v_dbcsr_a)
      DEALLOCATE (metric_dbcsr)
      CALL cp_fm_release(v_fm)
      ! CALL cp_fm_release(metric_fm(1,1))
      CALL cp_fm_release(metric_fm)
      ! DEALLOCATE(metric_fm)
      CALL cp_fm_release(work_fm)
      CALL cp_fm_release(metric_inv_fm)
      DEALLOCATE(mat_RI_RI%matrix)
      DEALLOCATE(para_env_tensor)
      DEALLOCATE(dbcsr_dist_sub)
      DEALLOCATE(row_blk_sizes)
      PRINT *, 'create_2c_coulomb_rep_v ended' 

   END SUBROUTINE create_2c_coulomb_rep_v 
END MODULE rt_bse_utils