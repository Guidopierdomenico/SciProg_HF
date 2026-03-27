module core_hamiltonian
    use compute_integrals
    use diagonalization
    use molecular_structure
    use ao_basis
    implicit none
    private
    public construct_core_hamiltonian

contains

    subroutine construct_core_hamiltonian(molecule, ao_basis, n_AO, core_hamiltonian, S, C, eps)
        ! Variable containing the molecular structure
        type(molecular_structure_t), intent(in) :: molecule
        ! Variable containing the atomic orbital basis
        type(basis_set_info_t), intent(in) :: ao_basis
        !number atomic orbitals
        integer, intent(in)  :: n_AO
        !passed in from main already allocated
        real(8), intent(inout) :: core_hamiltonian(:,:), S(:,:), C(:,:), eps(:)
        !local arrays
        real(8) :: V(n_AO,n_AO),T(n_AO,n_AO)

        ! Compute the overlap matrix
        call   compute_1e_integrals ("OVL",ao_basis,ao_basis,S)
        ! Compute the kinetic matrix
        call   compute_1e_integrals ("KIN",ao_basis,ao_basis,T)
        ! Compute the potential matrix
        call   compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)
        ! Compute the core Hamiltonian matrix (the potential is positive, we scale with -e = -1 to get to the potential energy matrix)
        core_hamiltonian = T - V
        ! Diagonalize the core hamiltonian matrix
        call solve_genev (core_hamiltonian,S,C,eps)
        
    end subroutine

end module