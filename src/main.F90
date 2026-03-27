program HartreeFock

  ! Guido Pierdomenico, march 2026
  use molecular_structure
  use ao_basis
  use compute_integrals
  use diagonalization
  use LinearSystem
  use SCF_loop
  use DIIS
  use core_hamiltonian
  use InputOutput
  implicit none

  ! Variable containing the molecular structure
  type(molecular_structure_t) :: molecule
  ! Variable containing the atomic orbital basis
  type(basis_set_info_t) :: ao_basis
  ! Variable naming as in the description of the exercise
  integer  :: n_AO, n_occ
  integer  :: kappa, lambda, mu, nu
  real(8)  :: E_HF
  real(8), allocatable :: core_hamiltonian(:,:), S(:,:), C(:,:), eps(:)
  !array to hold atomic orbitals integrals
  real(8), allocatable :: ao_integrals (:,:,:,:)
  !variable for the number of atoms in the molecule
  integer :: number_atoms
  !variable to hold the filename
  character(100) :: filename
  !variable to hold the answer for DIIS use
  integer :: DIIS_use
  !variable for number of alpha orbitals and number of beta orbitals
  integer :: n_alpha, n_beta, n_total
  !variable for orbital energies
  real(8), allocatable :: eps_alpha(:), eps_beta(:)


  !asking user which file they want to analyze and whether they want to use DIIS or not
  call Input_file(filename, DIIS_use)
  

  ! Definition of the molecule
  call define_molecule(molecule, filename, number_atoms)
  ! Definition of the GTOs
  call define_basis(ao_basis, molecule, number_atoms)
  n_AO = ao_basis%nao
  ! Definition of the number of  number of spin up and spin down orbitals
  n_total = nint(sum(molecule%charge))
  n_beta = n_total/2
  n_alpha = n_total - n_beta

  !allocate arrays
  allocate(core_hamiltonian(n_AO,n_AO))
  allocate(S(n_AO,n_AO))
  allocate(C(n_AO,n_AO))
  allocate(eps((n_AO)))

  !construct core hamiltonian
  call construct_core_hamiltonian(molecule, ao_basis, n_AO, core_hamiltonian, S, C, eps)
  !print orbital energies for the core hamiltonian
  call Print_core_hamiltoian_orbital_energies(eps)

  ! Compute all 2-electron integrals
  allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
  call generate_2int (ao_basis,ao_integrals)

  !allocating orbital energies for later printing
  allocate(eps_alpha(n_AO))
  allocate(eps_beta(n_AO))

  !calculate hartree fock energy using SCF loop, with or without DIIS accelerator depending on user's choice
  if (DIIS_use == 1) then
    call DIIS_loop(molecule, n_AO, number_atoms, n_alpha, n_beta, core_hamiltonian, S, C, ao_integrals, E_HF, eps_alpha, eps_beta)
  elseif (DIIS_use == 2) then
    call Calculate_Hartree_Fock_Energy(molecule, n_AO, number_atoms, n_alpha, n_beta, core_hamiltonian, S, C, ao_integrals, E_HF, eps_alpha, eps_beta)
  else
    print '("you have entered the wrong number for the DIIS choice")'
  endif
  
  !print orbital energies and hartree fock eneregy
  call Print_orbital_energies_after_SCF(eps_alpha, eps_beta)
  call Print_HFenergy(E_HF)
end program