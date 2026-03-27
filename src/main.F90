program HartreeFock

  ! Guido Pierdomenico, march 2026
  use molecular_structure
  use ao_basis
  use compute_integrals
  use diagonalization
  use LinearSystem
  use SCF_loop
  use DIIS
  implicit none

  ! Variable containing the molecular structure
  type(molecular_structure_t) :: molecule
  ! Variable containing the atomic orbital basis
  type(basis_set_info_t) :: ao_basis
  ! Variable naming as in the description of the exercise
  integer  :: n_AO, n_occ
  integer  :: kappa, lambda, mu, nu
  real(8)  :: E_HF
  real(8), allocatable :: core_hamiltonian(:,:), V(:,:),T(:,:),S(:,:), C(:,:), eps(:)
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


  !asking user which file they want to analyze
  print *, "please enter the name of the file"
  print *, "---------------------"
  read *, filename
  print *, "you entered ", filename
  print *, "---------------------"

  !asking user whether they want to use DIIS or not
  print *, "please enter you want to use the DIIS accellerator"
  print *, "---------------------"
  print *, "please enter 1 for yes, i do want to and 2 for no, I don't want to"
  print *, "---------------------"
  read *, DIIS_use
  print *, "you entered ", DIIS_use
  print *, "---------------------"

  ! Definition of the molecule
  call define_molecule(molecule, filename, number_atoms)
  ! Definition of the GTOs
  call define_basis(ao_basis, molecule, number_atoms)
  n_AO = ao_basis%nao

  ! Definition of the number of  number of spin up and spin down orbitals
  n_total = nint(sum(molecule%charge))
  n_beta = n_total/2
  n_alpha = n_total - n_beta

  ! Compute the overlap matrix
  allocate (S(n_AO,n_AO))
  call   compute_1e_integrals ("OVL",ao_basis,ao_basis,S)
  ! Compute the kinetic matrix
  allocate (T(n_AO,n_AO))
  call   compute_1e_integrals ("KIN",ao_basis,ao_basis,T)
  ! Compute the potential matrix
  allocate (V(n_AO,n_AO))
  call   compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)
  ! Compute the core Hamiltonian matrix (the potential is positive, we scale with -e = -1 to get to the potential energy matrix)
  allocate(core_hamiltonian(n_AO,n_AO))
  core_hamiltonian = T - V
  ! Diagonalize the core hamiltonian matrix
  allocate (C(n_AO,n_AO))
  allocate (eps(n_AO))
  call solve_genev (core_hamiltonian,S,C,eps)
  print '("Orbital energies for the core Hamiltonian:")'
  print '((F12.4), " Hartrees")', eps
  print *, "---------------------"


  ! Compute all 2-electron integrals
  allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
  call generate_2int (ao_basis,ao_integrals)

  if (DIIS_use == 1) then
    call DIIS_loop(molecule, n_AO, number_atoms, n_alpha, n_beta, core_hamiltonian, S, C, ao_integrals)
  elseif (DIIS_use == 2) then
    call Calculate_Hartree_Fock_Energy(molecule, n_AO, number_atoms, n_alpha, n_beta, core_hamiltonian, S, C, ao_integrals)
  else
    print '("you have entered the wrong number for the DIIS choice")'
  endif
  

  
end program