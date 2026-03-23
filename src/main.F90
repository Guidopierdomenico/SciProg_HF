program HartreeFock

  ! Demonstration program that can be used as a starting point
  ! Lucas Visscher, March 2022

  use molecular_structure
  use ao_basis
  use compute_integrals
  use diagonalization

  implicit none

  ! Variable containing the molecular structure
  type(molecular_structure_t) :: molecule
  ! Variable containing the atomic orbital basis
  type(basis_set_info_t) :: ao_basis


  ! Variable naming as in the description of the exercise
  integer  :: n_AO, n_occ
  integer  :: kappa, lambda, mu, nu
  real(8)  :: E_HF
  real(8), allocatable :: core_hamiltonian(:,:), F(:,:),V(:,:),T(:,:),S(:,:), C(:,:), eps(:), D(:,:)





  ! The following large array can be eliminated when Fock matrix contruction is implemented
  real(8), allocatable :: ao_integrals (:,:,:,:)




  !variable to keep track of the number of iterations of the scf loop
  integer :: number_iterations_scf_loop
  !variable to temporarily store the calculated energy so that it can be confronted with the following iteration to check for convergenve
  real(8) :: previous_energy
  !variable to set the convergence treshold
  real(8) :: convergence_treshold
  !variable to keep track of convergence
  real(8) :: convergence
  !variable to calculate distance between a pair of atoms
  real(8) :: distance_ij
  !variable to calculate nuclear repulsion
  real(8) :: nuclear_repulsion
  !variable for the number of atoms in the molecule
  integer :: n_atoms
  !loop variables
  integer :: i, j


  ! Definition of the molecule
  call define_molecule(molecule)

  ! Definition of the GTOs
  call define_basis(ao_basis)
  n_AO = ao_basis%nao

  ! Definition of the number of occupied orbitals
  n_occ = 3 ! hardwired for this demonstration program, should be set via input

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
  print*, "Orbital energies for the core Hamiltonian:",eps


  !initializing variables for scf loop
  convergence_treshold = 1.0d-9
  number_iterations_scf_loop = 0
  previous_energy = 0.0_8

  ! Compute all 2-electron integrals
  allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
  call generate_2int (ao_basis,ao_integrals)

  !allocating memory for Fock matrix and Density matrix
  allocate (F(n_AO,n_AO))
  allocate (D(n_AO,n_AO))


  !SCF loop
  do

  number_iterations_scf_loop = number_iterations_scf_loop + 1
  ! Form the density matrix
  do lambda = 1, n_AO
      do kappa = 1, n_AO
        D(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))
    end do
  end do

  !Construct  the Fock matrix
  F = core_hamiltonian
  do nu = 1, n_AO
    do mu = 1, n_AO
      F = F + (2.D0 * ao_integrals(:,:,mu,nu) - ao_integrals(:,nu,mu,:))* D(mu, nu)
    enddo
  enddo

  ! Compute the Hartree-Fock energy
  E_HF = sum((core_hamiltonian+F) * D)
 

  !calculate convergence
  convergence = E_HF - previous_energy
  !exit loop if convergence is lower than treshold
  if (abs(convergence) < convergence_treshold) exit
  !update previous energy
  previous_energy = E_HF

  ! Diagonalize the Fock matrix to get new coefficients
  call solve_genev (F,S,C,eps)

    

  enddo

  n_atoms = 2

  !add nuclear repulsion
  nuclear_repulsion = 0.0_8
  do j = 1 + 1, n_atoms
    do i = 1, n_atoms - 1
      !calculate distance between atom i and atom j
      distance_ij = sqrt( sum( (molecule%coord(:,i)-molecule%coord(:,j))**2 ) )
      !add repulsion energy of the pair to the total
      nuclear_repulsion = nuclear_repulsion + (molecule%charge(i)*molecule%charge(j))/distance_ij
    enddo
  enddo
  E_HF = E_HF + nuclear_repulsion

  print*, "The Hartree-Fock energy:    ", E_HF
  !print '((F12.4))', 

end


   subroutine define_molecule(molecule)
     ! This routine should be improved such that an arbitrary molecule can be given as input
     ! the coordinates below are for a be-he dimer oriented along the x-axis with a bond length of 2 au
     use molecular_structure
     type(molecular_structure_t), intent(inout) :: molecule
     real(8) :: charge(2),coord(3,2)
     charge(1)   = 4.D0
     charge(2)   = 2.D0
     coord       = 0.D0
     coord(1,2)  = 2.D0
     call add_atoms_to_molecule(molecule,charge,coord)
   end subroutine

   subroutine define_basis(ao_basis)
    ! This routine can be extended to use better basis sets 
    ! The coordinates of the shell centers are the nuclear coordinates
    ! Think of a refactoring of define_molecule and define_basis to ensure consistency 
     use ao_basis
     type(basis_set_info_t), intent(inout) :: ao_basis
     type(basis_func_info_t) :: gto
     ! Be:  2 uncontracted s-funs:    l      coord          exp      
     call add_shell_to_basis(ao_basis,0,(/0.D0,0.D0,0.D0/),4.D0)
     call add_shell_to_basis(ao_basis,0,(/0.D0,0.D0,0.D0/),1.D0)
     ! He:  1 uncontracted s-fun:     l      coord          exp      
     call add_shell_to_basis(ao_basis,0,(/2.D0,0.D0,0.D0/),1.D0)
   end subroutine

   

