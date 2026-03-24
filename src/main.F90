program HartreeFock

  ! Guido Pierdomenico, march 2026

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
  real(8), allocatable :: core_hamiltonian(:,:), V(:,:),T(:,:),S(:,:), C(:,:), eps(:)





  ! The following large array can be eliminated when Fock matrix contruction is implemented
  real(8), allocatable :: ao_integrals (:,:,:,:)




  !variable to keep track of the number of iterations of the scf loop
  integer :: number_iterations_scf_loop
  !variable to temporarily store the calculated energy so that it can be confronted with the following iteration to check for convergenve
  real(8) :: previous_energy
  !variable to set the convergence treshold
  real(8) :: convergence_treshold
  !variables to keep track of convergence for both Density matrices
  real(8) :: convergence_alpha, convergence_beta
  !variable to calculate distance between a pair of atoms
  real(8) :: distance_ij
  !variable to calculate nuclear repulsion
  real(8) :: nuclear_repulsion
  !variable for the number of atoms in the molecule
  integer :: number_atoms
  !loop variables
  integer :: i, j
  !variable to hold the filename
  character(100) :: filename
  !arrays to hold Fock matrices for both alpha spin and beta spin
  real(8), allocatable  :: F_alpha(:,:), F_beta(:,:)
  !arrays to hold coefficients and orbitals energies for both alpha spin and beta spin
  real(8), allocatable  :: C_alpha(:,:), C_beta(:,:), eps_alpha(:), eps_beta(:)
  !array to hold density matrix calculated at the previous iteration, used to check for convergence
  real(8), allocatable  :: D_alpha_previous(:,:), D_beta_previous(:,:)
  !array to hold density matrices for alpha spin and for beta spin and total density matrix
  real(8), allocatable  :: D_alpha(:,:), D_beta(:,:), D_total(:,:)
  !variable for number of alpha orbitals and number of beta orbitals
  integer :: n_alpha, n_beta, n_total
  !variable for the maximum amount of iterations
  integer :: maximum_number_iterations
  

  print *, "please enter the name of the file"
  print *, "---------------------"
  read *, filename
  print *, "you entered ", filename
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
  print '((F12.4))', eps
  print *, "---------------------"

  !initializing variables for scf loop
  convergence_treshold = 1.0d-9
  number_iterations_scf_loop = 0
  previous_energy = 0.0_8

  ! Compute all 2-electron integrals
  allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
  call generate_2int (ao_basis,ao_integrals)

  !allocating memory for arrays used in SCF loop
  allocate (D_alpha(n_AO,n_AO))
  allocate (D_beta(n_AO,n_AO))
  allocate (D_total(n_AO,n_AO))
  allocate (D_alpha_previous(n_AO,n_AO))
  allocate (D_beta_previous(n_AO,n_AO))
  allocate (C_alpha(n_AO,n_AO))
  allocate (C_beta(n_AO,n_AO))
  allocate (eps_alpha(n_AO))
  allocate (eps_beta(n_AO))
  !initializing D_previous matrix to 0 for the first iteration
  D_alpha_previous = 0.0_8
  D_beta_previous = 0.0_8

  !setting the maximum number of iterations
  maximum_number_iterations = 100

  !setting initial value of coeficcients matrices
  C_alpha = C
  C_beta = C
  !SCF loop
  do
  number_iterations_scf_loop = number_iterations_scf_loop + 1
  
    !if statement to exit loop if the maximum number of iterations has been exceded
    if (number_iterations_scf_loop > maximum_number_iterations) then
      print *, "error: maximum number of iterations reached"
      exit
    endif

  ! Form the density matrix
    do lambda = 1, n_AO
        do kappa = 1, n_AO
          D_alpha(kappa,lambda) = sum(C_alpha(kappa,1:n_alpha)*C_alpha(lambda,1:n_alpha))
          D_beta(kappa,lambda) = sum(C_beta(kappa,1:n_beta)*C_beta(lambda,1:n_beta))
      end do
    end do
    D_total = D_alpha + D_beta
    !Construct  the Fock matrix
    F_alpha = core_hamiltonian
    F_beta = core_hamiltonian
    do nu = 1, n_AO
      do mu = 1, n_AO
        F_alpha = F_alpha + (ao_integrals(:,:,mu,nu))* D_total(mu, nu) - ao_integrals(:,nu,mu,:)*D_alpha(mu, nu)
        F_beta = F_beta + (ao_integrals(:,:,mu,nu))* D_total(mu, nu) - ao_integrals(:,nu,mu,:)*D_beta(mu, nu)
      enddo
    enddo
    ! Compute the Hartree-Fock energy
    E_HF = 0.5D0 * sum((core_hamiltonian + F_alpha) * D_alpha) + 0.5D0 * sum((core_hamiltonian + F_beta) * D_beta)
    !calculate convergence
    convergence_alpha = sqrt( sum ( abs( (D_alpha-D_alpha_previous)**2)))
    convergence_beta = sqrt( sum ( abs( (D_beta-D_beta_previous)**2)))
    !exit loop if convergence is lower than treshold
    if ((abs(convergence_alpha) < convergence_treshold) .and. (abs(convergence_beta) < convergence_treshold)) then
      print *, "convergence reached"
      print *, "---------------------"
      print *, "Number of steps needed to converge: ", number_iterations_scf_loop
      print *, "---------------------"
      exit
    endif



    !update previous density
    D_alpha_previous = D_alpha
    D_beta_previous = D_beta
    ! Diagonalize the Fock matrix to get new coefficients
    call solve_genev (F_alpha,S,C_alpha,eps_alpha)
    call solve_genev (F_beta,S,C_beta,eps_beta)


  enddo


  !compute nuclear repulsion
  nuclear_repulsion = 0.0_8
  do j = 2, number_atoms
    do i = 1, j-1
      !calculate distance between atom i and atom j
      distance_ij = sqrt( sum( (molecule%coord(:,i)-molecule%coord(:,j))**2 ) )
      !add repulsion energy of the pair to the total
      nuclear_repulsion = nuclear_repulsion + (molecule%charge(i)*molecule%charge(j))/distance_ij
    enddo
  enddo
  !add nuclear repulsion
  E_HF = E_HF + nuclear_repulsion

  !printing hartree fock energy
  print '("The Hartree-Fock energy:    ", (F12.4), " Hartrees")', E_HF 
  print *, "---------------------"
  

end program


  subroutine define_molecule(molecule, filename, number_atoms)
    use molecular_structure
    type(molecular_structure_t), intent(inout) :: molecule
    character(len=*), intent(in)               :: filename
    integer, intent(out)                       :: number_atoms
    real(8), allocatable                       :: charge(:),coord(:,:)
    integer                                    :: i

    !reading the number of atoms from the input file
    open(10, file = filename, status='old')
    read(10, *) number_atoms
    !allocating charges and coordinates arrays using number of atoms just read
    allocate(charge(number_atoms))
    allocate(coord(3, number_atoms))
    !loop to read the charges and coordinates from the file
    do i = 1, number_atoms
      read(10, *) charge(i), coord(:,i)
    enddo
    close(10)

    call add_atoms_to_molecule(molecule,charge,coord)
  end subroutine

  subroutine define_basis(ao_basis, molecule, number_atoms) 
    use ao_basis
    use molecular_structure
    type(basis_set_info_t), intent(inout)   :: ao_basis
    type(basis_func_info_t)                 :: gto
    integer, intent(in)                     :: number_atoms
    type(molecular_structure_t), intent(in) :: molecule
    integer                                 :: i
    !do loop to place the basis function on molecule coordinates
    do i = 1, number_atoms
      !for hydrogen 3 s functions
      if (molecule%charge(i) == 1) then
        call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),3.D0)
        call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),1.D0)  
        call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),1.D-1)
      else
      !for all non-hydrogen elements 5 s functions, 3 p functions and 1 d function
        !s functions
        call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),1.D-1)
        call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),35.D-2)
        call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),1.D0)
        call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),3.D0)
        call add_shell_to_basis(ao_basis,0,molecule%coord(:,i),1.D1)
        !p functions
        call add_shell_to_basis(ao_basis,1,molecule%coord(:,i),2.D-1)
        call add_shell_to_basis(ao_basis,1,molecule%coord(:,i),1.D0)
        call add_shell_to_basis(ao_basis,1,molecule%coord(:,i),5.D0)
        !d function
        call add_shell_to_basis(ao_basis,2,molecule%coord(:,i),1.D0)
      endif
    enddo
  end subroutine

   

