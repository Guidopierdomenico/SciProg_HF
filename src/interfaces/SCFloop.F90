module SCF_loop
    use molecular_structure
    use Diagonalization
    implicit none
    private
    public Calculate_Hartree_Fock_Energy

contains
    subroutine Calculate_Hartree_Fock_Energy(molecule, n_AO, number_atoms, n_alpha, n_beta, core_hamiltonian, S, C, ao_integrals)
        ! Variable containing the molecular structure
        type(molecular_structure_t), intent(in) :: molecule
        ! Variable naming as in the description of the exercise
        real(8), intent(in) :: core_hamiltonian(:,:), S(:,:), C(:,:)
        integer, intent(in)  :: n_AO
        integer  :: kappa, lambda, mu, nu
        real(8)  :: E_HF
        real(8), intent(in) :: ao_integrals (:,:,:,:)
        !variable for number of alpha orbitals and number of beta orbitals
        integer, intent(in) :: n_alpha, n_beta
        !variable to keep track of the number of iterations of the scf loop
        integer :: number_iterations_scf_loop
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
        !arrays to hold Fock matrices for both alpha spin and beta spin
        real(8)  :: F_alpha(n_AO,n_AO), F_beta(n_AO,n_AO)
        !arrays to hold coefficients and orbitals energies for both alpha spin and beta spin
        real(8)  :: C_alpha(n_AO,n_AO), C_beta(n_AO,n_AO), eps_alpha(n_AO), eps_beta(n_AO)
        !array to hold density matrix calculated at the previous iteration, used to check for convergence
        real(8)  :: D_alpha_previous(n_AO,n_AO), D_beta_previous(n_AO,n_AO)
        !array to hold density matrices for alpha spin and for beta spin and total density matrix
        real(8)  :: D_alpha(n_AO,n_AO), D_beta(n_AO,n_AO), D_total(n_AO,n_AO)
        !variable for the maximum amount of iterations
        integer :: maximum_number_iterations

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

        !printing alpha and beta orbital energies
        print '("Orbital energies for the alpha orbitals:")'
        print '((F12.4), " Hartrees")', eps_alpha
        print *, "---------------------"
        print '("Orbital energies for the beta orbitals:")'
        print '((F12.4), " Hartrees")', eps_beta
        print *, "---------------------"
        !printing hartree fock energy
        print '("The Hartree-Fock energy:    ", (F12.4), " Hartrees")', E_HF 
        print *, "---------------------"
    end subroutine

end module