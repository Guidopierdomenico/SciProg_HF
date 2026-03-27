module DIIS
    use molecular_structure
    use Diagonalization
    implicit none
    private
    public DIIS_loop

contains

    subroutine DIIS_loop(molecule, n_AO, number_atoms, n_alpha, n_beta, core_hamiltonian, S, C, ao_integrals)
        !Arguments passed in from Main)
        type(molecular_structure_t), intent(in) :: molecule
        integer, intent(in)                     :: n_AO, number_atoms, n_alpha, n_beta
        real(8), intent(in)                     :: core_hamiltonian(n_AO,n_AO), S(n_AO,n_AO)
        real(8), intent(in)                     :: C(n_AO,n_AO)
        real(8), intent(in)                     :: ao_integrals(n_AO,n_AO,n_AO,n_AO)

        !DIIS Parameter
        integer, parameter :: max_DIIS = 10

        !Standard Local Memory
        real(8) :: F_alpha(n_AO,n_AO), F_beta(n_AO,n_AO)
        real(8) :: C_alpha(n_AO,n_AO), C_beta(n_AO,n_AO)
        real(8) :: D_alpha(n_AO,n_AO), D_beta(n_AO,n_AO), D_total(n_AO,n_AO)
        real(8) :: D_alpha_previous(n_AO,n_AO), D_beta_previous(n_AO,n_AO)
        real(8) :: eps_alpha(n_AO), eps_beta(n_AO)
        
         !IIS-Specific Memory (Automatically sized using n_AO and max_DIIS)
        real(8) :: F_alpha_new(n_AO,n_AO), F_beta_new(n_AO,n_AO)
        real(8) :: error_matrix_alpha(n_AO,n_AO), error_matrix_beta(n_AO,n_AO)
        real(8) :: F_alpha_memory(n_AO,n_AO, max_DIIS), F_beta_memory(n_AO,n_AO, max_DIIS)
        real(8) :: error_matrix_alpha_memory(n_AO,n_AO, max_DIIS), error_matrix_beta_memory(n_AO,n_AO, max_DIIS)
        real(8) :: B_DIIS(max_DIIS, max_DIIS), A_DIIS(max_DIIS+1, max_DIIS+1)
        real(8) :: rhs_DIIS(max_DIIS+1), coeffs_DIIS(max_DIIS+1)

        ! Local Variables
        integer :: number_iterations_scf_loop, maximum_number_iterations
        integer :: kappa, lambda, mu, nu, i, j, m
        real(8) :: E_HF, convergence_treshold, convergence_alpha, convergence_beta
        real(8) :: distance_ij, nuclear_repulsion

        !initializing D_previous matrix to 0 for the first iteration
        D_alpha_previous = 0.0_8
        D_beta_previous = 0.0_8

        !setting the maximum number of iterations
        maximum_number_iterations = 100

        !setting initial value of coeficcients matrices
        C_alpha = C
        C_beta = C
        
        !setting max amount of matrices to be held by DIIS memory arrays
        max_DIIS = 10

        do
            number_iterations_scf_loop = number_iterations_scf_loop + 1
            !setting the loop variable for the DIIS loops
            m = min(number_iterations_scf_loop, max_DIIS)
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

            !calculating error matrix to check for DIIS convergence
            error_matrix_alpha = matmul(matmul(F_alpha, D_alpha), S) - matmul(matmul(S, D_alpha), F_alpha)
            error_matrix_beta = matmul(matmul(F_beta, D_beta), S) - matmul(matmul(S, D_beta), F_beta) 
            !creating a memory of the SCF loop matrices to be used for the DIIS routine
            !filling first all of the available slots
            if (number_iterations_scf_loop <= max_DIIS) then 
                !Fock matrices memory
                F_alpha_memory(:,:,number_iterations_scf_loop) = F_alpha
                F_beta_memory(:,:,number_iterations_scf_loop) = F_beta
                !Error matrices memory
                error_matrix_alpha_memory(:,:,number_iterations_scf_loop) = error_matrix_alpha
                error_matrix_beta_memory(:,:,number_iterations_scf_loop) = error_matrix_beta
                !when the memory array is full move everything down by one slot, making space for the matrix of the current iteration
            else
                !Fock matrices memory
                F_alpha_memory(:,:,1:max_DIIS-1) = F_alpha_memory(:,:,2:max_DIIS)
                F_alpha_memory(:,:,max_DIIS) = F_alpha
                F_beta_memory(:,:,1:max_DIIS-1) = F_beta_memory(:,:,2:max_DIIS)
                F_beta_memory(:,:,max_DIIS) = F_beta
                !Error matrices memory
                error_matrix_alpha_memory(:,:,1:max_DIIS-1) = error_matrix_alpha_memory(:,:,2:max_DIIS)
                error_matrix_alpha_memory(:,:,max_DIIS) = error_matrix_alpha
                error_matrix_beta_memory(:,:,1:max_DIIS-1) = error_matrix_beta_memory(:,:,2:max_DIIS)
                error_matrix_beta_memory(:,:,max_DIIS) = error_matrix_beta
            endif 

            !Constructing B matrix for DIIS routine
            !loop until the current iteration until memory has been filled
            do j = 1, m
                do i = 1, m
                    B_DIIS(i,j) = sum(error_matrix_alpha_memory(:,:,i) * error_matrix_alpha_memory(:,:,j)) + sum(error_matrix_beta_memory(:,:,i) * error_matrix_beta_memory(:,:,j))
                enddo
            enddo
            !Constructing A matrix for DIIS routine
            A_DIIS = 0.0_8
            A_DIIS(:m, :m) = B_DIIS(:m,:m)
            A_DIIS(:m, m+1) = -1.0_8
            A_DIIS(m+1, :m) = -1.0_8 
            !Constructing rhs vector for DIIS routine
            rhs_DIIS = 0.0_8
            rhs_DIIS(m+1) = -1.0_8

            !Solviing the system of equations to find tge DIIS coefficients
            call linearsolver(A_DIIS(:m+1, :m+1), rhs_DIIS(1:m+1), coeffs_DIIS(1:m+1))

            !Calculate a new Fock matrix using the DIIS coefficients and the Fock matrix memory
            F_alpha_new = 0.0_8
            do i = 1, m
                F_alpha_new = F_alpha_new + coeffs_DIIS(i)*F_alpha_memory(:,:, i)
            enddo
            F_beta_new = 0.0_8
            do i = 1, m
                F_beta_new = F_beta_new + coeffs_DIIS(i)*F_beta_memory(:,:, i)
            enddo

        

            !exit loop if convergence is lower than treshold
            if (( maxval(abs(error_matrix_alpha)) < convergence_treshold) .and. ( maxval(abs(error_matrix_beta)) < convergence_treshold)) then
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
            call solve_genev (F_alpha_new,S,C_alpha,eps_alpha)
            call solve_genev (F_beta_new,S,C_beta,eps_beta)


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