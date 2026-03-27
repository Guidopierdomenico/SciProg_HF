module InputOutput

    implicit none
    private
    public Input_file, Print_core_hamiltoian_orbital_energies, Print_orbital_energies_after_SCF, Print_HFenergy

contains 

    subroutine Input_file(filename, DIIS_use)
        !variable to hold the filename
        character(100), intent(out) :: filename
        !variable to hold the answer for DIIS use
        integer, intent(out) :: DIIS_use

        !asking user which file they want to analyze
        print *, "please enter the name of the file"
        print *, "---------------------"
        read *, filename
        print *, "you entered ", filename
        print *, "---------------------"

        !asking user whether they want to use DIIS or not
        print *, "please enter you whether want to use the DIIS accellerator"
        print *, "please enter 1 for (yes, i do want to), and 2 for (no, I don't want to)"
        print *, "---------------------"
        read *, DIIS_use
        print *, "you entered ", DIIS_use
        print *, "---------------------"
    end subroutine

    subroutine Print_core_hamiltoian_orbital_energies(eps)
        real(8), intent(in) :: eps(:)

        print '("Orbital energies for the core Hamiltonian:")'
        print '((F12.4), " Hartrees")', eps
        print *, "---------------------"
    end subroutine

    subroutine Print_orbital_energies_after_SCF(eps_alpha, eps_beta)
        real(8), intent(in) :: eps_alpha(:), eps_beta(:)

        !printing alpha and beta orbital energies
        print '("Orbital energies for the alpha orbitals after SCF:")'
        print '((F12.4), " Hartrees")', eps_alpha
        print *, "---------------------"
        print '("Orbital energies for the beta orbitals after SCF:")'
        print '((F12.4), " Hartrees")', eps_beta
        print *, "---------------------"
    end subroutine

    subroutine Print_HFenergy(E_HF)
        real(8), intent(in) :: E_HF
        !printing hartree fock energy
        print '("The Hartree-Fock energy:    ", (F12.4), " Hartrees")', E_HF 
        print *, "---------------------"
    end subroutine

end module