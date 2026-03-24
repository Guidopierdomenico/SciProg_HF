**Summary program**

The program takes as input a text file that the user wants to analyse, containing the description of a molecule: the number of atoms, the charges of the atoms and the coordinates of these atoms.



The program then reads the file and places basis sets on the atoms:
For hydrogen 3 s functions with exponents 3, 1 and 0.1.

For every other atom 5 s functions with exponents (0.1, 0.35, 1, 3, 10), 3 p functions (with exponents 0.2, 1, 5)

and (1 d function with exponent 1).



Then the overlap, the kinetic and the potential matrix are computed using the atomic orbitals basis sets.

Then the core Hamiltonian is computed as the difference of the kinetic and potential energy matrices.

This is then Diagonalized to get the orbital energies.



Then the program start the self consistent field (SCF) loop to get the Hartree-Fock energy for the given molecule.

The SCF loop contains the following steps:

* takes count of the number of iterations
* checks if number of iterations has exceeded the threshold, in which case it exits the loop and prints an error message
* Forms the density matrix
* Constructs the Fock matrix
* Computes the Hartree-Fock energy
* Calculates the convergence: if the convergence is reached it exits the loop  and prints that the convergence was reached and how many iterations it took to reach it; if it wasn't it goes on with the next steps
* updates previous density, which will be used to check for convergence
* diagonalizes the Fock matrix to get new coefficients
* Starts new iteration



After the SCF loop, the nuclear repulsion is computed and added to the energy calculated in the loop to get the complete Hartree-Fock energy, which is then printed.

