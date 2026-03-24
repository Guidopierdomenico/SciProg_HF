**Detailed description of modules and subroutines**



define\_molecule subroutine (present in main program):

define\_molecule takes the file that's been inputted and reads the number of atoms, the charges of the atoms and the coordinates of the atoms and uses this as the input of the add\_atoms\_to\_molecule subroutine which updates the molecule.



define\_basis subroutine (present in main program):

define\_basis takes as input the molecule and the number of atoms and places basis sets on the atoms using the add\_shell\_to\_basis subroutine:
For hydrogen 3 s functions with exponents 3, 1 and 0.1.

For every other atom 5 s functions with exponents (0.1, 0.35, 1, 3, 10), 3 p functions (with exponents 0.2, 1, 5)

and (1 d function with exponent 1).





Molecular structure module:

molecular structure centers around the derived type "molecular\_structure\_t", which contains the number of atoms, the charge and the coordinates of the molecule.

It then contains the add\_atoms\_to\_molecule subroutine, which takes the molecule that needs to be updated and the atoms that need to be added to the molecule, it then computes new number of atoms of the molecule and also adds the charge and coordinates of those molecules to the original molecule.



ao\_basis module:

ao\_basis centers aronund two derived types, which are dependent on each other: basis\_func\_info\_t, which stores the physical properties of a single basis function shell(such as its center coordinates, orbital momentum, exponents and contraction coefficients), and basis\_set\_info\_t which holds an array of the individual shells to represent the entire molecular basis set.

Then the add\_shell\_to\_basis subroutine takes in the angular momentum number, the exponents and coefficients of the basis functions and adds them to the basis set.



Compute\_integrals module:

Compute integrals acts as a wrapper for exteranl integral libraries. It contains compute\_1e\_integrals(which interfaces with gen1int to calculate Overla, Kinetic and Potential matrices) and generate\_2int(which interfaces with interest to calculate the 4D array of two-electron integrals).



Diagonalization module:

Diagonalization contains two subroutines: solve\_genev and diagonalize.

Solve\_genev solves the generalized eigenvalue problem using Lodwins transformation to orthonormal basis.

Diagonalize on the other hand diagonalizes the matrix that is inputted using Lapack's DSYEV routine.















