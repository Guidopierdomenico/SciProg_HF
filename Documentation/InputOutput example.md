**Input/Output example**



As example for input file I created two text files, one for an hydrogen molecule "h2molecule.txt" and one for a water molecule "water.txt" (I have put both of these files in the general folder of the program so that after cloning my repository the program can be immediately run without needing to change anything).

The first line of both files contains the number of atoms of the molecule while the following lines contains the charge of the atom as the first value and the XYZ coordinates in Bohr as the following three values in the line.



Example input file "water.txt":

3

8.0 0.00000 0.00000 0.22000

1.0 0.00000 1.43000 -0.88000

1.0 0.00000 -1.43000 -0.88000



In order to run the program first the "make" command needs to be given to the terminal after which the program will be compiled and can be run using the "./myapp" command in the terminal.

When the program is run the terminal will show a "please enter the name of the file" prompt, after which the full name of the file needs to be inputted(e.g. water.txt, in the case of the water molecule). After this the program will print the name of the file that has been inputted and will print the orbital energies, whether convergence has been reached and the steps needed to reach convergence and the HartreeFock energy of the molecule.



Example Terminal output:

please enter the name of the file

&#x20;---------------------

water.txt

&#x20;you entered water.txt

&#x20;---------------------

&#x20;Orbital energies for the core Hamiltonian:

&#x20;    -32.0000    -1.5000    -1.2000    -0.8000 ... \[etc]

&#x20;---------------------

&#x20;convergence reached

&#x20;---------------------

&#x20;Number of steps needed to converge: 15

&#x20;---------------------

&#x20;The Hartree-Fock energy:     -63.9000 Hartrees

&#x20;---------------------

