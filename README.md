# Unrestricted Hartree-Fock (UHF) Solver with DIIS Acceleration
## Project Overview
This repository contains an implementation of the Unrestricted Hartree-Fock (UHF) method for calculating the electronic structure of molecules. Starting from a basic Restricted Hartree-Fock (RHF) skeleton, the codebase has been extensively refactored into a modular, high-performance system capable of handling open-shell systems and complex molecular geometries.
## Key Features & Contributions
* **Modular Architecture:** Refactored the original script into distinct Fortran modules to improve maintainability and scalability.
* **DIIS Convergence Accelerator:** Implemented a Direct Inversion in the Iterative Subspace (DIIS) module to significantly reduce the number of SCF iterations required for convergence.
* **UHF Formalism:** Migrated the logic to a split-spin formalism, calculating independent density and Fock matrices for $\alpha$ and $\beta$ electrons to allow for the study of radicals and ions.
* **Interactive I/O System:** Developed a dedicated I/O module that isolates terminal interactions from computational physics, featuring professional matrix formatting and user-driven job configuration.
* **External Library Integration:** Interfaces with industry-standard libraries including **LAPACK** (for diagonalization) and specialized integral libraries like **gen1int** and **interest**.

## Repository Structure
* `src/`: Contains all Fortran modules and the main program.
* `water.txt`: Example molecular input file (see also  `h2molecule.txt`, `hydroxyl.txt`).
* `documentation/`: Summary program, Detailed description subroutines and modules
* `Makefile`: Automated build system for easy compilation.

##  Installation and Usage
1.  **Compile:** Run `make` in the terminal to build the executable.
2.  **Execute:** Run `./myapp`.
3.  **Configure:** Input the molecule file (e.g., `water.txt`) and choose whether to enable the DIIS accelerator when prompted.

## Sample Output
The program provides a core Hamiltonian energies, convergence progress, final molecular orbital energies and final Hartree-Fock energy:

```text
please enter the name of the file
 ---------------------
water.txt
 you entered water.txt                                                                                  

 ---------------------
 please enter you whether you want to use the DIIS accellerator
 please enter 1 for (yes, i do want to), and 2 for (no, I don't want to)
 ---------------------
1
 you entered            1
 ---------------------
Orbital energies for the core Hamiltonian:
    -26.6612 Hartrees
     -8.8455 Hartrees
     -8.8439 Hartrees
     -8.7633 Hartrees
     -8.1705 Hartrees
     -4.9817 Hartrees
     -4.8241 Hartrees
     -4.4668 Hartrees
     -4.3660 Hartrees
     -4.3613 Hartrees
     -4.2587 Hartrees
     -4.1742 Hartrees
     -4.0000 Hartrees
     -3.9983 Hartrees
     -2.9131 Hartrees
     -2.7271 Hartrees
     -2.6050 Hartrees
     -2.5115 Hartrees
     -2.3727 Hartrees
     -2.2687 Hartrees
     -1.8108 Hartrees
     -1.8093 Hartrees
     -1.0660 Hartrees
      0.9058 Hartrees
      1.1443 Hartrees
      5.2145 Hartrees
 ---------------------
 convergence reached
 ---------------------
 Number of steps needed to converge:           16
 ---------------------
Orbital energies for the alpha orbitals after SCF:
    -15.6895 Hartrees
     -1.1420 Hartrees
     -0.7594 Hartrees
     -0.5864 Hartrees
     -0.5473 Hartrees
      0.1338 Hartrees
      0.1830 Hartrees
      0.3673 Hartrees
      0.7982 Hartrees
      0.9139 Hartrees
      0.9220 Hartrees
      0.9921 Hartrees
      1.1670 Hartrees
      1.7106 Hartrees
      2.4230 Hartrees
      2.4367 Hartrees
      2.4674 Hartrees
      2.8437 Hartrees
      3.0602 Hartrees
      5.1177 Hartrees
      5.8961 Hartrees
      6.3090 Hartrees
      6.4753 Hartrees
      6.6755 Hartrees
      6.7365 Hartrees
     14.3644 Hartrees
 ---------------------
Orbital energies for the beta orbitals after SCF:
    -15.6895 Hartrees
     -1.1420 Hartrees
     -0.7594 Hartrees
     -0.5864 Hartrees
     -0.5473 Hartrees
      0.1338 Hartrees
      0.1830 Hartrees
      0.3673 Hartrees
      0.7982 Hartrees
      0.9139 Hartrees
      0.9220 Hartrees
      0.9921 Hartrees
      1.1670 Hartrees
      1.7106 Hartrees
      2.4230 Hartrees
      2.4367 Hartrees
      2.4674 Hartrees
      2.8437 Hartrees
      3.0602 Hartrees
      5.1177 Hartrees
      5.8961 Hartrees
      6.3090 Hartrees
      6.4753 Hartrees
      6.6755 Hartrees
      6.7365 Hartrees
     14.3644 Hartrees
 ---------------------
The Hartree-Fock energy:        -63.9369 Hartrees
 ---------------------
```
