!This module defines a molecular structure data type
       module molecular_structure

        implicit none
        private

        public add_atoms_to_molecule, define_molecule
        integer, parameter, private:: REALD=8 ! default real for this module

!       Molecule:
        type, public:: molecular_structure_t
         integer:: num_atoms=0                    ! number of atoms
         real(REALD), allocatable :: charge(:)    ! nuclear charges
         real(REALD), allocatable :: coord(:,:)   ! coordinates of the atoms
        end type molecular_structure_t

        contains

        subroutine add_atoms_to_molecule (molecule,add_charge,add_coord)
          type(molecular_structure_t), intent(inout) :: molecule   ! the molecule to be updated
          real(REALD), intent(in) :: add_charge(:), add_coord(:,:) ! the additional atoms
 
          real(REALD), allocatable :: charge_new(:)
          real(REALD), allocatable :: coord_new(:,:)
          integer :: num_atoms_new, add_atoms
          
          add_atoms = size(add_charge)
          if (size(add_coord,1) /= 3)         stop "first coordinate dimension should be 3"
          if (size(add_coord,2) /= add_atoms) stop "inconsistency in charges and coordinates"

          ! Make the extended molecule
          num_atoms_new = molecule%num_atoms + add_atoms
          allocate(charge_new(num_atoms_new)) 
          allocate(coord_new(3,num_atoms_new)) 
          charge_new (molecule%num_atoms+1:num_atoms_new) = add_charge
          coord_new(:,molecule%num_atoms+1:num_atoms_new) = add_coord
          if (molecule%num_atoms /= 0) then
             charge_new(1:molecule%num_atoms)  = molecule%charge
             coord_new(:,1:molecule%num_atoms) = molecule%coord
             deallocate (molecule%charge)
             deallocate (molecule%coord)
          end if

          ! and copy it into the molecule type
          allocate (molecule%charge(num_atoms_new))
          allocate (molecule%coord(3,num_atoms_new))
          molecule%charge = charge_new
          molecule%coord  = coord_new
          molecule%num_atoms = num_atoms_new

          ! local allocatable arrays should be destroyed automatically, but just to make sure
          deallocate (charge_new)
          deallocate (coord_new)

        end subroutine

        subroutine define_molecule(molecule, filename, number_atoms)
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
  
end module molecular_structure

       