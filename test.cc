#include "atom.h"
#include <iostream>
#include "molecule.h"

int main() {
	molecule mol("./test.txt", "CH4");
	mol.set_molecule_coord("./test_coords_1.txt");
	mol.print_coordinates();
	std::cout << mol.bond_length_derivative(2, 1, "x");
};
