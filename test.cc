#include "atom.h"
#include <iostream>
#include "molecule.h"
#include "eigen-eigen-5a0156e40feb/Eigen/"

int main() {
	molecule mol("./test.txt", "CH4");
	mol.set_molecule_coord("./test_coords_1.txt");
	mol.print_coordinates();
	std::cout << mol.dihedral_angle_derivative(6,4, "y");
};
