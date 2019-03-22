#include "atom.h"
#include <iostream>
#include "molecule_2.h"
#include "eigen-eigen-5a0156e40feb/Eigen/Dense"

int main() {
	molecule mol("./h2co_Z.txt", "H2CO");
	mol.set_molecule_coord(0,"./h2co_eq.txt");
	//mol.set_molecule_coord_Z(1, "./h20_coord.txt");
	std::cout << "test 1" << std::endl;
	mol.set_L_matrix("./h2co-4D-nc.pot");
	std::cout << " test 2 " << std::endl;
	mol.set_grid_coeffs("./h2co-4D-nc.pot", 6);
	std::cout << " test 3 " << std::endl;
	//std::cout << mol.energy() << std::endl; 
	mol.print_cart_coords(0);
	mol.print_coordinates(0);
};
