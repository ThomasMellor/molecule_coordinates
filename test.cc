#include "atom.h"
#include <iostream>
#include "molecule_2.h"
#include "eigen-eigen-5a0156e40feb/Eigen/Dense"

int main() {
	molecule mol("./test.txt", "CH4");
	mol.set_molecule_coord(0 ,"./h20eq.txt");
	mol.set_molecule_coord_Z(1,"./test_coords_2.txt");
	mol.print_cart_coords(0);
//	mol.print_coordinates(0);
//	mol.print_coordinates(1);
	mol.print_cart_coords(1);
	mol.set_coefficients("./H2O.1.testC.out");
	mol.set_L_matrix("./h20_l.txt");
//	std::cout << mol.normal_coordinates() << std::endl;
	
	std::cout << "energy " << mol.energy() << std::endl;	
	mol.print_cart_coords(1);
};
