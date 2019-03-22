#include "atom.h"
#include <iostream>
#include "molecule_2.h"
#include "eigen-eigen-5a0156e40feb/Eigen/Dense"

int main() {
	molecule mol("./h20_Z.txt", "H2O");
	mol.set_molecule_coord(0,"./h2o_eq.txt");
	mol.set_molecule_coord_Z(1, "./h20_coord.txt");
	mol.set_L_matrix("./h2o_1i.xpot");
	mol.set_molecule_coord_norm("./h2o_norm.txt"); 
	mol.set_grid_coeffs("./h2o_1i.xpot", 6);
	std::cout << mol.energy() << std::endl; 
	//mol.print_cart_coords(0);
	//mol.print_cart_coords(1);
};
