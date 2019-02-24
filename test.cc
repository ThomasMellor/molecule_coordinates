#include "atom.h"
#include <iostream>
#include "molecule_2.h"
#include "eigen-eigen-5a0156e40feb/Eigen/Dense"

int main() {
	molecule mol("./h20_Z.txt", "H2O");
	mol.set_molecule_coord(0,"./h2o_eq.txt");
	mol.set_molecule_coord(1, "./h20_coord.txt");
	mol.set_L_matrix("./h2o.xpot");
	mol.print_coordinates(0);
	mol.print_coordinates(1);
	std::cout << mol.L_mat << std::endl;
	std::cout << mol.M_mat << std::endl;
	Eigen::MatrixXd A = mol.Amat();
	std::cout << A << std::endl;
	Eigen::MatrixXd T = molecule::Tmat(A);
	mol.rotate_coords(T);
	std::cout << "Eckart" << std::endl << mol.Eckart_cond() <<std::endl;
	std::cout << "normal" << std::endl << mol.normal_coordinates()  << std::endl;
};
