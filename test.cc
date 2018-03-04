#include "atom.h"
#include <iostream>
#include "molecule.h"
#include "eigen-eigen-5a0156e40feb/Eigen/Dense"

int main() {
	molecule mol("./test_2.txt", "CH4");
	mol.set_molecule_coord("./test_coords_2.txt");
	mol.print_coordinates();
	std::vector<std::vector<double>> test_coord= {{},{1},{1,2}};
	Eigen::MatrixXd mat = mol.empty_matrix();
	Eigen::VectorXd vec = mol.empty_vector();
	Eigen::VectorXd x = mol.empty_vector();
	mol.derivative_vector(vec, test_coord);
	std::cout << vec << std::endl;
	mol.derivative_matrix(mat);
	std::cout << mat << std::endl;	
	std::cout << mol.bond_length_derivative(2,2, "x") << std::endl;
	std::cout << mol.angle_derivative(3,2, "x") << std::endl;
	std::cout << mol.bond_length_derivative(3,2,"x") << std::endl;
	int counter = 0;
//while((mol.coord_difference(test_coord) > 0.001) || (counter < 10000)) {
//		mol.derivative_matrix(mat);
		//mol.derivative_vector(vec, test_coord);
		x = mat.colPivHouseholderQr().solve(vec);
		std::cout << x << std::endl;
		mol.update_molecule_coord(x);	
		//std::cout << counter << " " << mol.coord_difference(test_coord) << std::endl;
		//mol.print_coordinates();
	//};
	mol.print_coordinates();
};
