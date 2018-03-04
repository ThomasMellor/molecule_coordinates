#include "atom.h"
#include <iostream>
#include "molecule.h"
#include "eigen-eigen-5a0156e40feb/Eigen/Dense"

int main() {
	molecule mol("./test.txt", "CH4");
	mol.set_molecule_coord("./test_coords_1.txt");
	mol.print_coordinates();
	std::vector<std::vector<double>> test_coord= {{},{1},{1,1}, {1,2,0.1},{1.4,2,3},{1,2,0.1},{1,2,2},{2,2,-2}};
	Eigen::MatrixXd mat = mol.empty_matrix();
	Eigen::VectorXd vec = mol.empty_vector();
	Eigen::VectorXd x = mol.empty_vector();
	mol.derivative_vector(vec, test_coord);
	mol.derivative_matrix(mat);
	int counter = 0;
while((mol.coord_difference(test_coord) > 0.00000001) || (counter < 5)) {
		mol.derivative_matrix(mat);
		mol.derivative_vector(vec, test_coord);
		x = mat.colPivHouseholderQr().solve(vec);
		x =-x;
		counter++;
		mol.update_molecule_coord(x);	
		std::cout << counter << " " << mol.coord_difference(test_coord) << std::endl;
		mol.print_coordinates();
	};
	mol.print_coordinates(); 
};
