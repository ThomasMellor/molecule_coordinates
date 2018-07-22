#include "atom.h"
#include <iostream>
#include "molecule_2.h"
#include "eigen-eigen-5a0156e40feb/Eigen/Dense"

int main() {
	molecule mol("./test.txt", "CH4");
	mol.set_molecule_coord(0,"./test_coords_1.txt");	
 	mol.set_molecule_coord_Z(1,"./test_coords_2.txt");

		mol.set_L_matrix("glyoxal_d2.pot");

//	std::cout << mol.get_atom_from_num(4).get_dihedral_angle_atom();	
//	std::cout << " hello " << mol.L_mat << std::endl;
//	std::cout << "yo" << mol.M_mat << std::endl;

//	mol.print_cart_coords(0);
	mol.print_coordinates(1);

	//mol.print_cart_coords(1);
	//mol.print_coordinates(1);
//	std::cout << "com " << mol.centre_of_mass(0) << std::endl; 	
//	std::cout << "com 2 " << mol.centre_of_mass(1) << std::endl;
	Eigen::MatrixXd A = mol.Amat();
//	std::cout << std::endl;
//	std::cout << A << std::endl;
//	std::cout << "test" << std::endl;
	Eigen::MatrixXd ATA = molecule::ATAmat(A);
	Eigen::MatrixXd T = molecule::Tmat(A,ATA);	
	mol.rotate_coords(T);
	
	std::cout << "eckart" << mol.Eckart_cond() << std::endl;
	
	
//	mol.print_cart_coords(0);
//	mol.print_cart_coords(1);
	std::cout << "yo" << std::endl;
	std::cout << mol.normal_coordinates() << std::endl;
	/*	std::vector<std::vector<double>> test_coord= {{},{1},{1,1}, {1,2,0.1},{1.4,2,3},{1,2,0.1},{1,2,2},{2,2,-2}};
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
	*/
};
