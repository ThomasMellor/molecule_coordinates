#include "matrix_op.h"
#include <iostream>


int main() {
	Eigen::MatrixXd mat_1(3,4);
	Eigen::MatrixXd mat_2(2,2);
	mat_1 << 3,5,2,5,4,4,5,5,4,4,6,4;
	mat_2 << 3,8,7,3; 
	Eigen::MatrixXd mat_3 = outer_product(mat_1, mat_2);
	std::cout << mat_1 << std::endl << mat_2 << std::endl << mat_3 << std::endl;
};
