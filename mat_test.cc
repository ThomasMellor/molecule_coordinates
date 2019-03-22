#include "matrix_op.h"
#include <iostream>

int main() {
	Eigen::MatrixXd mat_1(6,3);
	Eigen::MatrixXd mat_2(6,3);
	Eigen::MatrixXd mat_3(6,3);
	mat_1 << 1,2,4,
		1,4,16,
		1,5,15,
		3,42,54,
		4,423,2,
		4,99,8;
		
	mat_2 << 1,6,36,
		1,7,49,
		1,8,64,
		4,4,53,
		4,9,8,
		10,9,8;
	

	mat_3 << 1,9,81,
		1,10,100,
		1,11,121,
		9,8,9,
		10,9,89,
		988,7,6;
	
	Eigen::MatrixXd vec_1(27,1);
		
	for(int i = 0; i < vec_1.rows(); i++) {
		vec_1(i,0) = i;
	};
	std::cout << " vec 1 " << vec_1 << std::endl;
	Eigen::MatrixXd vec_2 = vec_1;
	vec_2.resize(3,9);
	std::cout << " vec 2 " << vec_2 << std::endl;
	std::cout << "first " << outer_product(outer_product(mat_1,mat_2),mat_3)*vec_1 << std::endl;
	Eigen::MatrixXd vec_3 = (mat_3*vec_2*outer_product(mat_1.transpose(), mat_2.transpose()));
	std::cout << "vec 3 " << vec_3 << std::endl;
	Eigen::MatrixXd vec_4 = vec_3;
	vec_4.resize(6*6*6,1);
	std::cout << "rows " <<  vec_4.rows() << std::endl;
	std::cout << "second "  << vec_4 << std::endl;
};
