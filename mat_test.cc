#include "matrix_op.h"
#include <iostream>

int main() {
	Eigen::MatrixXd mat_1(4,3);
	Eigen::MatrixXd mat_2(2,2);
	mat_1 << 1,2,4,
		1,3,9,
		1,4,16,
		1,5,25;
	mat_2 << 1,5,
		1,6; 
	Eigen::MatrixXd mat_3 = outer_product(mat_1, mat_2);
	std::cout << mat_1 << std::endl << mat_2 << std::endl << mat_3 << std::endl;
	std::vector<int> vec = {4,3,0,4,5,6,4,3,3};
	sort_vec(vec);
	for(int a : vec) {
		std::cout << a << " ";
	};	
	std::cout << std::endl;

	Eigen::MatrixXd mat_1i = mat_1.completeOrthogonalDecomposition().pseudoInverse();
	Eigen::MatrixXd mat_2i = mat_2.completeOrthogonalDecomposition().pseudoInverse();
	std::cout << mat_1i << std::endl;
	std::cout << " right \n " << mat_1.transpose()*(mat_1*mat_1.transpose()).inverse() << std::endl;
	std::cout << " left \n " << (mat_1.transpose()*mat_1).inverse()*mat_1.transpose() << std::endl;

	std::cout << mat_1i*mat_1 << std::endl;
	std::cout << mat_2i*mat_2 << std::endl;
	std::cout << outer_product(mat_2i, mat_1i)*mat_3;
	int order =  2;
	int num_grid = 3;
	for(int i = 0; i < 9; i++) {
		std::cout << number_convert(i, order, num_grid) << std::endl;
	};
};
