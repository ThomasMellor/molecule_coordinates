#include "matrix_op.h"

Eigen::MatrixXd outer_product(const Eigen::MatrixXd& mat_a, const Eigen::MatrixXd& mat_b) {
	
	int row_a = mat_a.rows(); int row_b = mat_b.rows(); 
	int col_a = mat_a.cols(); int col_b = mat_b.rows();
	Eigen::MatrixXd mat_ab = Eigen::MatrixXd::Zero(row_a*row_b, col_a*col_b); 
	
	for(int i = 0; i < row_a; i++ ) {
		for(int j = 0; j < col_a; j++) {
			mat_ab.block(i*row_b, j*col_b, row_b, col_b) = mat_b*mat_a(i,j);
		};
	};
	return mat_ab;
};
