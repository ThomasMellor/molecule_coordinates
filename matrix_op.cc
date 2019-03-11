#include "matrix_op.h"
#include <bits/stdc++.h>
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

std::vector<int> sort_vec(std::vector<int>& vec) {
	sort(vec.begin(), vec.end());	
	return vec;
};
bool contains_num(const std::vector<int>& vec, int x) {
	return (std::find(vec.begin(), vec.end(), x) != vec.end());
};

bool contains_all_nums(const std::vector<int>& vec, const  std::vector<int>& sub_vec) {
	for(int x : sub_vec) {
		if(!contains_num(vec, x)) {
			return false;
		};
	};
	return true; 
};
