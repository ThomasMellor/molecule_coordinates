#ifndef MATRIX_OP_H
#define MATRIX_OP_H


#include "eigen-eigen-5a0156e40feb/Eigen/Dense"
#include "eigen-eigen-5a0156e40feb/Eigen/QR"
#include <algorithm> 

Eigen::MatrixXd outer_product(const Eigen::MatrixXd& mat_a, const Eigen::MatrixXd& mat_b);
std::vector<int> sort_vec(std::vector<int>& vec); 



bool contains_num(const std::vector<int>& vec, int x);
bool contains_all_nums(const std::vector<int>& vec, const std::vector<int>& sub_vec); 

#endif 

