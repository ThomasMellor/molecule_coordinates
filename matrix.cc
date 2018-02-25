#include "matrix.h"
#include <iostream>

unsigned int matrix::size() const {
	return dimension;
}

int matrix::point_1D(int i, int j, int N) {
	return N*i+j;
};

double matrix::point(int x, int y) const {
	return (*this).matrix_points[point_1D(x, y, dimension)]; 
};

void matrix::set(int x, int y, double val) {
	if( x < 0 || x >= dimension || y < 0 || y >= dimension) {
		std::cerr << x << ", " << y << " out of matrix bounds." << std::endl;
		exit(1);
	};	
	this->matrix_points[point_1D(x, y, dimension)] = val; 
};

/*
 * Constructors for matrices
 */
matrix::matrix(unsigned int s) : dimension(s) {};

matrix::matrix(const matrix& mat) : dimension(mat.size()) {
	    this-> matrix_points = mat.matrix_points;
};

matrix& matrix::operator=(const matrix& mat) {
	if( this -> size() != mat.size()) {
		throw std::invalid_argument("matrices not the same size");
	};
        this -> matrix_points = mat.matrix_points;
        return *this;
};

matrix matrix::invert(const matrix& mat) const {
	int N = mat.size()
	mat inverted_mat(N);
	for(int p = 0; p < N; p++) {
		if(mat
	};
};

