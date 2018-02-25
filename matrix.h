#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include "math.h"

class matrix {
	private:
		const unsigned int dimension;
		std::vector<double> matrix_points = std::vector<double>(pow(dimension, 2)); //array of points
		static int point_1D(int, int, int);
	public:
		matrix(unsigned int);
		matrix(const matrix&);
		matrix& operator=(const matrix&);
		~matrix() {};

		double point(int, int) const;
		void set(int, int, double);
		unsigned int size() const;
		matrix invert(const matrix) const;		
};

#endif
