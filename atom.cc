#include "atom.h"

atom::atom(unsigned int num, std::string nm) : number(num), name(nm) {};

std::string atom::get_name() const {
	return name;
};	

int atom::get_number() const {
	return number;
};

std::vector<double> atom::get_cart_coords() const {
	return cart_coords;	
};

void atom::set_cart_coords(double x, double y, double z) {
	cart_coords = {x, y, z};
	return;
};
