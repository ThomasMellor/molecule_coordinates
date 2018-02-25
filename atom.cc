#include "atom.h"

atom::atom(unsigned int num, std::string nm) : number(num), name(nm) {};

std::string atom::get_name() const {
	return name;
};	

int atom::get_number() {
	return number;
};

std::vector<double> atom::get_cart_coord() const {
	return cart_coord;	
};

void atom::set_cart_coord(double x, double y, double z) {
	cart_coord = {x, y, z};
	return;
};

void atom::set_bond_length_atom(int connected_to) {
	bond_length_atom = connected_to;
};

void atom::set_angle_atom(int atom_1) {
	angle_atom = atom_1; 
};

void atom::set_dihedral_angle_atom( int atom_1) {
	dihedral_angle_atom = atom_1; 
};

int atom::get_bond_length_atom() {
	return bond_length_atom;
};
int atom::get_angle_atom() {
	return angle_atom;
};

int atom::get_dihedral_angle_atom() {
	return dihedral_angle_atom;
};
