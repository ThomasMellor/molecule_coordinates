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

void atom::set_bond_length_atom(atom const& connected_to) {
	bond_length_atom = &connected_to;
};

void atom::set_angle_atoms(atom const& atom_1, atom const& atom_2) {
	angle_atoms = {&atom_1, &atom_2};
};

void atom::set_dihedral_angle_atoms(atom const& atom_1, atom const& atom_2, atom const& atom_3) {
	dihedral_angle_atoms = {&atom_1, &atom_2, &atom_3};
};
