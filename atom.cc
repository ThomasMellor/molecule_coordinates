#include "atom.h"
#include <limits>
#include <iostream>

atom::atom(int num, std::string nm) :  name(nm) , bond_length_atom(-1),
	angle_atom(-1), dihedral_angle_atom(-1), number(check_sign(num)) {};

atom::atom(int num, std::string nm, int second_atom) :  name(nm), number(check_sign(num)),
	bond_length_atom(check_sign(second_atom)), angle_atom(-1), dihedral_angle_atom(-1) {};

atom::atom(int num, std::string nm, int second_atom,int third_atom) :  name(nm), number(num),
	bond_length_atom(check_sign(second_atom)), angle_atom(check_sign(third_atom)), dihedral_angle_atom(-1) {};
		
atom::atom(int num, std::string nm, int second_atom,
	int third_atom, int fourth_atom) : name(nm), number(num),
	bond_length_atom(check_sign(second_atom)), angle_atom(check_sign(third_atom)),
   	dihedral_angle_atom(check_sign(fourth_atom)) {};		

std::string atom::get_name() const {
	return name;
};	

int atom::get_number() const {
	return number;
};

std::vector<double> atom::get_cart_coord() const {
	return cart_coord;	
};

void atom::set_cart_coord(double x, double y, double z) {
	cart_coord = {x, y, z};
	return;
};

int atom::get_bond_length_atom() const {
	return check_value(bond_length_atom);
};

int atom::get_angle_atom() const{
	return check_value(angle_atom);
};

int atom::get_dihedral_angle_atom() const {
	return check_value(dihedral_angle_atom);
};

std::vector<int> atom::get_connected_atoms() const {
	std::vector<int> atoms;
	atoms.push_back(get_number());
	if(bond_length_atom != -1) {
		atoms.push_back(bond_length_atom);
	};
	if(angle_atom != -1) {
		atoms.push_back(angle_atom);
	};
	if(dihedral_angle_atom != -1) {
		atoms.push_back(dihedral_angle_atom); 
	};
	return atoms;
};

int atom::check_sign(int atom_num) {
	if(atom_num < 0) {
		std::cerr << "Atom number invalid" << std::endl;
		exit(1);
	} else {
		return atom_num;
	};
};

int atom::check_value(int atom_num) {
	if(atom_num != -1) {
		return atom_num;
	} else {
		std::cerr << "Atom does not have this connection" << std::endl;
		exit(1);
	};
};


