#include "molecule.h"
#include <fstream>
#include <sstream>
#include <iostream>

int molecule::get_num_atoms() {
	return num_atoms;
};

atom const& molecule::get_atom_from_num(unsigned int num) {
	return atoms[num-1];	
};

molecule::molecule(std::string z_matrix_file, std::string molecule_name) : name(molecule_name) {
	std::ifstream stream(z_matrix_file);
	if(!stream) {
		std::cerr << "Error opening file " + z_matrix_file << std::endl;
		exit(1);
	};
	std::string line;

	if(!getline(stream, line)){
		file_error_message();
		exit(1);	
	};
	std::istringstream iss_1(line);
	std::string atom_name;
	int number_of_atoms = 0;
	if(!(iss_1 >> atom_name)){
		file_error_message();
		exit(1);
	};
	number_of_atoms++;
	atoms.push_back(atom(number_of_atoms, atom_name));
	
	if(!getline(stream, line)) {
		file_error_message();
		exit(1);
	};
	std::istringstream iss_2(line);
	if(!(iss_2 >> atom_name)){
		file_error_message();
		exit(1);
	};
	int num_bond_atom;
	if(!(iss_2 >> num_bond_atom)){
		file_error_message();
		exit(1);
	};
	if(num_bond_atom > number_of_atoms) {
		file_error_message();
		exit(1);
	};
	number_of_atoms++;
	atom second_atom(number_of_atoms, atom_name);
	second_atom.set_bond_length_atom((*this).get_atom_from_num(num_bond_atom));
	atoms.push_back(second_atom);

	if(!getline(stream, line)) {
		return;
	};
	std::istringstream iss_3(line);
	if(!(iss_3 >> atom_name)){
		file_error_message();
		exit(1);
	};
	if(!(iss_3 >> num_bond_atom)){
		file_error_message();
		exit(1);
	};
	if(num_bond_atom > number_of_atoms) {
		file_error_message();
		exit(1);
	};
	int num_angle_atom;
	if(!(iss_3 >> num_angle_atom)){
		file_error_message();
		exit(1);
	};
	if(num_angle_atom > number_of_atoms) {
		file_error_message();
		exit(1);
	};
	number_of_atoms++;
	atom third_atom(number_of_atoms, atom_name);
	third_atom.set_bond_length_atom(this-> get_atom_from_num(num_bond_atom));
	third_atom.set_angle_atoms(this -> get_atom_from_num(num_bond_atom), 
			this-> get_atom_from_num(num_angle_atom));
	atoms.push_back(third_atom);

	while(getline(stream, line)) {
		std::istringstream iss_4(line);
		if(!(iss_4 >> atom_name)){
			file_error_message();
			exit(1);
		};
		if(!(iss_4 >> num_bond_atom)){
			file_error_message();
			exit(1);
		};
		if(num_bond_atom > number_of_atoms) {
			file_error_message();
			exit(1);
		};
		if(!(iss_4 >> num_angle_atom)){
			file_error_message();
			exit(1);
		};
		if(num_angle_atom > number_of_atoms) {
			file_error_message();
			exit(1);
		};	
		int num_dihedral_atom;
		if(!(iss_4 >> num_dihedral_atom)){
			file_error_message();
			exit(1);
		};
		number_of_atoms++;
		atom remaining_atom(number_of_atoms, atom_name);
		remaining_atom.set_bond_length_atom(this -> get_atom_from_num(num_bond_atom));
		remaining_atom.set_angle_atoms(this -> get_atom_from_num(num_bond_atom),
			this -> get_atom_from_num(num_angle_atom));
		remaining_atom.set_dihedral_angle_atoms( this -> get_atom_from_num(num_angle_atom),
			this -> get_atom_from_num(num_bond_atom),
			this -> get_atom_from_num(num_dihedral_atom));
		atoms.push_back(remaining_atom);		
	};
};

void file_error_message() {
	std::cerr << "File contains an error" << std::endl;	
};


