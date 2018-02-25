#include "molecule.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <functional>

typedef int (atom::*coord_func)();

int molecule::get_num_atoms() const {
	return num_atoms;
};

atom& molecule::get_atom_from_num(unsigned int num) {
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
	int num_bond_atom = 0;
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
	second_atom.set_bond_length_atom( num_bond_atom);
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
	int num_angle_atom = 0;
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
	third_atom.set_bond_length_atom(num_bond_atom);
	third_atom.set_angle_atom(num_angle_atom);
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
		remaining_atom.set_bond_length_atom(num_bond_atom);
		remaining_atom.set_angle_atom(num_angle_atom);
		remaining_atom.set_dihedral_angle_atom(num_dihedral_atom);
		atoms.push_back(remaining_atom);		
	};
	num_atoms = number_of_atoms;
};

void molecule::print_coordinates() {
	for(int i = 1; i <= num_atoms; i++) {
		if(i == 1) {
			std::cout << (*this).get_atom_from_num(i).get_name() << std::endl;
		} else if(i == 2) {
			atom& cur_atom = (*this).get_atom_from_num(i);
			std::cout << cur_atom.get_name() << " ";
			std::cout << cur_atom.get_bond_length_atom() << " ";
			std::cout << (this -> bond_length(i)) << std::endl;
		} else if(i == 3) {
			atom& cur_atom = (*this).get_atom_from_num(i);
			std::cout << cur_atom.get_name() << " ";
			std::cout << cur_atom.get_bond_length_atom() << " ";
			std::cout << (this -> bond_length(i)) << " ";
			std::cout << cur_atom.get_angle_atom() << " ";
			std::cout << (this -> angle(i)) << std::endl;;
		} else {
			atom& cur_atom = (*this).get_atom_from_num(i);
			std::cout << cur_atom.get_name() << " ";
			std::cout << cur_atom.get_bond_length_atom() << " ";
			std::cout << (this -> bond_length(i)) << " ";
			std::cout << cur_atom.get_angle_atom() << " ";
			std::cout << (this -> angle(i)) << " ";
			std::cout << cur_atom.get_dihedral_angle_atom() << " ";
			std::cout << (this -> dihedral_angle(i)) << std::endl;
		};
	};
};

void molecule::set_molecule_coord(std::string coord_file) {
	std::ifstream stream(coord_file);
	std::string line;
	int atom_counter = 0;
	while(getline(stream, line)) {
		atom_counter++;
		if(atom_counter > (*this).num_atoms) {
			std::cerr << "Too many atoms" << std::endl;
			exit(1);
		};
		std::istringstream iss(line);
		std::string name;
		iss >> name;
		if((this -> get_atom_from_num(atom_counter)).get_name() != name) {
			std::cout << "Incorrect atom name" << std::endl;
			exit(1);
		};
		double x, y, z;
		iss >> x; iss >> y; iss >> z;
		this -> set_atom_coord(atom_counter, x, y, z);
	};
	return;
};


void molecule::set_atom_coord(int atom_num, double x, double y, double z) {
	(this -> get_atom_from_num(atom_num)).set_cart_coord(x, y, z);
};	

std::vector<std::vector<double>> molecule::atom_and_connected_coord(int atom_num) {
	std::vector<coord_func> coordinate_functions = 
		{&atom::get_number, &atom::get_bond_length_atom, &atom::get_angle_atom, &atom::get_dihedral_angle_atom}; 
	std::vector<std::vector<double>> coordinates_of_atoms;
	atom& first_atom = this -> get_atom_from_num(atom_num);
	
	int num_connected = atom_num;
	if(num_connected > 4) {
		num_connected = 4;
	};	
	for(int i = 0; i < num_connected; i++) {
		coord_func current_coord_func = coordinate_functions[i];
		atom &current_atom = this -> get_atom_from_num( (first_atom.*current_coord_func)() );
		std::vector<double> vec = current_atom.get_cart_coord();
		coordinates_of_atoms.push_back(vec);
	};
	return coordinates_of_atoms;

};

double molecule::bond_length(int atom_num) {
	if(atom_num == 1) {
		std::cerr << "No bond length defined for atom one." << std::endl;
		exit(1);
	};
	/*atom& first_atom = this -> get_atom_from_num(atom_num);
	atom& second_atom = this -> get_atom_from_num(first_atom.get_bond_length_atom());
	std::vector<double> vec_1 = first_atom.get_cart_coord();
	std::vector<double> vec_2 = second_atom.get_cart_coord();
	*/
	std::vector<std::vector<double>> coordinates_of_atoms = atom_and_connected_coord(atom_num);
	return distance(coordinates_of_atoms[0], coordinates_of_atoms[1]); 
};

double molecule::angle(int atom_num) {
	if(atom_num < 3) {
		std::cerr << "No angle defined for this atom." << std::endl;
		exit(1);
	};
	/*atom& first_atom = (*this).get_atom_from_num(atom_num);	
	std::vector<double> vec_1 = first_atom.get_cart_coord();
	atom& second_atom = this -> get_atom_from_num(first_atom.get_bond_length_atom());
	std::vector<double> vec_2 = second_atom.get_cart_coord();	
	atom& third_atom = this -> get_atom_from_num(first_atom.get_angle_atom());
	std::vector<double> vec_3 = third_atom.get_cart_coord();
	*/
	std::vector<std::vector<double>> coordinates_of_atoms = atom_and_connected_coord(atom_num);
	return calculate_angle(coordinates_of_atoms[0], coordinates_of_atoms[1], coordinates_of_atoms[2]); 
};

double molecule::dihedral_angle(int atom_num) {
	if(atom_num < 4) {
		std::cerr << "No dihedral angle defined for this atom." << std::endl;
		exit(1);
	};
	/*atom& first_atom = (*this).get_atom_from_num(atom_num);
	std::vector<double> vec_1 = first_atom.get_cart_coord();
	atom& second_atom = this -> get_atom_from_num(first_atom.get_bond_length_atom());
	std::vector<double> vec_2 = second_atom.get_cart_coord();	
	atom& third_atom = this -> get_atom_from_num(first_atom.get_angle_atom());
	std::vector<double> vec_3 = third_atom.get_cart_coord();
	atom& fourth_atom = this -> get_atom_from_num(first_atom.get_dihedral_angle_atom());
	std::vector<double> vec_4 = fourth_atom.get_cart_coord();
	*/
	std::vector<std::vector<double>> coordinates_of_atoms = atom_and_connected_coord(atom_num);
	return calculate_dihedral_angle(coordinates_of_atoms[0], coordinates_of_atoms[1], coordinates_of_atoms[2],
			coordinates_of_atoms[3]); 
};

double molecule::dot_product(std::vector<double> coord_1, std::vector<double> coord_2) {
	if(!check_coord_size(coord_1) || !check_coord_size(coord_2)) {
		coord_length_error_message();
		exit(1);	
	};	
	return coord_1[0]*coord_2[0] + 
		coord_1[1]*coord_2[1] + 
		coord_1[2]*coord_2[2];
};

std::vector<double> molecule::cross_product(std::vector<double> coord_1, std::vector<double> coord_2) {
	if(!check_coord_size(coord_1) || !check_coord_size(coord_2)) {
		coord_length_error_message();
		exit(1);	
	}
	std::vector<double> result = {coord_1[1]*coord_2[2] - coord_1[2]*coord_2[1], 
		coord_1[2]*coord_2[0] - coord_1[0]*coord_2[2], coord_1[0]*coord_2[1] - coord_2[0]*coord_1[1]};
	return result;
};
std::vector<double> molecule::displacement_vector(std::vector<double> coord_1, std::vector<double> coord_2) {
	if(!check_coord_size(coord_1) || !check_coord_size(coord_2)) {
		coord_length_error_message();
		exit(1);	
	};
	std::vector<double> displacement = {coord_1[0] -coord_2[0], coord_1[1] - coord_2[1], coord_1[2] -coord_2[2]};
	return displacement;
};

double molecule::calculate_angle(std::vector<double> coord_1, std::vector<double> coord_2, std::vector<double> coord_3) {
	std::vector<double> vec_1 = displacement_vector(coord_1, coord_2);
	std::vector<double> vec_2 = displacement_vector(coord_3, coord_2);
	return calculate_vec_angle(vec_1, vec_2);	
};

double molecule::calculate_vec_angle(std::vector<double> vec_1, std::vector<double> vec_2) {
	double numerator = dot_product(vec_1, vec_2);
	double denominator = vec_distance(vec_1)*vec_distance(vec_2);
	double angle = acos(numerator/denominator);
	return angle;
};

double molecule::vec_distance(std::vector<double> vec) {
	return sqrt(dot_product(vec, vec));
};

std::vector<double> molecule::vec_normalised(std::vector<double> vec) {
	std::vector<double> n_vec;
	double length = vec_distance(vec);
	for(int i = 0; i < vec.size(); i++) {
		n_vec.push_back(vec[i]/length);
	};	
	return n_vec;
};

double molecule::calculate_dihedral_angle(std::vector<double> coord_1, std::vector<double> coord_2,
	std::vector<double> coord_3, std::vector<double> coord_4) {
		
	std::vector<double> vec_1 = displacement_vector(coord_1, coord_2);
	std::vector<double> vec_2 = displacement_vector(coord_3, coord_2);

	std::vector<double> vec_3 = displacement_vector(coord_4, coord_3);

	std::vector<double> new_vec_1 = vec_normalised(cross_product(vec_1, vec_2));
	std::vector<double> new_vec_2 = vec_normalised(cross_product(vec_2, vec_3));
	std::vector<double> new_vec_3 = vec_normalised(vec_2);
	std::vector<double> new_vec_4 = cross_product(new_vec_1, new_vec_3);

	double x = dot_product(new_vec_1, new_vec_2), y = dot_product(new_vec_4, new_vec_2);
	
	return atan2(y,x);
};
 
double molecule::distance(std::vector<double> coord_1, std::vector<double> coord_2) {
	std::vector<double> displacement_vec = displacement_vector(coord_1, coord_2);
	return sqrt(dot_product(displacement_vec, displacement_vec));
};

bool molecule::check_coord_size(std::vector<double> coord) {
	if(coord.size() !=3) {
		return false;
	} else {
		return true;
	};
};

void molecule::file_error_message() {
	std::cerr << "File contains an error" << std::endl;	
};

void molecule::coord_length_error_message() {
	std::cerr << "Coordinates incorrect length" << std::endl;
};

