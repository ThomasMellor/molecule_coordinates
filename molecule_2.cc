#include "molecule_2.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <functional>
#include <map>



typedef int (atom::*coord_func)() const;

/*Eigen::MatrixXd molecule::empty_matrix() {
	return Eigen::MatrixXd(3*num_atoms, 3*num_atoms); 
};

Eigen::MatrixXd& molecule::derivative_matrix(Eigen::MatrixXd& mat) {	
	std::vector<std::string> axes = {"x", "y", "z"};
	for(int i = 1; i <= num_atoms; i++) {
		for(std::string axis_i : axes) {
			for(int j = 1; j <= i; j++) {
				for(std::string axis_j : axes) {
					double val = 0;
					for(int k = 2; k <= num_atoms; k++) {
						atom& first_atom = this -> get_atom_from_num(k);
						std::vector<int> con_atoms = {first_atom.get_number(), first_atom.get_bond_length_atom()};
						
						if( check_derivative_atoms(con_atoms, i, j)){
							val += bond_length_derivative(k, i, axis_i)*bond_length_derivative(k, j, axis_j);
						};
					};
					for(int k = 3; k <= num_atoms; k++) {
						atom& first_atom = this -> get_atom_from_num(k);
						std::vector<int> con_atoms = {first_atom.get_number(), first_atom.get_bond_length_atom(), 
							first_atom.get_angle_atom()};
						if(check_derivative_atoms(con_atoms, i, j)){
							val += angle_derivative(k, i, axis_i)*angle_derivative(k, j, axis_j);
						};
					};
					for(int k = 4; k <= num_atoms; k++) {
						atom& first_atom = this -> get_atom_from_num(k);
						std::vector<int> con_atoms = first_atom.get_connected_atoms();
						if(check_derivative_atoms(con_atoms, i, j)){
							val += dihedral_angle_derivative(k, i, axis_i)*dihedral_angle_derivative(k, j, axis_j);
						};
					};
					mat(3*(i-1) + axes_name_to_num(axis_i), 3*(j-1) + axes_name_to_num(axis_j)) = val;
					mat(3*(j-1) + axes_name_to_num(axis_j), 3*(i-1) + axes_name_to_num(axis_i)) = val;
				};
			};
		};	
	};
	return mat;
};

bool molecule::check_derivative_atoms(std::vector<int> connected_atoms, int atom_1, int atom_2) {
	bool bool_a1 = false;
	bool bool_a2 = false;
	for(int con_atom : connected_atoms) {
		if(atom_1 == con_atom) {
			bool_a1 = true;
		};
		if(atom_2 == con_atom) {
			bool_a2 = true;
		};	
	};
	return (bool_a1 & bool_a2); 
};

Eigen::VectorXd& molecule::derivative_vector(Eigen::VectorXd& vec, std::vector<std::vector<double>> try_coord) {
	std::vector<std::string> axes = {"x", "y", "z"};
	for(int i = 1; i <= num_atoms; i++) {
		for(std::string axis : axes) {
			double val = 0;	
			for(int k = 2; k <= num_atoms; k++) {
				atom& first_atom = this -> get_atom_from_num(k);
				std::vector<int> con_atoms = {first_atom.get_number(), first_atom.get_bond_length_atom()};

				if( check_derivative_atoms(con_atoms, i, i)){
					val += bond_length_derivative(k, i, axis)*(bond_length(k) - try_coord[k-1][0]);
				};
			};
			for(int k = 3; k <= num_atoms; k++) {
				atom& first_atom = this -> get_atom_from_num(k);
				std::vector<int> con_atoms = {first_atom.get_number(), first_atom.get_bond_length_atom(), 
					first_atom.get_angle_atom()};
				if(check_derivative_atoms(con_atoms, i, i)){
				val += angle_derivative(k, i, axis)*(angle(k) - try_coord[k-1][1]);
				};
			};
			for(int k = 4; k <= num_atoms; k++) {
				atom& first_atom = this -> get_atom_from_num(k);
				std::vector<int> con_atoms = first_atom.get_connected_atoms();
				if(check_derivative_atoms(con_atoms, i, i)){
					val += dihedral_angle_derivative(k, i, axis)*(dihedral_angle(k) - try_coord[k-1][2]);
				};
			};
			vec(3*(i-1) + axes_name_to_num(axis)) = val;
		};
	};
	return vec;	
};

Eigen::VectorXd molecule::empty_vector() {
	return Eigen::VectorXd(3*num_atoms);
};

double molecule::coord_difference(std::vector<std::vector<double>> try_coord) {
	double val = 0;
	for(int i = 2; i <= num_atoms; i++) {
		val += pow(bond_length(i) - try_coord[i-1][0], 2);
	};
	for(int i = 3; i <= num_atoms; i++) {
		val += pow(angle(i) - try_coord[i-1][1],2);
	};
	for(int i = 4; i <= num_atoms; i++) {
		val += pow(dihedral_angle(i) - try_coord[i-1][2],2);
	};
	return val;
};
*/
int molecule::get_num_atoms() const {
	return num_atoms;
};

atom& molecule::get_atom_from_num(unsigned int num) {
	return atoms[num-1];	
};

void molecule::print_coords(std::vector<double> coords) {
	std::cout << "(";
	for(double i : coords) {
		std::cout << i << " ";
	};
	std::cout << ")" << std::endl;
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
	atom second_atom(number_of_atoms, atom_name, num_bond_atom);
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
	atom third_atom(number_of_atoms, atom_name, num_bond_atom, num_angle_atom);
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
		atom remaining_atom(number_of_atoms, atom_name, num_bond_atom, num_angle_atom, num_dihedral_atom);
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

void molecule::print_cart_coords() {
	for(int i = 1; i <= num_atoms; i++) {
		std::vector<double> coords = get_atom_coord(i);
		std::cout << "(";
		for(double j : coords) {
			std::cout << j << " ";
		};	
		std::cout << ")" << std::endl;
	};
	return;
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
			std::cerr << "Incorrect atom name" << std::endl;
			exit(1);
		};
		double x, y, z;
		iss >> x; iss >> y; iss >> z;
		this -> set_atom_coord(atom_counter, x, y, z);
	};
	return;
};

void molecule::set_molecule_coord_Z(std::string coord_file) {
	std::ifstream stream(coord_file);
	if(!stream) {
		std::cerr << "Error opening file " + coord_file << std::endl;
		exit(1);
	};
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
		if(!(iss >> name)){
			file_error_message();
			exit(1);
		};
		if((this -> get_atom_from_num(atom_counter)).get_name() != name) {
			std::cerr << "Incorrect atom name" << std::endl;
			exit(1);
		};
		int bond_atom;
		double r;
		int angle_atom;
		double angle;
		int di_angle_atom;
		double di_angle;
		if(atom_counter == 1) {
			this -> set_atom_coord(atom_counter, 0, 0, 0);
		} else if(atom_counter == 2) {

			if(!(iss >> bond_atom)){
				file_error_message();
				exit(1);
			};
			if(!(iss >> r)){
				file_error_message();
				exit(1);
			};

			this -> set_atom_coord(atom_counter, 0, 0, r);
		} else if(atom_counter == 3) {

			if(!(iss >> bond_atom)){
				file_error_message();
				exit(1);
			};
			if(!(iss >> r)){
				file_error_message();
				exit(1);
			};
			if(!(iss >> angle_atom)){
				file_error_message();
				exit(1);
			};
			if(!(iss >> angle)){
				file_error_message();
				exit(1);
			};

			int pm = 3  - 2*bond_atom;
			std::vector<double> bond_coords = get_atom_coord(bond_atom);
			set_atom_coord(atom_counter, 0, r*sin(angle), bond_coords[2] + r*pm*cos(angle));
		} else if(atom_counter > 3) {
			if(!(iss >> bond_atom)){
				file_error_message();
				exit(1);
			};
			if(!(iss >> r)){
				file_error_message();
				exit(1);
			};
			if(!(iss >> angle_atom)){
				file_error_message();
				exit(1);
			};
			if(!(iss >> angle)){
				file_error_message();
				exit(1);
			};
			if(!(iss >> di_angle_atom)){
				file_error_message();
				exit(1);
			};
			if(!(iss >> di_angle)){
				file_error_message();
				exit(1);
			};
			std::vector<double> angle_vec = get_atom_coord(angle_atom);
			std::vector<double> z_axis = vec_normalised( displacement_vector(get_atom_coord(bond_atom), angle_vec) );
			std::vector<double> y_axis = vec_normalised(cross_product(z_axis, displacement_vector(get_atom_coord(di_angle_atom),angle_vec)  ) );
			std::vector<double> x_axis = cross_product( y_axis, z_axis);
			std::vector<double> bond_coords = get_atom_coord(bond_atom);
			std::vector<double> added_vec = {0,0,0};
			for(int i = 0; i < added_vec.size(); i++) {
				added_vec[i] = r*(x_axis[i]*sin(angle)*cos(di_angle) 
						+ y_axis[i]*sin(angle)*sin(di_angle)
						- z_axis[i]*cos(angle));
			};
			set_atom_coord(atom_counter, bond_coords[0] + added_vec[0], bond_coords[1] + added_vec[1], 
				bond_coords[2] + added_vec[2]);	
		};
	};
	return;

};

void molecule::update_molecule_coord(const Eigen::VectorXd& vec) {
	if(vec.size() != 3*num_atoms) {
		std::cerr << "Incorrect vector size" << std::endl;
		exit(1);
	};
	for(int i = 1; i <= num_atoms; i++) {
		
		update_atom_coord(i, vec(3*(i-1)), vec(3*(i-1) + 1), vec(3*(i-1) + 2));
	};	
	return;  
};


void molecule::set_atom_coord(int atom_num, double x, double y, double z) {
	(this -> get_atom_from_num(atom_num)).set_cart_coord(x, y, z);
};	

void molecule::update_atom_coord(int atom_num, double dx, double dy, double dz) {
	(this -> get_atom_from_num(atom_num)).update_cart_coord(dx, dy, dz);
};

std::vector<std::vector<double>> molecule::atom_and_connected_coord(std::vector<int> connected_atoms) {
	std::vector<std::vector<double>> coordinates_of_atoms;
	for(int i : connected_atoms) {
		atom &current_atom = this -> get_atom_from_num(i);
		std::vector<double> vec = current_atom.get_cart_coord();
		coordinates_of_atoms.push_back(vec);
	};
	return coordinates_of_atoms;
};

std::vector<double> molecule::get_atom_coord(int atom_num) {
	atom &current_atom = this -> get_atom_from_num(atom_num);
	std::vector<double> vec = current_atom.get_cart_coord();
	return vec;
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
	atom& first_atom = this -> get_atom_from_num(atom_num);
	std::vector<int> connected_atoms = first_atom.get_connected_atoms();
	std::vector<std::vector<double>> coordinates_of_atoms = atom_and_connected_coord(connected_atoms);
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
	atom& first_atom = this -> get_atom_from_num(atom_num);
	std::vector<int> connected_atoms = first_atom.get_connected_atoms();
	std::vector<std::vector<double>> coordinates_of_atoms = atom_and_connected_coord(connected_atoms);
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
	atom& first_atom = this -> get_atom_from_num(atom_num);
	std::vector<int> connected_atoms = first_atom.get_connected_atoms();
	std::vector<std::vector<double>> coordinates_of_atoms = atom_and_connected_coord(connected_atoms);
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
	/*	
	std::vector<double> vec_1 = displacement_vector(coord_1, coord_2);
	std::vector<double> vec_2 = displacement_vector(coord_3, coord_2);
	std::vector<double> vec_3 = displacement_vector(coord_4, coord_3);

	std::vector<double> new_vec_1 = vec_normalised(cross_product(vec_1, vec_2));
	std::vector<double> new_vec_2 = vec_normalised(cross_product(vec_2, vec_3));
	std::vector<double> new_vec_3 = vec_normalised(vec_2);
	std::vector<double> new_vec_4 = cross_product(new_vec_1, new_vec_3);

	double x = dot_product(new_vec_1, new_vec_2), y = dot_product(new_vec_4, new_vec_2);
	
	return atan2(y,x);
	*/
	std::vector<double> x1 = displacement_vector(coord_1, coord_2);
	std::vector<double> z = displacement_vector(coord_2, coord_3);
	std::vector<double> x2 = displacement_vector(coord_4, coord_3);
	std::vector<double> first = cross_product(cross_product(z, x1), z);
	std::vector<double> second = cross_product(cross_product(z, x2),z);
	return calculate_vec_angle(first, second);	
};
 
double molecule::distance(std::vector<double> coord_1, std::vector<double> coord_2) {
	std::vector<double> displacement_vec = displacement_vector(coord_1, coord_2);
	return sqrt(dot_product(displacement_vec, displacement_vec));
};
/*
double molecule::bond_length_derivative(int atom_num, int second_atom_num, std::string axis) {
	atom& first_atom = this -> get_atom_from_num(atom_num);
	std::vector<int> connected_atoms = first_atom.get_connected_atoms();
	std::vector<std::vector<double>> coordinates_of_atoms = atom_and_connected_coord(connected_atoms);
	double q = coordinate_value(second_atom_num, axis);
	double dq = derivative_increment(q);
	int pos = connected_atom_pos(connected_atoms, second_atom_num);
	coordinates_of_atoms[pos][axes_name_to_num(axis)] += dq;
	double R_plus = distance(coordinates_of_atoms[0], coordinates_of_atoms[1]);
	coordinates_of_atoms[pos][axes_name_to_num(axis)] += -2*dq;
	double R_minus = distance(coordinates_of_atoms[0], coordinates_of_atoms[1]);
	return derivative_value(R_plus, R_minus, dq);	
};

double molecule::angle_derivative(int atom_num, int second_atom_num, std::string axis) {
	atom& first_atom = this -> get_atom_from_num(atom_num);
	std::vector<int> connected_atoms = first_atom.get_connected_atoms();
	
	std::vector<std::vector<double>> coordinates_of_atoms = atom_and_connected_coord(connected_atoms);
	
	double q = coordinate_value(second_atom_num, axis);
	double dq = derivative_increment(q);
	int pos = connected_atom_pos(connected_atoms, second_atom_num);	
	coordinates_of_atoms[pos][axes_name_to_num(axis)] += dq;
	double R_plus = calculate_angle(coordinates_of_atoms[0], coordinates_of_atoms[1], coordinates_of_atoms[2]);
	coordinates_of_atoms[pos][axes_name_to_num(axis)] += -2*dq;
	double R_minus = calculate_angle(coordinates_of_atoms[0], coordinates_of_atoms[1], coordinates_of_atoms[2]);
	return derivative_value(R_plus, R_minus, dq);
};

double molecule::dihedral_angle_derivative(int atom_num, int second_atom_num, std::string axis){
	atom& first_atom = this -> get_atom_from_num(atom_num);
	std::vector<int> connected_atoms = first_atom.get_connected_atoms();
	
	std::vector<std::vector<double>> coordinates_of_atoms = atom_and_connected_coord(connected_atoms);
	double q = coordinate_value(second_atom_num, axis);

	double dq = derivative_increment(q);
	int pos = connected_atom_pos(connected_atoms, second_atom_num);	
	coordinates_of_atoms[pos][axes_name_to_num(axis)] += dq;
	double R_plus = calculate_dihedral_angle(coordinates_of_atoms[0], coordinates_of_atoms[1], 
		coordinates_of_atoms[2], coordinates_of_atoms[3]);
	coordinates_of_atoms[pos][axes_name_to_num(axis)] += -2*dq;
	double R_minus = calculate_dihedral_angle(coordinates_of_atoms[0], coordinates_of_atoms[1],
		coordinates_of_atoms[2], coordinates_of_atoms[3]);
	return derivative_value(R_plus, R_minus, dq);
};

double molecule::derivative_increment(double q) {
	double dq = 0.0001;
	if(q < 0) {
		q = -q; 
	};
	if(q > pow(10,-8)) {
		dq = 0.001*q;
	} else {
		dq = 0.0001;
	};
	return dq;
};

double molecule::derivative_value(double R_plus, double R_minus, double dq) {
	return (R_plus - R_minus)/(2*dq);		
};
*/

double molecule::coordinate_value(int second_atom_num, std::string axis) {
	double q = (this -> get_atom_from_num(second_atom_num)).get_cart_coord()[axes_name_to_num(axis)];
	return q;
};	

int molecule::connected_atom_pos(std::vector<int> connected_atoms, int second_atom) {
	bool is_in = false;
	int pos = -1;
	for(int i : connected_atoms) {
		pos++;
		if( i == second_atom) {
			is_in = true;
			break;	
		}; 
	};
	if(!is_in) {
		std::cerr << "Second atom not connected to first" << std::endl;
		exit(1);
	};
	return pos;	
};

int molecule::axes_name_to_num(std::string axes) {
		if(axes == "x") {
			return 0;
		} else if (axes == "y") {
			return 1;
		} else if (axes == "z") {
			return 2;
		} else {
			std::cerr << "Axis not of valid type";
			exit(1);
		};
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

