#include "molecule_2.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <functional>
#include <map>



typedef int (atom::*coord_func)() const;

Eigen::MatrixXd molecule::Amat() {
	Eigen::MatrixXd A = Eigen::MatrixXd(3,3);
	for(int i = 1; i <= num_atoms; i++) {
		atom& cur_atom = get_atom_from_num(i);
		double mass = cur_atom.get_mass();
		Eigen::Vector3d coord = cur_atom.get_cart_coord(1);
		Eigen::Vector3d eq_coord = cur_atom.get_cart_coord(0);
		for(int j = 0; j < coord.size(); j++) {
			for(int k = 0; k < coord.size(); k++) {
				A(j,k) += mass*coord(j)*eq_coord(k);
			};
		};
	};
	return A;
};

Eigen::MatrixXd molecule::ATAmat(const Eigen::MatrixXd& A) {
	return A.transpose()*A;
};

Eigen::MatrixXd molecule::Tmat(const Eigen::MatrixXd& A, const Eigen::MatrixXd& ATA) {
	Eigen::EigenSolver<Eigen::MatrixXd> es;
	es.compute(ATA, true);

	Eigen::Vector3d vec_1 = es.eigenvectors().col(0).real();  
	
	vec_1 = vec_1/sqrt(vec_1.dot(vec_1));
	//vec_normalised(vec_1);	
	Eigen::Vector3d vec_2 = es.eigenvectors().col(1).real();
	//vec_normalised(vec_2);
	vec_2 = vec_2/sqrt(vec_2.dot(vec_2));
	Eigen::Vector3d vec_3 = vec_1.cross(vec_2);
	//vec_normalised(vec_3); 
	vec_3 = vec_3/sqrt(vec_3.dot(vec_3));
	std::vector<Eigen::Vector3d> first_set = {vec_1, vec_2, vec_3};

	Eigen::Vector3d wec_1 = A*vec_1;
	//vec_normalised(wec_1);
	wec_1 = wec_1/sqrt(wec_1.dot(wec_1));
	Eigen::Vector3d wec_2 = A*vec_2;
	//vec_normalised(wec_2);
	wec_2 = wec_2/sqrt(wec_2.dot(wec_2));
	Eigen::Vector3d wec_3 = wec_1.cross(wec_2);
	//vec_normalised(wec_3);
	wec_3 = wec_3/sqrt(wec_3.dot(wec_3));

	std::vector<Eigen::Vector3d> second_set = {wec_1, wec_2, wec_3}; 
	
	Eigen::MatrixXd T = Eigen::MatrixXd(3,3);
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			for(int k = 0; k < 3; k++) {
				T(i, j) += first_set[k][i]*second_set[k][j];
			};	
		};
	};
	return T;
};

Eigen::Vector3d vec_normalised(Eigen::Vector3d vec) {
	return vec/sqrt(vec.dot(vec));	
};

void molecule::rotate_coords(const Eigen::MatrixXd& T) {	
	for(int i = 1; i <= num_atoms; i++) {
		Eigen::Vector3d coords = get_atom_coord(1, i);
		Eigen::Vector3d values = {0,0,0};
		for(int j = 0; j < 3; j++) {
			for(int k = 0; k < 3; k++) {
				values[j] += coords[k]*T(j,k);	
			};
		};
		set_atom_coord(1, i, values[0], values[1], values[2]);
	};
	return;
};

Eigen::Vector3d molecule::Eckart_cond() {
	Eigen::Vector3d total = {0,0,0};
	for(int i = 1; i <= num_atoms; i++) {
		atom& cur_atom = get_atom_from_num(i);
		double cur_mass = cur_atom.get_mass();
		Eigen::Vector3d eq_coord = cur_atom.get_cart_coord(0);
		Eigen::Vector3d coord = cur_atom.get_cart_coord(1);
		Eigen::Vector3d cross = eq_coord.cross(coord);
		total += cur_mass*cross;
	};
	return total;
};	


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
	double mass;
	int number_of_atoms = 0;
	if(!(iss_1 >> atom_name)){
		file_error_message();
		exit(1);
	};
	if(!(iss_1 >> mass)) {
		file_error_message();
		exit(1);
	};
	number_of_atoms++;
	atoms.push_back(atom(number_of_atoms, mass, atom_name));

	if(!getline(stream, line)) {
		file_error_message();
		exit(1);
	};
	std::istringstream iss_2(line);
	if(!(iss_2 >> atom_name)){
		file_error_message();
		exit(1);
	};
	if(!(iss_2 >> mass)) {
		file_error_message();
		exit(1);
	}
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
	atom second_atom(number_of_atoms, mass,  atom_name, num_bond_atom);
	atoms.push_back(second_atom);

	if(!getline(stream, line)) {
		return;
	};
	std::istringstream iss_3(line);
	if(!(iss_3 >> atom_name)){
		file_error_message();
		exit(1);
	};
	if(!(iss_3 >> mass)) {
		file_error_message();
		exit(1);
	}
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
	atom third_atom(number_of_atoms, mass, atom_name, num_bond_atom, num_angle_atom);
	atoms.push_back(third_atom);
	while(getline(stream, line)) {
		std::istringstream iss_4(line);
		if(!(iss_4 >> atom_name)){
			file_error_message();
			exit(1);
		};
		if(!(iss_4 >> mass)) {
			file_error_message();
			exit(1);
		}
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
		atom remaining_atom(number_of_atoms, mass, atom_name, num_bond_atom, num_angle_atom, num_dihedral_atom);
		atoms.push_back(remaining_atom);		
	};
	num_atoms = number_of_atoms;
};

void molecule::print_coordinates(int type) {
	for(int i = 1; i <= num_atoms; i++) {
		if(i == 1) {
			std::cout << (*this).get_atom_from_num(i).get_name() << std::endl;
		} else if(i == 2) {
			atom& cur_atom = (*this).get_atom_from_num(i);
			std::cout << cur_atom.get_name() << " ";
			std::cout << cur_atom.get_bond_length_atom() << " ";
			std::cout << (this -> bond_length(type, i)) << std::endl;
		} else if(i == 3) {
			atom& cur_atom = (*this).get_atom_from_num(i);
			std::cout << cur_atom.get_name() << " ";
			std::cout << cur_atom.get_bond_length_atom() << " ";
			std::cout << (this -> bond_length(type, i)) << " ";
			std::cout << cur_atom.get_angle_atom() << " ";
			std::cout << (this -> angle(type, i)) << std::endl;;
		} else {
			atom& cur_atom = (*this).get_atom_from_num(i);
			std::cout << cur_atom.get_name() << " ";
			std::cout << cur_atom.get_bond_length_atom() << " ";
			std::cout << (this -> bond_length(type, i)) << " ";
			std::cout << cur_atom.get_angle_atom() << " ";
			std::cout << (this -> angle(type, i)) << " ";
			std::cout << cur_atom.get_dihedral_angle_atom() << " ";
			std::cout << (this -> dihedral_angle(type, i)) << std::endl;
		};
	};
};

void molecule::print_cart_coords(int type) {
	for(int i = 1; i <= num_atoms; i++) {
		Eigen::Vector3d coords = get_atom_coord(type, i);
		std::cout << coords.transpose()	<< std::endl;
	};
	return;
};

void molecule::set_molecule_coord(int type, std::string coord_file) {
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
		this -> set_atom_coord(type, atom_counter, x, y, z);
	};
	move_to_COM(type);
	return;
};

void molecule::set_molecule_coord_Z(int type, std::string coord_file) {
	std::ifstream stream(coord_file);
	if(!stream) {
		std::cerr << "Error opening file " + coord_file << std::endl;
		exit(1);
	};
	std::string line;
	int atom_counter = 0;
	int bond_atom;
	double r;
	int angle_atom;
	double angle;
	int di_angle_atom;
	double di_angle;
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
		atom& cur_atom = get_atom_from_num(atom_counter);

		if(atom_counter == 1) {
			cur_atom.set_cart_coord(type, 0, 0, 0);
		} else if(atom_counter == 2) {
			if(!(iss >> bond_atom)){
				file_error_message();
				exit(1);
			};
			if(!(iss >> r)){
				file_error_message();
				exit(1);
			};
	
			if(cur_atom.get_bond_length_atom() != bond_atom) {
				Z_coord_error();
				exit(1);
			};
		
			cur_atom.set_cart_coord(type, 0, 0, r);
		
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
			
			if( (cur_atom.get_bond_length_atom() != bond_atom) && 
			(cur_atom.get_angle_atom() != angle_atom) ) {
				Z_coord_error();
				exit(1); 
			};

			int pm = 3  - 2*bond_atom;
			Eigen::Vector3d bond_coords = get_atom_coord(type, bond_atom);
			cur_atom.set_cart_coord(type, 0, r*sin(angle), bond_coords[2] + r*pm*cos(angle));
		
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
	
			if( (cur_atom.get_bond_length_atom() != bond_atom) && 
			(cur_atom.get_angle_atom() != angle_atom) && 
		    (cur_atom.get_dihedral_angle_atom() != di_angle_atom) ){
				Z_coord_error();
				exit(1); 
			};
		
			Eigen::Vector3d angle_vec = get_atom_coord(type, angle_atom);
			Eigen::Vector3d z_axis =  get_atom_coord(type, bond_atom) - angle_vec;
			z_axis = z_axis/sqrt(z_axis.dot(z_axis));
			Eigen::Vector3d y_axis = z_axis.cross(get_atom_coord(type, di_angle_atom) - angle_vec);
			y_axis = y_axis/sqrt(y_axis.dot(y_axis));
			Eigen::Vector3d x_axis = y_axis.cross(z_axis);
			Eigen::Vector3d bond_coords = get_atom_coord(type, bond_atom);
			Eigen::Vector3d added_vec = {0,0,0};
				added_vec = r*(x_axis*sin(angle)*cos(di_angle) 
						+ y_axis*sin(angle)*sin(di_angle)
						- z_axis*cos(angle));
			cur_atom.set_cart_coord(type, bond_coords[0] + added_vec[0], bond_coords[1] + added_vec[1], 
				bond_coords[2] + added_vec[2]);	
		};
	};
	move_to_COM(type);
	return;

};

void molecule::move_to_COM(int type) {
	Eigen::Vector3d com = {0,0,0};
	double total_mass; 
	for(int i = 1; i <= num_atoms; i++) {
		atom& cur_atom = get_atom_from_num(i);
		Eigen::Vector3d coords = cur_atom.get_cart_coord(type);
		double cur_mass = cur_atom.get_mass();
		com += cur_mass*coords;
		total_mass += cur_mass;
	};		
	com  = com/total_mass;
	for(int i = 1; i <= num_atoms; i++) {
		update_atom_coord(type, i, -com[0], -com[1], -com[2]);
	};
};

void molecule::update_molecule_coord(const Eigen::VectorXd& vec) {

	if(vec.size() != 3*num_atoms) {
		std::cerr << "Incorrect vector size" << std::endl;
		exit(1);
	};
	for(int i = 1; i <= num_atoms; i++) {
		
		update_atom_coord(0, i, vec(3*(i-1)), vec(3*(i-1) + 1), vec(3*(i-1) + 2));
	};	
	return;  
};


void molecule::set_atom_coord(int type, int atom_num, double x, double y, double z) {
	(this -> get_atom_from_num(atom_num)).set_cart_coord(type, x, y, z);
};	

void molecule::update_atom_coord(int type, int atom_num, double dx, double dy, double dz) {
	(this -> get_atom_from_num(atom_num)).update_cart_coord(type,  dx, dy, dz);
};

std::vector<Eigen::Vector3d> molecule::atom_and_connected_coord(int type, std::vector<int> connected_atoms) {
	std::vector<Eigen::Vector3d> coordinates_of_atoms;
	for(int i : connected_atoms) {
		atom &current_atom = this -> get_atom_from_num(i);
		Eigen::Vector3d vec = current_atom.get_cart_coord(type);
		coordinates_of_atoms.push_back(vec);
	};
	return coordinates_of_atoms;
};


Eigen::Vector3d molecule::get_atom_coord(int type, int atom_num) {
	atom &current_atom = this -> get_atom_from_num(atom_num);
	Eigen::Vector3d vec = current_atom.get_cart_coord(type);
	return vec;
};

double molecule::bond_length(int type, int atom_num) {
	if(atom_num == 1) {
		std::cerr << "No bond length defined for atom one." << std::endl;
		exit(1);
	};
	/*atom& first_atom = this -> get_atom_from_num(atom_num);
	atom& second_atom = this -> get_atom_from_num(first_atom.get_bond_length_atom());
	std::vector<double> vec_1 = first_atom.get_cart_coord();
	std::vector<double> vec_2 = second_atom.get_cart_coord();
	*/
	atom& first_atom = this -> get_atom_from_num( atom_num);
	std::vector<int> connected_atoms = first_atom.get_connected_atoms();
	std::vector<Eigen::Vector3d> coordinates_of_atoms = atom_and_connected_coord(type, connected_atoms);
	return distance(coordinates_of_atoms[0], coordinates_of_atoms[1]); 
};

double molecule::angle(int type, int atom_num) {
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
	std::vector<Eigen::Vector3d> coordinates_of_atoms = atom_and_connected_coord(type, connected_atoms);
	return calculate_angle(coordinates_of_atoms[0], coordinates_of_atoms[1], coordinates_of_atoms[2]); 
};

double molecule::dihedral_angle(int type, int atom_num) {
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
	std::vector<Eigen::Vector3d> coordinates_of_atoms = atom_and_connected_coord(type, connected_atoms);
	return calculate_dihedral_angle(coordinates_of_atoms[0], coordinates_of_atoms[1], coordinates_of_atoms[2],
			coordinates_of_atoms[3]); 
};

//double molecule::dot_product(std::vector<double> coord_1, std::vector<double> coord_2) {
//	if(!check_coord_size(coord_1) || !check_coord_size(coord_2)) {
//		coord_length_error_message();
//		exit(1);	
//	};	
//	return coord_1[0]*coord_2[0] + 
//		coord_1[1]*coord_2[1] + 
//		coord_1[2]*coord_2[2];
//};

//std::vector<double> molecule::cross_product(std::vector<double> coord_1, std::vector<double> coord_2) {
//	if(!check_coord_size(coord_1) || !check_coord_size(coord_2)) {
//		coord_length_error_message();
//		exit(1);	
//	}
//	std::vector<double> result = {coord_1[1]*coord_2[2] - coord_1[2]*coord_2[1], 
//		coord_1[2]*coord_2[0] - coord_1[0]*coord_2[2], coord_1[0]*coord_2[1] - coord_2[0]*coord_1[1]};
//	return result;
//};
//std::vector<double> molecule::displacement_vector(std::vector<double> coord_1, std::vector<double> coord_2) {
//	if(!check_coord_size(coord_1) || !check_coord_size(coord_2)) {
//		coord_length_error_message();
//		exit(1);	
//	};
//	std::vector<double> displacement = {coord_1[0] -coord_2[0], coord_1[1] - coord_2[1], coord_1[2] -coord_2[2]};
//	return displacement;
//};

double molecule::calculate_angle(Eigen::Vector3d coord_1, Eigen::Vector3d coord_2, Eigen::Vector3d coord_3) {
	Eigen::Vector3d vec_1 = coord_1 - coord_2;
	Eigen::Vector3d vec_2 = coord_3 - coord_2;
	return calculate_vec_angle(vec_1, vec_2);	
};

double molecule::calculate_vec_angle(Eigen::Vector3d vec_1, Eigen::Vector3d vec_2) {
	double numerator = vec_1.dot(vec_2);
	double denominator = vec_distance(vec_1)*vec_distance(vec_2);
	double angle = acos(numerator/denominator);
	return angle;
};

double molecule::vec_distance(Eigen::Vector3d vec) {
	return sqrt(vec.dot(vec));
};

//std::vector<double> molecule::vec_normalised(std::vector<double> vec) {
//	std::vector<double> n_vec;
//	double length = vec_distance(vec);
//	for(int i = 0; i < vec.size(); i++) {
//		n_vec.push_back(vec[i]/length);
//	};	
//	return n_vec;
//};

double molecule::calculate_dihedral_angle(Eigen::Vector3d coord_1, Eigen::Vector3d coord_2,
	Eigen::Vector3d coord_3, Eigen::Vector3d coord_4) {
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
	Eigen::Vector3d x1 = coord_1 - coord_2;
	Eigen::Vector3d z = coord_2 - coord_3;
	Eigen::Vector3d x2 = coord_4 - coord_3;
	Eigen::Vector3d first = ( z.cross(x1) ).cross(z);
	Eigen::Vector3d second = ( z.cross(x2) ).cross(z);
	return calculate_vec_angle(first, second);	
};
 
double molecule::distance(Eigen::Vector3d coord_1, Eigen::Vector3d coord_2) {
	Eigen::Vector3d displacement_vec = coord_1 - coord_2;
	return vec_distance(displacement_vec);
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

double molecule::coordinate_value(int type, int second_atom_num, std::string axis) {
	double q = (this -> get_atom_from_num(second_atom_num)).get_cart_coord(type)[axes_name_to_num(axis)];
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


void molecule::Z_coord_error() {
	std::cerr << "The inputted Z matrix does not match what is stored in this molecule" << std::endl;
};
