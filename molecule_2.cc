#include "molecule_2.h"
#include "matrix_op.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <functional>
#include <map>



	typedef int (atom::*coord_func)() const;
	const double proton_mass = 1822.888;
	const double threshold = 0.0001;
	molecule::coefficient::coefficient(int num_modes, int num_order, int num_dimension) : modes(num_modes), 
	order(num_order), dimension(num_dimension), coeffs(pow(num_modes, num_dimension), std::vector<double>(pow(num_order,num_dimension))) {}; 

	molecule::coefficient::coefficient() {};

	void molecule::coefficient::set_coefficient_2D(int mode_1, int mode_2, int order_1, int order_2, int val) {
		coeffs[mode_1*modes+ mode_2][order_1*order+ order_2] = val;	
		};
	double molecule::coefficient::get_coefficient_2D(int mode_1, int mode_2, int order_1, int order_2){ 
		return coeffs[mode_1*modes + mode_2][order_1*order + order_2];
	};
	void molecule::coefficient::set_coefficient_1D(int mode, int order, double val) {
		coeffs[mode][order] = val;
	};
	double molecule::coefficient::get_coefficient_1D(int mode, int order) {
		return coeffs[mode][order]; 
	};

	void molecule::print_coefficients_1D() {
		coeffs_1D.print_coefficients();
	};

	void molecule::coefficient::print_coefficients() {
		for(std::vector<double> i  : coeffs) {
			for(double j : i) {	
				std::cout << j << " ";
			};	
			std::cout << std::endl;	
		};
	};

	int molecule::coefficient::get_modes() {
		return modes;
	};

	int molecule::coefficient::get_order() {
		return order;
	};
	
	int molecule::coefficient::get_col_num() {
		return col_num;
	};
	
	molecule::grid_coeffs::grid_coeffs(std::vector<int> num_modes, std::vector<std::string> input_labels, 
		Eigen::MatrixXd input_coeffs, int input_poly_order, int input_col) : 
			modes(num_modes), labels(input_labels), coeffs(input_coeffs),  poly_order(input_poly_order), column(input_col) {};

	std::vector<int> molecule::grid_coeffs::get_modes() {return modes;};
	
	int molecule::grid_coeffs::get_poly_order() {return poly_order;};
	
	Eigen::MatrixXd molecule::grid_coeffs::get_coeffs() {return coeffs;};
	
	std::vector<std::string> molecule::grid_coeffs::get_labels() {return labels;};
	
	int molecule::grid_coeffs::get_column() { return column;};

	Eigen::MatrixXd molecule::empty_matrix() {
		return Eigen::MatrixXd(3*num_atoms, 3*num_atoms); 
	};


	Eigen::MatrixXd molecule::Amat() {
		Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3,3);
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


	Eigen::MatrixXd molecule::Tmat(const Eigen::MatrixXd& A) {
		Eigen::MatrixXd ATA = A.transpose()*A;
		Eigen::MatrixXd AAT  = A*A.transpose();
		std::cout << ATA << std::endl;
		std::cout << AAT << std::endl;
		Eigen::EigenSolver<Eigen::MatrixXd> es;
		es.compute(ATA, true);
		Eigen::Vector3d vec_1 = es.eigenvectors().col(0).real();  
		vec_1 = vec_1/sqrt(vec_1.dot(vec_1));
		Eigen::Vector3d vec_2 = es.eigenvectors().col(1).real();
		vec_2 = vec_2/sqrt(vec_2.dot(vec_2));
		Eigen::Vector3d vec_3 = vec_1.cross(vec_2);
		vec_3 = vec_3/sqrt(vec_3.dot(vec_3));
		std::vector<Eigen::Vector3d> first_set = {vec_1, vec_2, vec_3};
		
		Eigen::EigenSolver<Eigen::MatrixXd> ef;
		ef.compute(AAT, true);
		int skip;
		Eigen::Vector3d wec_1;
		Eigen::Vector3d wec_2;
	
		std::cout << "ATA " << es.eigenvalues() << std::endl;
		std::cout << "AAT " << ef.eigenvalues() << std::endl;	
		for(int i = 0; i < 3; i++) {
			if(fabs(es.eigenvalues()(0) - ef.eigenvalues()(i)) < threshold) {
				wec_1 = ef.eigenvectors().col(i).real();
				wec_1 = wec_1/sqrt(wec_1.dot(wec_1));
				skip = i; 
				break;
			};
		};	
		for(int i = 0; i < 3; i++) {
			if( i == skip) {
				continue;
			};
			if(fabs(es.eigenvalues()(1) - ef.eigenvalues()(i)) < threshold) {
				wec_2 = ef.eigenvectors().col(i).real();
				wec_2 = wec_2/sqrt(wec_2.dot(wec_2));
				break;
			};
		};
		
		Eigen::Vector3d wec_3 = wec_1.cross(wec_2);
		wec_3 = wec_3/sqrt(wec_3.dot(wec_3));

		std::vector<Eigen::Vector3d> second_set = {wec_1, wec_2, wec_3}; 

		Eigen::MatrixXd T = Eigen::MatrixXd::Zero(3,3);
		for(int i = 0; i < 3; i++) {
			for(int j = 0; j < 3; j++) {
				for(int k = 0; k < 3; k++) {
					T(i, j) += first_set[k](i)*second_set[k](j);
				};	
			};
		};
		std::cout << T << std::endl;
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

	Eigen::VectorXd molecule::normal_coordinates() {
		Eigen::VectorXd normal(num_atoms*3);
		for(int i = 0; i < num_atoms; i++) {
			atom& cur_atom = get_atom_from_num(i+1);
			Eigen::Vector3d atom_displacement = cur_atom.get_cart_coord(1) - cur_atom.get_cart_coord(0);
			for(int j = 0; j < 3; j++) {
				normal(i*3 + j) = atom_displacement(j);	
			};
		};
		return L_mat.transpose()*sqrt(proton_mass)*M_mat*normal;
	};


	molecule::molecule(std::string z_matrix_file, std::string molecule_name) : 
		name(molecule_name) {
		std::vector<double> masses;
		std::ifstream stream(z_matrix_file);
		if(!stream) {
			std::cerr << "Error opening file " + z_matrix_file << std::endl;
			exit(1);
		};
		std::string line;

		if(!getline(stream, line)){
			file_error_message(z_matrix_file);
			exit(1);	
		};
		std::istringstream iss_1(line);
		std::string atom_name;
		double mass;
		int number_of_atoms = 0;
		if(!(iss_1 >> atom_name)){
			file_error_message(z_matrix_file);
			exit(1);
		};
		if(!(iss_1 >> mass)) {
			file_error_message(z_matrix_file);
			exit(1);
		};
		number_of_atoms++;
		atoms.push_back(atom(number_of_atoms, mass, atom_name));
		masses.insert(masses.end(), 3, sqrt(mass));

		if(!getline(stream, line)) {
			file_error_message(z_matrix_file);
			exit(1);
		};
		std::istringstream iss_2(line);
		if(!(iss_2 >> atom_name)){
			file_error_message(z_matrix_file);
			exit(1);
		};
		if(!(iss_2 >> mass)) {
			file_error_message(z_matrix_file);
			exit(1);
		}
		int num_bond_atom = 0;
		if(!(iss_2 >> num_bond_atom)){
			file_error_message(z_matrix_file);
			exit(1);
		};
		if(num_bond_atom > number_of_atoms) {
			file_error_message(z_matrix_file);
			exit(1);
		};
		number_of_atoms++;
		atom second_atom(number_of_atoms, mass,  atom_name, num_bond_atom);
		atoms.push_back(second_atom);
		
		masses.insert(masses.end(), 3, sqrt(mass));

		if(!getline(stream, line)) {
			return;
		};
		std::istringstream iss_3(line);
		if(!(iss_3 >> atom_name)){
			file_error_message(z_matrix_file);
			exit(1);
		};
		if(!(iss_3 >> mass)) {
			file_error_message(z_matrix_file);
			exit(1);
		};
		if(!(iss_3 >> num_bond_atom)){
			file_error_message(z_matrix_file);
			exit(1);
		};
		if(num_bond_atom > number_of_atoms) {
			file_error_message(z_matrix_file);
			exit(1);
		};
		int num_angle_atom = 0;
		if(!(iss_3 >> num_angle_atom)){
			file_error_message(z_matrix_file);
			exit(1);
		};
		if(num_angle_atom > number_of_atoms) {
			file_error_message(z_matrix_file);
			exit(1);
		};
		number_of_atoms++;
		atom third_atom(number_of_atoms, mass, atom_name, num_bond_atom, num_angle_atom);
		atoms.push_back(third_atom);
		masses.insert(masses.end(), 3, sqrt(mass));

		while(getline(stream, line)) {
			std::istringstream iss_4(line);
			if(!(iss_4 >> atom_name)){
				file_error_message(z_matrix_file);
				exit(1);
			};
			if(!(iss_4 >> mass)) {
				file_error_message(z_matrix_file);
				exit(1);
			}
			if(!(iss_4 >> num_bond_atom)){
				file_error_message(z_matrix_file);
				exit(1);
			};
			if(num_bond_atom > number_of_atoms) {
				file_error_message(z_matrix_file);
				exit(1);
			};
			if(!(iss_4 >> num_angle_atom)){
				file_error_message(z_matrix_file);
				exit(1);
			};
			if(num_angle_atom > number_of_atoms) {
				file_error_message(z_matrix_file);
				exit(1);
			};	
			int num_dihedral_atom;
			if(!(iss_4 >> num_dihedral_atom)){
				file_error_message(z_matrix_file);
				exit(1);
			};
			number_of_atoms++;
			atom remaining_atom(number_of_atoms, mass, atom_name, num_bond_atom, num_angle_atom, num_dihedral_atom);
			atoms.push_back(remaining_atom);
			masses.insert(masses.end(), 3, sqrt(mass));	
		};
		num_atoms = number_of_atoms;
		M_mat.resize(num_atoms*3, num_atoms*3); 
		M_mat = Eigen::MatrixXd::Zero(num_atoms*3, num_atoms*3);
		for(int i = 0; i < num_atoms*3; i++) {
			M_mat(i,i) = masses[i];  
		};
		return;
	};

	void molecule::set_L_matrix(std::string L_matrix_file) {
		std::ifstream stream(L_matrix_file);
		if(!stream) {
			std::cerr << "Error opening file " + L_matrix_file << std::endl;
			exit(1);
		};
		std::string line;
		std::string word;
		while(getline(stream, line)) {
			std::istringstream iss_1(line);
			iss_1 >> word;
			if(word != "###") {
				continue;
			} else {
				iss_1 >> word;
				if(word != "DISPLACEMENT") {
					continue;
				};
				break;
			};
		};
		while(getline(stream , line)) {
			if(line.length() == 0) {
				continue;
			} else {
				break; 
			};
		};
		double frequency;
		
		std::istringstream iss_2(line);
		iss_2 >> word;
		iss_2 >> word;

		while(iss_2 >> frequency) {
			frequencies.push_back(frequency);					
		};
		getline(stream , line);
		getline(stream , line);
		double component;
		std::istringstream iss_3(line);
		iss_3 >> word;
		L_mat = Eigen::MatrixXd::Zero(num_atoms*3, num_atoms*3 - 6);
		
		for(int i = 0; i < num_atoms*3; i++) {
			if(i!=0) {
				getline(stream, line);
				std::istringstream iss_4(line);
				int j = 0;
				while(iss_4 >> component) {
					L_mat(i,j) = component*sqrt(proton_mass);
					j++;
				};
			} else {
				std::istringstream iss_4(line);
				int j = 0;
				while(iss_4 >> component) {
					L_mat(i,j) = component*sqrt(proton_mass);
					j++;
				};
			};
		}; 
		L_mat = M_mat*L_mat;	
		std::cout << "L mat end " << L_mat << std::endl;	
		L_matrix_set = true;
		return; 	
	};

	void molecule::set_coefficients(std::string coeff_file) {
		std::ifstream stream(coeff_file);
		if(!stream) {
			std::cerr << "Error opening file" + coeff_file << std::endl;
			exit(1);
		};
		std::string line;
		std::string word_1;
		std::string word_2;
		while(getline(stream, line)) {
			std::istringstream iss_1(line);
			iss_1 >> word_1;
			iss_1 >> word_2; 
			if(word_1 != "1D" or word_2 != "polynomials") {
				continue;
			} else {
				break;
			};
		};
		getline(stream, line);
		getline(stream, line);
		std::istringstream iss_2(line);
		iss_2>> word_1; 
		double val;
		int order;
		while(iss_2 >> val) {
			order = val;
		};
		int modes = num_atoms*3 - 6;
		coeffs_1D = coefficient(modes, order, 1);
		getline(stream, line);
		
		int mode;
		for(int i = 0; i < modes; i++) {
			getline(stream, line);
			std::istringstream iss_3(line); 
			iss_3 >> mode;
			for(int j = 0; j < order; j++) {
				iss_3 >> val;
				coeffs_1D.set_coefficient_1D(mode - 1, j, val);
			};
		};

		while(getline(stream, line)) {
			std::istringstream iss_1(line);
			iss_1 >> word_1;
			iss_1 >> word_2; 
			if(word_1 != "2D" or word_2 != "polynomials") {
				continue;
			} else {
				break;
			};
		}
		getline(stream, line);
		getline(stream, line);
		std::istringstream iss_4(line);
		iss_4 >> word_1; 
		while(iss_4 >> val) {
			order = val;
		};
		coeffs_2D = coefficient(modes, order, 2);
		getline(stream, line);
		
		int mode_1; int mode_2; int order_1; int order_2;
		
		for(int k = 0; k < 3; k++) {
			for(int i = 0; i < order; i++) {
				getline(stream, line);
				std::istringstream iss_5(line); 
				iss_5 >> mode_1; iss_5 >> mode_2;
				iss_5 >> val;
				for(int j = 0; j < order; j++) {
					iss_5 >> val;
					coeffs_2D.set_coefficient_2D(mode_1 - 1, mode_2 - 1, j, i, val);
				};
			};
			getline(stream, line);
		};
		coefficients_set = true;
		return; 
	};
	
void molecule::set_grid_coeffs(std::string grid_file, int input_poly_order) {
	std::ifstream stream(grid_file);
	std::string line;
	if(!stream) {
		std::cerr << "Error opening file" + grid_file << std::endl;
		exit(1);
	};
	int order;
	std::string rest_of_sentence = find_line(stream, 5, "Order of the potential energy ");	
	std::istringstream iss_1(rest_of_sentence);
	std::string word;
	iss_1 >> word;
	iss_1 >> word;
	iss_1 >> word;
	iss_1 >> expansion_order;
	
	multi_level_mat = Eigen::MatrixXi::Zero(expansion_order, expansion_order);	
	grid_coeffs_vector.resize(expansion_order);
	
	rest_of_sentence = find_line(stream, 3, "### MULTI LEVEL ");
	for(int i = 0; i < expansion_order; i++) {
		getline(stream, line);
		std::istringstream iss_2(line);
		int mat_element;
		iss_2 >> word;
		dimension_labels.push_back(word);
		for(int j = 0; j < expansion_order; j++) {
			iss_2 >> mat_element; 
			multi_level_mat(i, j) = mat_element;
		};
	};
	
	for(int cur_order = 1; cur_order <= expansion_order; cur_order++) {
		std::stringstream ss;
		ss <<"### " << cur_order << "D SURFACES ";
		std::string sentence = ss.str();
		rest_of_sentence = find_line(stream, 3, sentence);
		
		getline(stream, line);
		std::istringstream iss_3(line);
		int num_columns;
		int num_surfaces;
		for(int j = 0; j < 5; j++) {
			iss_3 >> line;
		};	
		iss_3 >> num_columns;
		num_columns += cur_order; 
		getline(stream, line);
		getline(stream, line);
		
		std::istringstream iss_4(line);
		for(int j = 0; j < 3; j++) {
			iss_4 >> line;
		};
		iss_4 >> num_surfaces;
		
		for(int j = 0; j < num_surfaces; j++) {
			std::stringstream nss;
			nss << "# " << cur_order << "D SURFACE: MODES: ";  
			sentence = nss.str();
			rest_of_sentence = find_line(stream, 4, sentence);
			std::istringstream iss_5(rest_of_sentence);
			int mode;
			std::vector<int> modes;
			
			for(int k = 0; k < cur_order; k++) {
				iss_5 >> mode;
				modes.push_back(mode);
			};
			sort_vec(modes);
			getline(stream, line);
			std::istringstream iss_6(line);
			for(int k = 0; k < 8; k++) {
				iss_6 >> word;
			};
			int num_ab_points;
			if(!(num_ab_points = std::stoi(word))) {
				std::string num_points_str = word.substr(0, word.length()-1);
				num_ab_points = std::stoi(word);
			};
			rest_of_sentence = find_line(stream, 1, "Grid "); 	
			
			std::istringstream iss_7(rest_of_sentence);
			iss_7 >> word;	
			
			bool has_dipole = false;
			if((iss_7 >> word) && word == "Dipole") {
				has_dipole = true;
			};	
			getline(stream, line);
			std::istringstream iss_8(line);
			int k = 0;
			for(int l = 0; l < cur_order; l++) {
				k++;
				iss_8 >> word;
			};
			iss_8 >> word;
			std::vector<std::string> energy_labels;
			energy_labels.push_back(word);
			k++;
			while(k < num_columns) {
				if(has_dipole) {
					for(int l = 0; l < 3; l++) {
						iss_8 >> word;
						k++;
					};
				};
				if(!(k < num_columns)) {
					break;
				};
				iss_8 >> word;
				k++;
				energy_labels.push_back(word);		
			};
			
			int correct_col;
			if(has_dipole) {
				correct_col = (num_columns - cur_order)/4 + cur_order;
			} else {
				correct_col = num_columns;
			};
			Eigen::MatrixXd grid_points = Eigen::MatrixXd::Zero(num_ab_points, correct_col); 
			for(int l = 0; l < num_ab_points; l++) {
				double value;
				getline(stream, line);
				std::istringstream iss_9(line);
				k = 0;
				for(int m = 0; m < cur_order; m++) {
					iss_9 >> value;
					grid_points(l, k) = value;
					k++;
				};
				iss_9 >> value;
				grid_points(l, k) = value; 
				k++;
				while(k < correct_col) {
					if(has_dipole) {
						for(int n = 0; n < 3; n++) {
							iss_9 >> value;
						};
					};
					if(!(k < correct_col)) {
						break;
					};
					iss_9 >> value;
					grid_points(l, k) = value;
					k++;
				};
			};
			std::vector<Eigen::MatrixXd> inverted_design_mat = 
				molecule::inverted_design_matrices(grid_points, cur_order, input_poly_order);
		
			for(int col = 1; col < correct_col - cur_order; col++) {
				Eigen::MatrixXd pot = get_V(grid_points, expansion_order, col, input_poly_order, modes); 
				std::cout << pot << std::endl;
				Eigen::MatrixXd fit = fitting_coefficients(0, pot, inverted_design_mat);
				molecule::grid_coeffs coeffs_object(modes, energy_labels, fit, input_poly_order, col);
				grid_coeffs_vector[cur_order - 1].push_back(coeffs_object);
			};
			for(int m = 1; m <= cur_order; m++) {
				int num = correct_energy_col(cur_order, m);
			std::cout << "order = " << cur_order << " sub order = " <<  m << " col = "  << num << std::endl;
			};
		};
	};
};

Eigen::MatrixXd  molecule::fitting_coefficients(int level, const Eigen::MatrixXd& V, const std::vector<Eigen::MatrixXd>& inverted_design_mat) {
	
	int dim = inverted_design_mat.size();
	int num_grid = inverted_design_mat[0].rows(); int num_basis = inverted_design_mat[0].cols();
	std::cout << "num_grid " << num_grid << " dim " << dim << " V.rows() " <<V.rows() << std::endl;
	if(V.rows() != pow(num_grid, dim)) {
		std::cerr << "Not enough potenial energies" << std::endl;
		exit(1);
	};

	Eigen::MatrixXd reshaped_V = V;	
	reshaped_V.resize(pow(num_grid, dim-level-2), pow(num_grid, 2));
	Eigen::MatrixXd V1 = Eigen::MatrixXd::Zero(pow(num_basis, dim-level-2),pow(num_grid, 2));
	if((dim-level)==3) {
		V1 = inverted_design_mat[dim-1]*reshaped_V;
	} else if((dim-level)==4) {
		V1 = outer_product(inverted_design_mat[dim-2], inverted_design_mat[dim-1])*reshaped_V;
	} else {
		for(int i = 0; i < pow(num_grid, 2); i++) {
			Eigen::MatrixXd V1_tmp = fitting_coefficients(level + 2, reshaped_V.col(i), inverted_design_mat);
			V1.col(i) = V1_tmp;
		};
	};
	Eigen::MatrixXd coeffs = V1*outer_product( (inverted_design_mat[level]).transpose(), (inverted_design_mat[level+1]).transpose());
};

Eigen::MatrixXd molecule::get_V(const Eigen::MatrixXd& grid_points, int order, int column, int poly_order, const std::vector<int>& modes) {
	std::cout << "column " << column << " order " << order << std::endl;
	Eigen::MatrixXd V = grid_points.col(column - 1 + order - 1);
	std::cout << V << std::endl;
	for(int i = 0; i < order - 1; i++) {
		int correct_col = correct_energy_col(order, i + 1)  + column - 1;
		std::cout << "correct_col, i, order, column " << correct_col << " " << i << " " << order << " " << column << std::endl;
		for(grid_coeffs coeff : grid_coeffs_vector[i]) {
			if( (!contains_all_nums(modes, coeff.get_modes())) or (correct_col != coeff.get_column() ) ) {
				continue;
			};	
			for(int j = 0; j < V.rows() ; j++) {
				Eigen::MatrixXd coordinates = Eigen::MatrixXd::Zero(3*num_atoms-6 , 1);
				for(int mode : coeff.get_modes()) {
					for(int k = 0; k < modes.size(); k++) {	
						if(modes[j] == mode) {
							coordinates(mode-1) = grid_points(i, j);
						};
					};
				};
				std::cout << coordinates << std::endl;
				Eigen::MatrixXd poly_mat = poly_factor_mat(coordinates, poly_order);
				Eigen::MatrixXd factor_mat = poly_mat.row(0);
				for(int k = 0; k < poly_mat.rows() - 1; k++) {
					Eigen::MatrixXd temp_mat = outer_product(factor_mat, poly_mat.row(k + 1));
					factor_mat = temp_mat;
				};
				for(int k = 0; k < factor_mat.size(); k++) {
					factor_mat(k) *= coeff.get_coeffs()(k);
					V(j) -= factor_mat(k);
				};
			};	 
		}; 
			
	};
	
};

int molecule::correct_energy_col(int cur_order, int sub_order) {
	std::vector<int> dim_tally;
	dim_tally.push_back(1);
	for(int i = sub_order - 1; i < cur_order - 1; i++) {
		int mat_val = multi_level_mat(sub_order - 1, i);
		std::cout << "mat_val  "  << mat_val << " (i,j) = " << sub_order -1  << " "  << i << std::endl;
		std::cout << " back " << dim_tally.back() << std::endl;
		if( mat_val >= 0 ) {
			dim_tally.push_back(dim_tally.back() + mat_val);
		} else {
			dim_tally.push_back(dim_tally[-mat_val-sub_order] + 1);
		};
	};
	return dim_tally.back();
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
		atom& cur_atom = get_atom_from_num(i); 
		Eigen::Vector3d coords = cur_atom.get_cart_coord(type); 
		std::cout << cur_atom.get_name() << coords.transpose()	<< std::endl;
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
		if(!(iss >> x)) {
			file_error_message(coord_file);
			exit(0);		
		} 
		if(!(iss >> y)) {
			file_error_message(coord_file);
			exit(0);
		};	
		if(!(iss >> z)) {
			file_error_message(coord_file);
			exit(0);
		};	
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
			file_error_message(coord_file);
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
				file_error_message(coord_file);
				exit(1);
			};
			if(!(iss >> r)){
				file_error_message(coord_file);
				exit(1);
			};
			if(cur_atom.get_bond_length_atom() != bond_atom) {
				Z_coord_error();
				exit(1);
			};
		
			cur_atom.set_cart_coord(type, 0, 0, r);
		
		} else if(atom_counter == 3) {
			
			if(!(iss >> bond_atom)){
				file_error_message(coord_file);
				exit(1);
			};
			if(!(iss >> r)){
				file_error_message(coord_file);
				exit(1);
			};
			if(!(iss >> angle_atom)){
				file_error_message(coord_file);
				exit(1);
			};
			if(!(iss >> angle)){
				file_error_message(coord_file);
				exit(1);
			};
				
			if( (cur_atom.get_bond_length_atom() != bond_atom) or 
			(cur_atom.get_angle_atom() != angle_atom) ) {
				Z_coord_error();
				exit(1); 
			};

			int pm = 3 - 2*bond_atom;
			Eigen::Vector3d bond_coords = get_atom_coord(type, bond_atom);
			cur_atom.set_cart_coord(type, 0, r*sin(angle), bond_coords[2] + r*pm*cos(angle));
		
		} else if(atom_counter > 3) {
			if(!(iss >> bond_atom)){
				file_error_message(coord_file);
				exit(1);
			};
			if(!(iss >> r)){
				file_error_message(coord_file);
				exit(1);
			};
			if(!(iss >> angle_atom)){
				file_error_message(coord_file);
				exit(1);
			};
			if(!(iss >> angle)){
				file_error_message(coord_file);
				exit(1);
			};
			if(!(iss >> di_angle_atom)){
				file_error_message(coord_file);
				exit(1);
			};
			if(!(iss >> di_angle)){
				file_error_message(coord_file);
				exit(1);
			};
	
			if( (cur_atom.get_bond_length_atom() != bond_atom) or  
			(cur_atom.get_angle_atom() != angle_atom) or 
		    (cur_atom.get_dihedral_angle_atom() != di_angle_atom) ){
				Z_coord_error();
				exit(1); 
			};
	
			angle = M_PI - angle;

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
						+ z_axis*cos(angle));
			cur_atom.set_cart_coord(type, bond_coords[0] + added_vec[0], bond_coords[1] + added_vec[1], 
				bond_coords[2] + added_vec[2]);	
		};
	};
	move_to_COM(type);
	return;

};



double molecule::energy() {
	if(!coefficients_set) {
		std::cerr << "Haven't set the coefficients" << std::endl;
		exit(1);
	};
	if(!L_matrix_set) {
		std::cerr << "Haven't set the L matrix" << std::endl;
		exit(1);
	};
	
	Eigen::MatrixXd A = Amat();
	Eigen::MatrixXd T = Tmat(A);	
	std::cout << T << std::endl;
	rotate_coords(T);
	
	Eigen::Vector3d norm_coord = normal_coordinates(); 
	double total_energy = 0; 
	for(int i = 0; i < coeffs_1D.get_modes(); i++) {
		for(int j = 0; j < coeffs_1D.get_order(); j++) {
			total_energy +=  coeffs_1D.get_coefficient_1D(i, j)*pow(norm_coord(i), j+1);	
		};
	};	
	for(int i = 0; i < coeffs_2D.get_modes(); i++) {
		for(int j = 0; j < i; j++) {
			for(int k = 0; k < coeffs_2D.get_order(); k++) {
				for(int l = 0; l < coeffs_2D.get_order(); l++) {
					total_energy += coeffs_2D.get_coefficient_2D(i, j, k, l)*pow(norm_coord(i),k+1)*pow(norm_coord(j),l+1);
				};
			};
		};
	};
	return total_energy; 
};

Eigen::Vector3d molecule::centre_of_mass(int type) {
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
	return com; 
};

void molecule::move_to_COM(int type) {
	Eigen::Vector3d com = centre_of_mass(type);
	for(int i = 1; i <= num_atoms; i++) {
		update_atom_coord(type, i, -com[0], -com[1], -com[2]);
	};
	return; 
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
//	std::cout << coord_1 << std::endl << coord_2 << std::endl << coord_3 << std::endl;
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
	Eigen::Vector3d vec_1 = coord_1 - coord_2;
	Eigen::Vector3d z = coord_2 - coord_3;
	z = z/sqrt(z.dot(z));
	Eigen::Vector3d vec_2 = coord_4 - coord_3;
	Eigen::Vector3d x_1 = (z.cross(vec_1)).cross(z);
	x_1 = x_1/sqrt(x_1.dot(x_1));

	Eigen::Vector3d y = z.cross(vec_2);
	y = y/sqrt(y.dot(y));
	Eigen::Vector3d x_2 = y.cross(z);
	x_2 = x_2/sqrt(x_2.dot(x_2)); 

	double x_component = x_1.dot(x_2);
	double y_component = x_1.dot(y);
	
	double angle = atan2(y_component, x_component);
	if(angle < 0) {
		angle += 2*M_PI;
	};
	return angle;	
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

void molecule::file_error_message(std::string file) {
	std::cerr << "File " << file << " contains an error" << std::endl;	
};

void molecule::coord_length_error_message() {
	std::cerr << "Coordinates incorrect length" << std::endl;
};


void molecule::Z_coord_error() {
	std::cerr << "The inputted Z matrix does not match what is stored in this molecule" << std::endl;
};

std::string molecule::find_line(std::ifstream& stream, int num_words, const std::string& target_sentence) {
	std::string line;
	std::stringstream fss;
	while( getline(stream, line) ){
		std::istringstream iss_1(line);
		std::string word;
		std::stringstream ss;
		for(int i = 0; i < num_words; i++) {
			iss_1 >> word;
			ss << word << " ";
		};
		std::string sentence = ss.str();
		if(sentence == target_sentence) {
			while(iss_1 >> word) {
				fss << word << " ";
			};
			break;
		};
	};
	return fss.str();
};

std::vector<Eigen::MatrixXd> molecule::inverted_design_matrices(const Eigen::MatrixXd& grid_points, int dim, int poly_order) {
	std::vector<Eigen::MatrixXd> matrices;
	for(int i = 0; i < dim; i++) {
		Eigen::MatrixXd design_mat = Eigen::MatrixXd::Zero(grid_points.rows(), poly_order + 1);
		for(int j = 0; j <= poly_order; j++) {
			for(int k = 0; k < grid_points.rows(); k++) {
				design_mat(k, j) = pow(grid_points(k, i), j); 
			};
		};
		Eigen::MatrixXd inverse_mat = design_mat.completeOrthogonalDecomposition().pseudoInverse();
		matrices.push_back(inverse_mat);
	};
	return matrices;
};


Eigen::MatrixXd poly_factor_mat( const Eigen::MatrixXd& coordinates, int poly_order) {
	Eigen::MatrixXd poly_mat = Eigen::MatrixXd::Constant(coordinates.rows(), poly_order + 1, 1);
	for(int i = 1; i < poly_mat.cols(); i++) {
		for(int j = 0; j < poly_mat.rows(); j++) {
			poly_mat(i,j) = poly_mat(i-1,j)*coordinates(j);
		};
	};
	return poly_mat;
};



