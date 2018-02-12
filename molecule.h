#ifndef MOLECULE_H
#define MOLECULE_H

#include "atom.h"
#include <vector>
#include <string> 

class molecule {
	public:
		int num_atoms = 0;
		std::vector<atom> atoms;
		std::string const name;
		atom& get_atom_from_num(unsigned int n);
		
		static void file_error_message();
		static void coord_length_error_message();
		static double dot_product(std::vector<double> coords_1, std::vector<double> coords_2);
		static std::vector<double> displacement_vector(std::vector<double> coords_1, std::vector<double> coords_2);
		static double distance(std::vector<double>  coords_1, std::vector<double> coords_2);
		static double calculate_angle(std::vector<double> coord_1, std::vector<double> coord_2, std::vector<double> coord_3);
		static double calculate_dihedral_angle(std::vector<double>, std::vector<double>, std::vector<double>,
				std::vector<double>);
		static double calculate_vec_angle(std::vector<double>, std::vector<double>);
		static double vec_distance(std::vector<double>);
		static bool check_coord_size(std::vector<double> coord);
		static std::vector<double> cross_product(std::vector<double> coord_1, std::vector<double> coord_2);
		static std::vector<double> vec_normalised(std::vector<double> vec);
	public:
		void print_coordinates();
		int get_num_atoms();
		molecule(std::string z_matrix_file, std::string molecule_name);	
		double bond_length(int atom_num);
		double angle(int atom_num);
		double dihedral_angle(int atom_num);
		void set_atom_coord(int atom_num, double x, double y, double z );
		void set_molecule_coord(std::string coord_file);
};
#endif
