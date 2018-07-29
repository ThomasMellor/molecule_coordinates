#ifndef MOLECULE_H
#define MOLECULE_H

#include "atom.h"
#include <vector>
#include <string> 
#include "eigen-eigen-5a0156e40feb/Eigen/Dense"

class molecule {

	class coefficient {
		public:
			int modes;
			int order;
			int dimension;
			std::vector<std::vector<double>> coeffs; 
		public:
			int get_modes();
			int get_order();
			coefficient(int num_modes, int num_order, int num_dimension);
			coefficient(); 
			void set_coefficient_2D(int mode_1, int mode_2, int order_1, int order_2, int val);
			double get_coefficient_2D(int mode_1, int mode_2, int order_1, int order_2);
			void set_coefficient_1D(int mode, int order, double val); 
			double get_coefficient_1D(int mode, int order);
			void print_coefficients(); 
	};

	public:
		int num_atoms = 0;
		bool coefficients_set = false;
		bool L_matrix_set = false;
		std::vector<atom> atoms;
		std::string const name;
		atom& get_atom_from_num(unsigned int n);
		
		std::vector<double> frequencies;
		Eigen::MatrixXd L_mat;
		Eigen::MatrixXd M_mat;
		molecule::coefficient coeffs_1D;
		molecule::coefficient coeffs_2D;

		static void file_error_message(std::string file);
		static void coord_length_error_message();

		//static double dot_product(Eigen::Vector3d coords_1, Eigen::Vector3d coords_2);
		//static std::vector<double> displacement_vector(std::vector<double> coords_1, std::vector<double> coords_2);
		static double distance(Eigen::Vector3d  coords_1, Eigen::Vector3d coords_2);
		static double calculate_angle(Eigen::Vector3d coord_1, Eigen::Vector3d coord_2, Eigen::Vector3d coord_3);
		static double calculate_dihedral_angle(Eigen::Vector3d coord_1, Eigen::Vector3d coord_2, Eigen::Vector3d coord_3, Eigen::Vector3d coord_4);
		static double calculate_vec_angle(Eigen::Vector3d coord_1, Eigen::Vector3d coord_2);
		static double vec_distance(Eigen::Vector3d coord);
		static bool check_coord_size(std::vector<double> coord);
		static Eigen::Vector3d cross_product(Eigen::Vector3d coord_1, Eigen::Vector3d coord_2);
		//static std::vector<double> vec_normalised(std::vector<double> vec);
		std::vector<Eigen::Vector3d> atom_and_connected_coord(int type, std::vector<int> connected_atoms);
		
		double bond_length_derivative(int atom_num, int second_atom_num, std::string axis);
		double angle_derivative(int atom_num, int second_atom_num, std::string axis);
		double dihedral_angle_derivative(int atom_num, int second_atom_num, std::string axis);
		
		static int axes_name_to_num(std::string axes);
		double coordinate_value(int type, int second_atom_num, std::string axis);
		static double derivative_increment(double q);
		static int connected_atom_pos(std::vector<int> connected_atoms, int second_atom);	
		static double derivative_value(double R_plus, double R_minus, double dq);
		static void print_coords(std::vector<double> coord);
	
		Eigen::MatrixXd empty_matrix();	
		Eigen::MatrixXd Amat();
		static Eigen::MatrixXd Tmat(const Eigen::MatrixXd& A);
		static Eigen::MatrixXd ATAmat(const Eigen::MatrixXd& A);
		void rotate_coords(const Eigen::MatrixXd& T);
		Eigen::Vector3d Eckart_cond();

		static Eigen::MatrixXd& derivative_matrix(Eigen::MatrixXd& mat);
		static Eigen::Vector3d vec_normalised(Eigen::Vector3d vec);
		static bool check_derivative_atoms(std::vector<int> connected_atoms, int atom_1, int atom_2);
		Eigen::VectorXd& derivative_vector(Eigen::VectorXd& vec, std::vector<std::vector<double>> try_coord);
		Eigen::VectorXd empty_vector();
		static void Z_coord_error(); 		
		Eigen::VectorXd normal_coordinates(); 


		double coord_difference(std::vector<std::vector<double>> try_coord);
		void move_to_COM(int type);
		Eigen::Vector3d centre_of_mass(int type); 
	public:
		void print_coordinates(int type);
		void print_cart_coords(int type);
		int get_num_atoms() const;
		molecule(std::string z_matrix_file, std::string molecule_name);	
		void set_L_matrix(std::string L_matrix_file);
		void set_coefficients(std::string coefficient_file);
		double bond_length(int type, int atom_num);
		double angle(int type, int atom_num);
		double dihedral_angle(int type, int atom_num);
		void set_atom_coord(int type, int atom_num, double x, double y, double z );
		Eigen::Vector3d get_atom_coord(int type, int atom_num);
		void update_atom_coord(int type, int atom_num, double dx, double dy, double dz);
		void set_molecule_coord(int type, std::string coord_file);
		void set_molecule_coord_Z(int type, std::string coord_file);
		void update_molecule_coord(const Eigen::VectorXd& vec);
		void print_coefficients_1D(); 
		double energy(); 
};
#endif
