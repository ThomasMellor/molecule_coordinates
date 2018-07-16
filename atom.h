#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include "eigen-eigen-5a0156e40feb/Eigen/Dense"

class atom {
	private:	
		const int number;
		const std::string name;
		const double mass;
		Eigen::Vector3d cart_coord = Eigen::Vector3d(0,0,0);
		Eigen::Vector3d eq_coord = Eigen::Vector3d(0,0,0);

		const int bond_length_atom;
		const int angle_atom;
		const int dihedral_angle_atom;
		static double check_sign(double val);
		static int check_value(int atom_num);
	public:
		atom(int num, double m, std::string name);	
		atom(int num, double m, std::string name, int second_atom);
		atom(int num, double m, std::string name, int second_atom,
			int third_atom);
		atom(int num, double m, std::string name, int second_atom,
			int third_atom, int fourth_atom);
		 
		std::string get_name() const; 
		int get_number() const;
		double get_mass() const;
		Eigen::Vector3d get_cart_coord(int type) const;
		void set_cart_coord(int type, double x, double y, double z); 
		void update_cart_coord(int type, double dx, double dy, double dz);
		int get_bond_length_atom() const;
		int get_angle_atom() const;
		int get_dihedral_angle_atom() const;	
		std::vector<int> get_connected_atoms() const;
};

#endif
