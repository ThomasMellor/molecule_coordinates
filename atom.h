#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>

class atom {
	private:	
		const int number;
		const std::string name;
		std::vector<double> cart_coord = {0,0,0};
		const int bond_length_atom;
		const int angle_atom;
		const int dihedral_angle_atom;
		static int check_sign(int atom_num);
		static int check_value(int atom_num);
	public:
		atom(int num, std::string name);	
		atom(int num, std::string name, int second_atom);
		atom(int num, std::string name, int second_atom,
			int third_atom);
		atom(int num, std::string name, int second_atom,
			int third_atom, int fourth_atom);
		 
		std::string get_name() const; 
		int get_number() const;
		std::vector<double> get_cart_coord() const;
		void set_cart_coord(double x, double y, double z); 
		void update_cart_coord(double dx, double dy, double dz);
		int get_bond_length_atom() const;
		int get_angle_atom() const;
		int get_dihedral_angle_atom() const;	
		std::vector<int> get_connected_atoms() const;
};

#endif
