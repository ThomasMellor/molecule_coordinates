#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>

class atom {
	private:	
		const unsigned int number;
		const std::string name;
		std::vector<double> cart_coord = {0,0,0};
		int bond_length_atom;
		int angle_atom;
		int dihedral_angle_atom;
	public:
		atom(unsigned int num, std::string name);	
		std::string get_name() const; 
		int get_number();
		std::vector<double> get_cart_coord() const;
		void set_cart_coord(double x, double y, double z); 
		void set_bond_length_atom(int connected_to);
		void set_angle_atom(int atom_1);
		void set_dihedral_angle_atom(int atom_1);
		int get_bond_length_atom();
		int get_angle_atom();
		int get_dihedral_angle_atom();	
};

#endif
