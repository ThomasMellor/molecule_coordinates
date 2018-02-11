#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>

class atom {
	private:
		const unsigned int number;
		const std::string name;
		std::vector<double> cart_coords = {0,0,0};
		const atom *bond_length_atom;
		std::vector<const atom*> angle_atoms;
		std::vector<const atom*> dihedral_angle_atoms;
	public:
		atom(unsigned int num, std::string name);	
		std::string get_name() const; 
		int get_number() const;
		std::vector<double> get_cart_coords() const;
		void set_cart_coords(double x, double y, double z); 
		void set_bond_length_atom(atom const& connected_to);
		void set_angle_atoms(atom const& atom_1, atom const& atom_2);
		void set_dihedral_angle_atoms(atom const& atom_1, atom const& atom_2, atom const& atom_3);
};

#endif
