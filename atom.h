#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <memory>

class atom {
	private:
		const unsigned int number;
		const std::string name;
		std::vector<double> cart_coords = {0,0,0};
		std::shared_ptr<atom> bond_length_atom;
		std::vector<std::shared_ptr<atom>> angle_atoms;
		std::vector<std::shared_ptr<atom>> dihedral_angle_atoms;
	public:
		atom(unsigned int num, std::string name);	
		std::string get_name() const; 
		int get_number() const;
		std::vector<double> get_cart_coords() const;
		void set_cart_coords(double x, double y, double z); 
};

#endif
