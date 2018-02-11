#include "atom.h"
#include <iostream>

int main() {
	atom c_atom(2,"C");
	{	
		atom h_atom(1,"H");
		h_atom.set_bond_length_atom(c_atom);
		atom h_atom_2(h_atom);
		std::cout << h_atom_2.bond_length_atom -> get_name() << std::endl;		
	}
	std::cout << c_atom.get_name() << std::endl;
};
