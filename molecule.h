#ifndef MOLECULE_H
#define MOLECULE_H

#include "atom.h"
#include <vector>
#include <string> 

class molecule {
	private:
		int num_atoms;
		std::vector<atom> atoms;
		std::string const name;
		atom const& get_atom_from_num(unsigned int n);
	public:
		int get_num_atoms();
		molecule(std::string z_matrix_file, std::string molecule_name);	
};

void file_error_message();

#endif
