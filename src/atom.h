#pragma once
#include <vector>
#include <cstdint>

struct AtomRecord {
	
	std::vector<float> x, y, z; //positions
	
	std::vector<uint8_t> atom_type; 
	std::vector<uint8_t> locant; //relative position in protein
	std::vector<uint8_t> residue_type;

	std::vector<size_t> atom_index;
	std::vector<int> atom_degree;

	std::vector<float> b_factor;

	void reserve(size_t n) {

		x.reserve(n);
		y.reserve(n);
		z.reserve(n);

		atom_type.reserve(n);
		locant.reserve(n);
		residue_type.reserve(n);

		atom_index.reserve(n);
		atom_degree.reserve(n);

		b_factor.reserve(n);
	}
};