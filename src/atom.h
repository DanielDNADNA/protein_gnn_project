#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <cstdint>
#include <cassert>
#include <filesystem>

namespace fs = std::filesystem;

struct AtomRecord {
	
	std::vector<float> x, y, z; //positions
	
	std::vector<uint8_t> atom_type; 
	std::vector<uint8_t> locant; //relative position in protein
	std::vector<uint8_t> residue_type;

	std::vector<uint32_t> atom_index;
	std::vector<uint16_t> atom_degree;

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
	void resize(size_t n) {
		x.resize(n);
		y.resize(n);
		z.resize(n);

		atom_type.resize(n);
		locant.resize(n);
		residue_type.resize(n);

		atom_index.resize(n);
		atom_degree.resize(n);

		b_factor.resize(n);
	}

	void validate() const {

    	size_t n = x.size();

    	assert(y.size() == n);
    	assert(z.size() == n);

    	assert(atom_type.size() == n);
    	assert(locant.size() == n);
    	assert(residue_type.size() == n);

    	assert(atom_index.size() == n);
    	assert(atom_degree.size() == n);

    	assert(b_factor.size() == n);
	}
	
	void save(const fs::path& filepath) const {

		//write bytes into the file
		std::ofstream out(filepath, std::ios::binary);
		if (!out) {
            throw std::runtime_error("Cannot open file for writing");
        }

		size_t n = x.size();

		//interprets the file output as bytes (currently; size_t)
		out.write(reinterpret_cast<const char*>(&n), sizeof(n));

		auto writeVector = [&](const auto& vec) {

            out.write(
                reinterpret_cast<const char*>(vec.data()),
				//how many bytes to write
				vec.size() * sizeof(vec[0])
            );
		};

		for (size_t i = 0; i < atom_index.size(); i++) {
    		if (atom_index[i] != i) {
        		throw std::runtime_error("Attempting to save reordered AtomRecord");
    		}
		}	

		writeVector(x);
        writeVector(y);
        writeVector(z);

        writeVector(atom_type);
        writeVector(locant);
        writeVector(residue_type);

        writeVector(atom_index);
        writeVector(atom_degree);

        writeVector(b_factor);
	}

	AtomRecord load(const std::string& filename) {

    	AtomRecord atoms;
    	std::ifstream in(filename, std::ios::binary);

    	if (!in) {
        	throw std::runtime_error("Cannot open file for reading");
    	}

    	size_t n;

    	in.read(reinterpret_cast<char*>(&n), sizeof(n));

    	atoms.resize(n);

    	auto readVector = [&](auto& vec) {
        	in.read(reinterpret_cast<char*>(vec.data()), vec.size() * sizeof(vec[0]));
    	};

    	readVector(atoms.x);
    	readVector(atoms.y);
    	readVector(atoms.z);

    	readVector(atoms.atom_type);
    	readVector(atoms.locant);
    	readVector(atoms.residue_type);

    	readVector(atoms.atom_index);
    	readVector(atoms.atom_degree);

    	readVector(atoms.b_factor);

    	return atoms;
	}
	
};