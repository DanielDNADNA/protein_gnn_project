#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <filesystem>
#include <cstdint>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <gemmi/mmread_gz.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmcif.hpp>
#include "protein_gnn_project.h"

namespace py = pybind11;
namespace fs = std::filesystem;

struct AtomRecord {

	py::ssize_t count;

	std::vector<float> x, y, z; //positions
	
	std::vector<uint8_t> atom_type; 
	std::vector<uint8_t> locant; //relative position in protein
	std::vector<uint8_t> residue_type;

	std::vector<int> atom_index;
	std::vector<int> atom_degree;

	void reserve(py::ssize_t count) {
		size_t n = static_cast<size_t>(count);

		x.reserve(n);
		y.reserve(n);
		z.reserve(n);

		atom_type.reserve(n);
		locant.reserve(n);
		residue_type.reserve(n);

		atom_index.reserve(n);
		atom_degree.reserve(n);
	}
};

class Parse_mmCIF {

public: 

	py::ssize_t num_atoms = 0;
	gemmi::Structure st;

	std::unordered_map<std::string, int> atom_map = {
		{"C", 0}, 
		{"N", 1},
		{"O", 2},
		{"S", 3},
		{"P", 4}
	};
	
	std::unordered_map<std::string, int> residue_map = {
		{"ALA", 0},  {"ARG", 1},  {"ASN", 2},  {"ASP", 3},
		{"CYS", 4},  {"GLU", 5},  {"GLN", 6},  {"GLY", 7},
		{"HIS", 8},  {"ILE", 9},  {"LEU", 10}, {"LYS", 11},
		{"MET", 12}, {"PHE", 13}, {"PRO", 14}, {"SER", 15},
		{"THR", 16}, {"TRP", 17}, {"TYR", 18}, {"VAL", 19},
	};

	AtomRecord atoms;

	Parse_mmCIF(const std::string& filePath) {
		openFile(filePath);
		countAtoms();
	}

	//atom.element.name()
	int getAtomID(const std::string& elem) {
		auto it = atom_map.find(elem);
		if (it != atom_map.end()) return it->second;
		else return 5; //unknown
	}

	//res.name
	int getResidueID(const std::string& res) {
		auto it = residue_map.find(res); 
		if (it != residue_map.end()) return it->second;
		else return 20; //unknown
	}

	int getLocantDepth(const std::string& atom_name) {
		//these are already position 0
		if (atom_name == "N" || atom_name == "C" || atom_name == "O" ||
			atom_name == "OXT" || atom_name == "CA") {
			return 0;
		}
	
		//7 is chosen as "other"
		if (atom_name.length() < 2) return 7;
		if (!std::isalpha(atom_name[1])) return 7;

		char locant = atom_name[1];

		switch (locant) {
			case 'B': return 1; break; //Beta
			case 'G': return 2; break; //Gamma
			case 'D': return 3; break; //Delta
			case 'E': return 4; break; //Epsilon
			case 'Z': return 5; break; //Zeta
			case 'H': return 6; break; //Eta
			default: return 7; break;
		}
	}

	void openFile(const std::string& filePath) {
		st = gemmi::read_structure_gz(filePath, gemmi::CoorFormat::Detect);
	}

	void countAtoms() {
		py::ssize_t count = 0;
		//structure -> model -> chain -> residue -> atom

		for (const gemmi::Model& model : st.models) {
			for (const gemmi::Chain& chain : model.chains) {
				std::string chain_id = chain.name;
				for (const gemmi::Residue& res : chain.residues) {

					if (res.is_water())
						continue;

					for (const gemmi::Atom& atom : res.atoms) {

						if (atom.element == gemmi::El::H)
							continue;

						count++;

					}
				}
			}
		}
		num_atoms = count;
		return;
	}
	
	void constructAtoms() {

		atoms.reserve(num_atoms);
		int i = 0;

		for (const gemmi::Model& model : st.models) {
			for (const gemmi::Chain& chain : model.chains) {

				for (const gemmi::Residue& res : chain.residues) {

					if (res.is_water())
						continue;

					std::string res_name = res.name;
					int res_id = getResidueID(res.name);

					for (const gemmi::Atom& atom : res.atoms) {

						if (atom.element == gemmi::El::H)
							continue;
						
						float x_cord = atom.pos.x;
						float y_cord = atom.pos.y;
						float z_cord = atom.pos.z;

						int atom_type = getAtomID(atom.element.name());
						int atom_locant = getLocantDepth(atom.name);

						//placeholder
						int atom_degree = 0;

						atoms.x.push_back(x_cord);
						atoms.y.push_back(y_cord);
						atoms.z.push_back(z_cord);

						atoms.atom_type.push_back(atom_type);
						atoms.locant.push_back(atom_locant);
						atoms.residue_type.push_back(res_id);

						atoms.atom_index.push_back(i);
						atoms.atom_degree.push_back(atom_degree);

						i++;
					}
				}
			}
		}
	}

	py::array_t<float> getCoords() {

		py::array_t<float> coords(
			py::array::ShapeContainer{ num_atoms, (py::ssize_t)3 }
		);

		auto view = coords.mutable_unchecked<2>();

		for (size_t i = 0; i < atoms.x.size(); i++) {

			view(i, 0) = atoms.x[i];
			view(i, 1) = atoms.y[i];
			view(i, 2) = atoms.z[i];
		}

		return coords;
	}

	py::array_t<float> constructNodes() {

		//0-6: Atom type, 7-14: Atom locant, 15: atom degree, 16-36: residue type
		py::array_t<float> nodes(
			py::array::ShapeContainer{ num_atoms, (py::ssize_t)36 }
		);

		//initialize to 0
		std::fill(nodes.mutable_data(), nodes.mutable_data() + nodes.size(), 0.0f);

		auto view = nodes.mutable_unchecked<2>();

		for (size_t i = 0; i < atoms.x.size(); i++) {
			view(i, atoms.atom_type[i]) = 1.0f;
			view(i, atoms.locant[i] + 6) = 1.0f;
			view(i, 14) = atoms.atom_degree[i];
			view(i, 15 + atoms.residue_type[i]) = 1.0f;
		}

		return nodes;
	}
};

class ConstructEdges {

};


PYBIND11_MODULE(spatial_engine, m) {
	m.doc() = "Parses mmCIF file";
	py::class_<Parse_mmCIF>(m, "Parse_mmCIF")

		.def(py::init<const std::string&>())

		.def_readwrite("atom_map", &Parse_mmCIF::atom_map)

		.def("getAtomID", &Parse_mmCIF::getAtomID)
		.def("getResidueID", &Parse_mmCIF::getResidueID)
		.def("getLocantDepth", &Parse_mmCIF::getLocantDepth)
		.def("openFile", &Parse_mmCIF::openFile)
		.def("countAtoms", &Parse_mmCIF::countAtoms)
		.def("constructAtoms", &Parse_mmCIF::constructAtoms)
		.def("getCoords", &Parse_mmCIF::getCoords)
		.def("constructNodes", &Parse_mmCIF::constructNodes);

}
