#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <filesystem>
#include <cstdint>
#include <string_view>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <gemmi/mmread_gz.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/chemcomp.hpp>
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

	std::vector<float> b_factor;

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

		b_factor.reserve(n);
	}
};

struct Edges {

	py::ssize_t count;
	std::vector<int> sourceNode; 
	std::vector<int> sinkNode;
	std::vector<int> edgeType;

	void reserve(py::ssize_t count) {
		size_t n = static_cast<size_t>(count);
		sourceNode.reserve(3 * n);
		sinkNode.reserve(3 * n);
		edgeType.reserve(3 * n);
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
	Edges edges;

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
		edges.reserve(num_atoms);

		int i = 0;

		//store atom name and its index
		std::unordered_map<std::string, int> local_atoms;
		//reverse
		std::unordered_map<int, std::string> atoms_from_index;
		local_atoms.reserve(64);
		atoms_from_index.reserve(64);

		std::vector<int>indices;
		indices.reserve(local_atoms.size());

		for (const gemmi::Model& model : st.models) {
			for (const gemmi::Chain& chain : model.chains) {
				
				int prev_C_global_index = -1;

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

						std::string atom_name = atom.name;
						int atom_locant = getLocantDepth(atom_name);

						//placeholder
						int atom_degree = 0;

						atoms.x.push_back(x_cord);
						atoms.y.push_back(y_cord);
						atoms.z.push_back(z_cord);

						atoms.atom_type.push_back(atom_type);
						atoms.locant.push_back(atom_locant);
						atoms.residue_type.push_back(res_id);

						atoms.atom_index.push_back(i);
						local_atoms[atom_name] = i;
						atoms_from_index[i] = atom_name;

						atoms.atom_degree.push_back(atom_degree);

						atoms.b_factor.push_back(atom.b_iso);

						i++;
					}

					//checks and grabs the N index, to be connected with the C of the previous peptide
					if (prev_C_global_index != -1 && local_atoms.count("N") > 0) {
						int current_N_index = local_atoms["N"];
						constructInterResidueEdges(prev_C_global_index, current_N_index);
					}

					constructIntraResidueEdges(local_atoms, atoms_from_index, res_name, indices);

					if (local_atoms.count("C") > 0) {
						prev_C_global_index = local_atoms["C"];
					}
					else {
						prev_C_global_index = -1; 
					}

					local_atoms.clear();
					atoms_from_index.clear();

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

		for (py::ssize_t i = 0; i < atoms.x.size(); i++) {
			view(i, atoms.atom_type[i]) = 1.0f;
			view(i, atoms.locant[i] + 6) = 1.0f;
			view(i, 14) = atoms.atom_degree[i];
			view(i, 15 + atoms.residue_type[i]) = 1.0f;
		}

		return nodes;
	}

	py::array_t<float> constructTargetVector() {
		py::array_t<float>targetNode(
			py::array::ShapeContainer{num_atoms, (py::ssize_t)1 }
		);

		auto view = targetNode.mutable_unchecked<2>();

		for (py::ssize_t i = 0; i < atoms.b_factor.size(); i++) {
			view(i, 0) = atoms.b_factor[i];
		}

		return targetNode;
	}

	//use of AI to find all the rules
	// Returns: 0 (Single), 1 (Double), 2 (Triple), 3 (Aromatic/Delocalized)
	int getBondType(const std::string& res_name, const std::string& atom1, const std::string& atom2) {

		// Lambda to check pairs regardless of order (e.g., C-O is the same as O-C)
		auto is_pair = [&](std::string_view a, std::string_view b) {
			return (atom1 == a && atom2 == b) || (atom1 == b && atom2 == a);
			};

		if (is_pair("C", "O")) return 1;

		// --- 2. DELOCALIZED / AROMATIC RINGS ---
		if (res_name == "PHE" || res_name == "TYR") {
			if (is_pair("CG", "CD1") || is_pair("CD1", "CE1") || is_pair("CE1", "CZ") ||
				is_pair("CZ", "CE2") || is_pair("CE2", "CD2") || is_pair("CD2", "CG")) {
				return 3;
			}
		}
		
		// --- 3. SPECIFIC SIDECHAIN DOUBLE BONDS ---
		if (res_name == "ASN" && is_pair("CG", "OD1")) return 1;
		if (res_name == "GLN" && is_pair("CD", "OE1")) return 1;

		if (res_name == "ARG" && (is_pair("CZ", "NH1") || is_pair("CZ", "NH2"))) return 3;
		if (res_name == "ASP" && (is_pair("CG", "OD1") || is_pair("CG", "OD2"))) return 3;
		if (res_name == "GLU" && (is_pair("CD", "OE1") || is_pair("CD", "OE2"))) return 3;

		// If it's not in the list, it's a standard single bond.
		return 0;
	}

	void constructIntraResidueEdges(const std::unordered_map<std::string, int>& local_atoms, 
		const std::unordered_map<int, std::string>& atoms_from_index, const std::string& res_name, std::vector<int>& indices) {
		indices.clear();
		
		for (const auto& pair : local_atoms) {
			indices.push_back(pair.second);
		}

		for (size_t i = 0; i < indices.size(); i++) {
			for (size_t j = i + 1; j < indices.size(); j++) {
				int src = indices[i];
				int sink = indices[j];

				float dx = atoms.x[src] - atoms.x[sink];
				float dy = atoms.y[src] - atoms.y[sink];
				float dz = atoms.z[src] - atoms.z[sink];
				float dist_sq = (dx * dx) + (dy * dy) + (dz * dz);

				//assumes that covalent bonds do not extend past 1.9 A
				if (dist_sq < 3.61f) {

					std::string atom1_name = atoms_from_index.at(src);
					std::string atom2_name = atoms_from_index.at(sink);

					int edge_type = getBondType(res_name, atom1_name, atom2_name);

					edges.sourceNode.push_back(src);
					edges.sinkNode.push_back(sink);
					edges.edgeType.push_back(edge_type);

					edges.sourceNode.push_back(sink);
					edges.sinkNode.push_back(src);
					edges.edgeType.push_back(edge_type);
				}
			}
		}
		return; 
	}

	void constructInterResidueEdges(int prev_C_global_index, int current_N_index) {

		edges.sourceNode.push_back(prev_C_global_index);
		edges.sinkNode.push_back(current_N_index);
		edges.edgeType.push_back(0);

		edges.sourceNode.push_back(current_N_index);
		edges.sinkNode.push_back(prev_C_global_index);
		edges.edgeType.push_back(0);

		return; 
	}
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
