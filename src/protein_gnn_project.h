#pragma once
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <unordered_map>
#include <string>
#include <tuple>
#include <gemmi/mmread_gz.hpp>
#include <gemmi/model.hpp>
#include "atom.h"
#include "edges.h"

namespace py = pybind11;

class Parse_mmCIF {

public: 

	py::ssize_t num_atoms;
	float radius = 5.0f;
	gemmi::Structure st;
	int count = 0;

	std::unordered_map<std::string, int> atom_map = 
	{{"C", 0}, {"N", 1}, {"O", 2}, {"S", 3}, {"P", 4}};
	
	std::unordered_map<std::string, int> residue_map = {{"ALA", 0}, {"ARG", 1}, {"ASN", 2}, {"ASP", 3},
        {"CYS", 4}, {"GLN", 5}, {"GLU", 6}, {"GLY", 7}, {"HIS", 8}, {"ILE", 9}, {"LEU", 10}, {"LYS", 11},
        {"MET", 12}, {"PHE", 13}, {"PRO", 14}, {"SER", 15}, {"THR", 16}, {"TRP", 17}, {"TYR", 18}, {"VAL", 19}};

	AtomRecord atoms;
	Edges edges;

	//fills in the edges and atoms structs
	Parse_mmCIF(const std::string& filePath);

	//returns: coordinates, node features, target vector
	std::tuple<py::array_t<float>, py::array_t<float>, py::array_t<float>, 
		py::array_t<int64_t>, py::array_t<float>> exportTensors();

	std::tuple<py::array_t<int64_t>, py::array_t<float>> constructEdges();

	//atom.element.name()
	int getAtomID(const std::string& elem);

	//res.name
	int getResidueID(const std::string& res);

	int getLocantDepth(const std::string& atom_name);

	void openFile(const std::string& filePath);

	void countAtoms();

	void constructAtoms();

	void getNodeDegree();

	py::array_t<float> getCoords();
    
	py::array_t<float> constructNodeFeatures();
	//use of AI to find all the rules
	// Returns: 0 (Single), 1 (Double), 2 (Triple), 3 (Aromatic/Delocalized)

    py::array_t<float> constructTargetVector();
	int getBondType(const std::string& res_name, const std::string& atom1, 
        const std::string& atom2);

	void constructIntraResidueEdges(const std::unordered_map<std::string, int>& local_atoms, 
		const std::unordered_map<int, std::string>& atoms_from_index, const std::string& res_name, 
        std::vector<int>& indices);

	void constructInterResidueEdges(int prev_C_global_index, int current_N_index);

	void degree();
	
};