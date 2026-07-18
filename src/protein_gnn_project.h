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

	std::unordered_map<std::string, int> atom_map;
	
	std::unordered_map<std::string, int> residue_map;

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

};