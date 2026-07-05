#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <filesystem>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <gemmi/mmread_gz.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmcif.hpp>
#include "protein_gnn_project.h"

namespace py = pybind11;
namespace fs = std::filesystem;

class Parse_mmCIF {

public: 

	py::ssize_t num_atoms = 0;
	gemmi::Structure st;

	std::unordered_map<std::string, int> atom_type_id;
	

	Parse_mmCIF(const std::string& filePath) {
		openFile(filePath);
		countAtoms();
	}

	void openFile(const std::string& filePath) {
		st = gemmi::read_structure_gz(filePath, gemmi::CoorFormat::Detect);
	}

	void countAtoms() {
		py::ssize_t count = 0;
		//structure -> model -> chain -> residue -> atom

		for (const gemmi::Model& model : st.models) {
			for (const gemmi::Chain& chain : model.chains) {
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
	
	py::array_t<float> getCoords() {

		py::array_t<float> coords(
			py::array::ShapeContainer{ num_atoms, (py::ssize_t)3 }
		);

		py::ssize_t i = 0;
		auto view = coords.mutable_unchecked<2>();

		for (const gemmi::Model& model : st.models) {
			for (const gemmi::Chain& chain : model.chains) {
				for (const gemmi::Residue& res : chain.residues) {

					if (res.is_water())
						continue;

					for (const gemmi::Atom& atom : res.atoms) {

						if (atom.element == gemmi::El::H)
							continue;

						view(i, 0) = static_cast<float>(atom.pos.x);
						view(i, 1) = static_cast<float>(atom.pos.y);
						view(i, 2) = static_cast<float>(atom.pos.z);
						i++;
					}
				}
			}
		}

		return coords;
	}

	void assignIDs() {

	}
};




PYBIND11_MODULE(spatial_engine, m) {
	m.doc() = "Parses mmCIF file";
	py::class_<Parse_mmCIF>(m, "Parse_mmCIF")


		.def(py::init<const std::string&>())

		.def_readwrite("atom_type_id", &Parse_mmCIF::atom_type_id)
		
		.def("openFile", &Parse_mmCIF::openFile)
		.def("countAtoms", &Parse_mmCIF::countAtoms)
		.def("getCoords", &Parse_mmCIF::getCoords)
		.def("assignIDs", &Parse_mmCIF::assignIDs);

}
