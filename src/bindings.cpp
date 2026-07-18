#include <pybind11/pybind11.h>

#include "protein_gnn_project.h"
#include "kd_tree.h"

namespace py = pybind11;

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
		.def("constructNodeFeatures", &Parse_mmCIF::constructNodeFeatures);

    
}

