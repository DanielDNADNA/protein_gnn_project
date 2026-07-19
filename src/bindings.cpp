#include <pybind11/pybind11.h>

#include "protein_gnn_project.h"
#include "kd_tree.h"

namespace py = pybind11;

PYBIND11_MODULE(spatial_engine, m) {
	m.doc() = "Parses mmCIF file";
	py::class_<Parse_mmCIF>(m, "Parse_mmCIF")

		.def(py::init<const std::string&>())
		.def("exportTensors", &Parse_mmCIF::exportTensors);
}

