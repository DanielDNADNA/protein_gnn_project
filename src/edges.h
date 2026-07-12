#pragma once
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

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