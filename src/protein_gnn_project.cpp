#include <iostream>
#include <vector>
#include <string>
#include <pybind11/pybind11.h>
#include "protein_gnn_project.h"

void helloWorld() {
	std::cout << "Hello World!" << '\n';
	return;
}

PYBIND11_MODULE(spatial_engine, m) {
	m.doc() = "My first spatial engine plugin";
	m.def("hello_world", &helloWorld, "A function that prints to the console");
}
