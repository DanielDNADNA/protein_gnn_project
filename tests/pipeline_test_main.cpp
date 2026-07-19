#include <pybind11/embed.h>
#include "protein_gnn_project.h"
#include "kd_tree.h"
#include "atom.h"
#include <iostream>
#include <vector>
#include <chrono>


namespace py = pybind11;

int main() {
    
    py::scoped_interpreter guard{};
    
    std::string secondsType = "microseconds";
    std::string fileName = "1CRN.cif.gz";

    auto start = std::chrono::high_resolution_clock::now();

    Parse_mmCIF protein{fileName};

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Atoms and Edges struct construction time: " << elapsed.count() << " " << secondsType << "\n";
    
    start = std::chrono::high_resolution_clock::now();

    auto [coordinates, features, targetVector, edgeIndex, edgeAttr] = protein.exportTensors();

    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Tensor construction time: " << elapsed.count() << " " << secondsType << "\n";;

    AtomRecord atoms = protein.atoms;

    float avgDeg = 0; 
    for (size_t i = 0; i< atoms.atom_degree.size(); i++) {
        avgDeg += atoms.atom_degree[i];
    }

    avgDeg /= atoms.atom_degree.size();


    std::cout << "average degree: " << avgDeg << "\n";

    Edges edges = protein.edges;
    
    for (size_t i = 0; i < 10; i++) {
        std::cout << "index " << i << '\n';
        std::cout << "coords : \n";
        std::cout << "x: " << atoms.x[i] << "   ";
        std::cout << "y: " << atoms.y[i] << "   ";
        std::cout << "z: " << atoms.z[i] << "   " << '\n';
        std::cout << "b-factor: " << atoms.b_factor[i] << '\n'; 
    }
    
    return 0;
}