
#include "protein_gnn_project.h"
#include "kd_tree.h"
#include "atom.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <pybind11/embed.h>

namespace py = pybind11;

int main() {

    py::scoped_interpreter guard{};
    
    std::string secondsType = "microseconds";
    std::string fileName = "3J3Q.cif.gz";

    Parse_mmCIF protein{fileName};

    auto [coordinates, features, targetVector, edgeIndex, edgeAttr] = protein.exportTensors();

    AtomRecord atoms = protein.atoms;

    float avgDeg = 0; 
    for (size_t i = 0; i< atoms.atom_degree.size(); i++) {
        avgDeg += atoms.atom_degree[i];
    }

    return 0;
}

void bruteForceSearch(float qx, float qy, float qz, float radius, 
    const AtomRecord& atoms, std::vector<int>& neighbors) {
    
    neighbors.clear();
    
    for (size_t j = 0; j < atoms.x.size(); j++) {
        float dx = qx - atoms.x[j];
        float dy = qy - atoms.y[j];
        float dz = qz - atoms.z[j];

        if (dx * dx + dy * dy + dz * dz <= radius * radius) {
            neighbors.push_back(j);
        }
    }

    return;
}