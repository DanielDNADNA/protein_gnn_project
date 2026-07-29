#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <utility>
#include <cassert>

#include "atom.h"
#include "protein_gnn_project.h"
#include "cell_list.h"
#include "kd_tree.h"
#include "evaluation.h"

Evaluation::Evaluation(AtomRecord& atoms, const std::string& fileName) {

    std::cout << "For protein " << fileName << " with " << atoms.x.size() << " atoms: \n";
    assert(atoms.x.size() == atoms.y.size());
    assert(atoms.y.size() == atoms.z.size());

    //used for bruteforce
    atoms2 = atoms;

    cutoff = 5.0f;
    cutoff_sq = cutoff*cutoff;
}

void Evaluation::buildCellList(AtomRecord& atoms2) {
    auto start = std::chrono::steady_clock::now();

    cellList.emplace(atoms2, cutoff);

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << "Cell list build " << ms << " milliseconds\n";
}

void Evaluation::buildTree(AtomRecord& atoms2) {
    
    auto start = std::chrono::steady_clock::now();

    kdTree.emplace(atoms2);

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << "KD tree build " << ms << " milliseconds\n";
}

void Evaluation::getBruteForceNeighbors() {
    
    auto start = std::chrono::steady_clock::now();
    
    brute_force_neighbors.resize(atoms2.x.size());

    for (size_t i = 0; i < atoms2.x.size(); i++) {
        float x = atoms2.x[i];
        float y = atoms2.y[i];
        float z = atoms2.z[i];

        for (size_t j = 0; j < atoms2.x.size(); j++) {
            float dx = atoms2.x[j] - x;
            float dy = atoms2.y[j] - y;
            float dz = atoms2.z[j] - z;
            if (dx*dx + dy*dy + dz*dz <= cutoff_sq) {
                brute_force_neighbors[i]++;
            }
        }
    }

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << "Brute force search: " << ms << " milliseconds\n";
}

void Evaluation::getTreeNeighbors() {
    
    auto start = std::chrono::steady_clock::now();

    if (!kdTree) {
        throw std::runtime_error("KDTree not initialized");
    }

    tree_neighbors.clear();
    tree_neighbors.resize(atoms2.x.size());

    std::vector<int> neighbors;
	
    neighbors.reserve(64);

	for (size_t i = 0; i < atoms2.x.size(); i++) {

		int nodeDegree = kdTree->getNodeDegree (atoms2.x[i], atoms2.y[i], atoms2.z[i], 5.0f, 
			neighbors);

		tree_neighbors[i] = (static_cast<uint16_t>(nodeDegree));
	}

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << "KD-tree search: " << ms << " milliseconds\n";
}

void Evaluation::getCellNeighbors() {

    auto start = std::chrono::steady_clock::now();

    if (!cellList) {
        throw std::runtime_error("Cell list not initialized");
    }

    reorderedNeighbors = cellList->numNeighbors();

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << "Cell list search: " << ms << " milliseconds\n";

    //purely to make pairwise comparison possible; not part of native cell list usage
    reorderedNeighbors = cellList->degreesInOriginalOrder();
}

bool Evaluation::validate(const std::vector<uint16_t>& base_neighbor,
    const std::vector<uint16_t>& compare_neighbor) {

    assert (compare_neighbor.size() == base_neighbor.size());

    for (size_t i = 0; i < compare_neighbor.size(); i++) {
        if (compare_neighbor[i] != base_neighbor[i]) {
            std::cerr << "Mismatch at index " << i << '\n';
            return false;
        }
    }
    return true;
}
