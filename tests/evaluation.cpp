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


Evaluation::Evaluation(AtomRecord& atoms, const std::string& fileName) : atoms(atoms) {

    std::cout << "For protein " << fileName << " with " << atoms.x.size() << " atoms: \n";
    
    assert(atoms.x.size() == atoms.y.size());
    assert(atoms.y.size() == atoms.z.size());

    cutoff = 5.0f;
    cutoff_sq = cutoff*cutoff;
}

void Evaluation::buildCellList() {

    atoms2 = atoms;

    auto start = std::chrono::steady_clock::now();

    cellList.emplace(atoms2, cutoff);

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    cell_list_build_time = static_cast<double>(ms);
}

void Evaluation::buildTree() {
    
    auto start = std::chrono::steady_clock::now();

    kdTree.emplace(atoms);

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    tree_build_time = static_cast<double>(ms);
}

std::vector<uint16_t> Evaluation::getBruteForceNeighbors() {
    
    std::vector<uint16_t>brute_force_neighbors(atoms.x.size());

    auto start = std::chrono::steady_clock::now();

    for (size_t i = 0; i < atoms.x.size(); i++) {
        float x = atoms.x[i];
        float y = atoms.y[i];
        float z = atoms.z[i];

        for (size_t j = 0; j < atoms.x.size(); j++) {
            float dx = atoms.x[j] - x;
            float dy = atoms.y[j] - y;
            float dz = atoms.z[j] - z;
            if (dx*dx + dy*dy + dz*dz <= cutoff_sq) {
                brute_force_neighbors[i]++;
            }
        }
    }

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    brute_force_time = static_cast<double>(ms);

    return brute_force_neighbors;
}

std::vector<uint16_t> Evaluation::getTreeNeighbors() {

    if (!kdTree) {
        throw std::runtime_error("KDTree not initialized");
    }

    auto start = std::chrono::steady_clock::now();

    std::vector<uint16_t> tree_neighbors(atoms.x.size());

    std::vector<int> neighbors;
	
    neighbors.reserve(64);

	for (size_t i = 0; i < atoms.x.size(); i++) {

		int nodeDegree = kdTree->getNodeDegree (atoms.x[i], atoms.y[i], atoms.z[i], 5.0f, 
			neighbors);

		tree_neighbors[i] = (static_cast<uint16_t>(nodeDegree));
	}

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    tree_search_time = static_cast<double>(ms);

    return tree_neighbors;
}

std::vector<uint16_t> Evaluation::getCellNeighbors() {

    if (!cellList) {
        throw std::runtime_error("Cell list not initialized");
    }

    auto start = std::chrono::steady_clock::now();

    std::vector<uint16_t> reorderedNeighbors = cellList->numNeighbors();

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    //purely to make pairwise comparison possible; not part of native cell list usage
    reorderedNeighbors = cellList->degreesInOriginalOrder();
    return reorderedNeighbors;
}

bool Evaluation::validate(const std::vector<uint16_t>& base_neighbor,
    const std::vector<uint16_t>& compare_neighbor) {

    assert (compare_neighbor.size() == base_neighbor.size());

    for (size_t i = 0; i < compare_neighbor.size(); i++) {
        if (compare_neighbor[i] != base_neighbor[i]) {
            std::cout << "Mismatch at index " << i << '\n';
            return false;
        }
    }
    return true;
}

void Evaluation::printCellTime() {

    std::cout << "Cell list build: " << cell_list_build_time << " milliseconds\n";
    std::cout << "Cell list search: " << cell_list_search_time << " milliseconds\n";
}

void Evaluation::printTreeTime() {
    std::cout << "KD-tree build: " << tree_build_time << " milliseconds\n";
    std::cout << "KD-tree search: " << tree_search_time << " milliseconds\n";
}

void Evaluation::printBruteForceTime() {
    std::cout << "Brute force search: " << brute_force_time << "milliseconds\n";
}

