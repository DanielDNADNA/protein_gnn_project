#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <utility>
#include <cassert>
#include <stdexcept>

#include "atom.h"
#include "protein_gnn_project.h"
#include "cell_list.h"
#include "kd_tree.h"
#include "evaluation.h"

Evaluation::Evaluation(AtomRecord& atoms, const std::string& fileName) : atoms(atoms) {

    std::cout << "For protein " << fileName << " with " << atoms.x.size() << " atoms: \n";

    cutoff = 5.0f;
    cutoff_sq = cutoff*cutoff;
}

void Evaluation::evaluateAlgorithms (std::initializer_list<Algorithm> algorithms) {

    for (const auto algorithm : algorithms) {
        switch (algorithm) {
            case Algorithm::BruteForce: {
                std::vector<uint16_t> neighbors = getBruteForceNeighbors();
                printBruteForceTime();
                brute_force_neighbors = std::move(neighbors);
                break;
            }

            case Algorithm::KDTree: {
                buildTree();
                std::vector<uint16_t> neighbors = getTreeNeighbors();
                printTreeTime();
                tree_neighbors = std::move(neighbors);
                break; 
            }

            case Algorithm::CellList: {
                buildCellList();
                std::vector<uint16_t> neighbors = getCellNeighbors();
                printCellTime();
                cell_neighbors = std::move(neighbors);
                break;
            }

            case Algorithm::CellListOMP: {
                buildCellList();
                std::vector<uint16_t> neighbors = getOMPCellNeighbors();
                printOMPCellTime();
                OMP_cell_neighbors = std::move(neighbors);
                break;
            }

            case Algorithm::CudaCellList: {
                throw std::runtime_error("CUDA Cell List not implemented yet");
                break;
            }
        }
    }
}

void Evaluation::validateAlgorithms(std::initializer_list<Algorithm> algorithms) {
    if (algorithms.size() < 2) 
        throw std::runtime_error("Must have 2 or more algorithms to validate");

    auto baseline = *(algorithms.begin());

    std::vector<uint16_t> baseline_vector = match_neighbor(baseline);

    size_t i = 1; 

    for (auto next = std::next(algorithms.begin()); next != algorithms.end(); next++) {
        
        std::vector<uint16_t> comparison = match_neighbor(*next);

        if (!validate(baseline_vector, comparison)) {
            std::cerr << "Algorithm index " << i << " mismatches baseline\n"; 
            throw std::runtime_error("Error: mismatch");
        }

        i++;
    }
}

std::vector<uint16_t> Evaluation::match_neighbor(Algorithm algorithm) {
    
    std::vector<uint16_t> neighbors;

    switch (algorithm) {
    
        case Algorithm::BruteForce: 
        neighbors = brute_force_neighbors;
        break;

        case Algorithm::KDTree:
        neighbors = tree_neighbors;
        break;

        case Algorithm::CellList: 
        neighbors = cell_neighbors;
        break;

        case Algorithm::CellListOMP: 
        neighbors = OMP_cell_neighbors;
        break;

        case Algorithm::CudaCellList:
        throw std::runtime_error("CUDA Cell List not implemented yet");
        break;
    }

    if (neighbors.empty()) 
        throw std::runtime_error("The neighbors vector must be filled");

    return neighbors;
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
    
    std::vector<uint16_t> brute_force_neighbors(atoms.x.size());

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
    std::cout << " ";

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

    std::vector<uint16_t> neighbors = cellList->numNeighbors();

    auto end = std::chrono::steady_clock::now();

    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    cell_list_search_time = ms;

    //purely to make pairwise comparison possible; not part of native cell list usage
    std::vector<uint16_t>reorderedNeighbors = cellList->degreesInOriginalOrder(neighbors);
    return reorderedNeighbors;
}

std::vector<uint16_t> Evaluation::getOMPCellNeighbors() {
    if (!cellList) {
        throw std::runtime_error("Cell list not initialized");
    }

    auto start = std::chrono::steady_clock::now();

    std::vector<uint16_t> neighbors = cellList->parallelNumNeighbors();

    auto end = std::chrono::steady_clock::now();

    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    OMP_list_search_time = ms;

    //purely to make pairwise comparison possible; not part of native cell list usage
    std::vector<uint16_t>reorderedNeighbors = cellList->degreesInOriginalOrder(neighbors);
    return reorderedNeighbors;
}

bool Evaluation::validate(const std::vector<uint16_t>& base_neighbor,
    const std::vector<uint16_t>& compare_neighbor) {

    assert (compare_neighbor.size() == base_neighbor.size());

    for (size_t i = 0; i < compare_neighbor.size(); i++) {
        if (compare_neighbor[i] != base_neighbor[i]) {
            std::cout << "Mismatch at index " << i << '\n';
            std::cout << "Base case: " << base_neighbor[i] << "\n";
            std::cout << "Compare case: " << compare_neighbor[i] << "\n";
            return false;
        }
    }
    return true;
}

void Evaluation::printCellTime() {
    std::cout << "Cell list:      ";
    std::cout << "Build: " << cell_list_build_time << " milliseconds    ";
    std::cout << "Search: " << cell_list_search_time << " milliseconds\n";
}

void Evaluation::printOMPCellTime() {
    std::cout << "OMP Cell list  ";
    std::cout << "Build: " << cell_list_build_time << " milliseconds    ";
    std::cout << "Search: " << OMP_list_search_time << " milliseconds\n";
}

void Evaluation::printTreeTime() {
    std::cout << "KD-tree:       ";
    std::cout << "Build: " << tree_build_time << " milliseconds    ";
    std::cout << "Search: " << tree_search_time << " milliseconds\n";
}

void Evaluation::printBruteForceTime() {
    std::cout << "Brute force     ";
    std::cout << brute_force_time << " milliseconds\n";
}
