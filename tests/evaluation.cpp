#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <utility>
#include <cassert>
#include <stdexcept>
#include <nlohmann/json.hpp>
#include <filesystem>

#include "atom.h"
#include "protein_gnn_project.h"
#include "cell_list.h"
#include "kd_tree.h"
#include "evaluation.h"
#include "algorithm.h"

Evaluation::Evaluation(AtomRecord& atoms, const std::string& fileName, 
    std::vector<Algorithm> algorithms) : atoms(atoms), algorithms(algorithms) {

    std::cout << "For protein " << fileName << " with " << atoms.x.size() << " atoms: \n";

    protein_name = fileName;

    cutoff = 5.0f;
    cutoff_sq = cutoff*cutoff;

    for (auto it = algorithms.begin(); it != algorithms.end(); it++) {
        auto algorithm = *it;

        switch(algorithm) {
            case Algorithm::BruteForce: {
                getBruteForceNeighbors();
                break;
            }

            case Algorithm::KDTree: {
                buildTree();
                getTreeNeighbors();
                break;
            }

            case Algorithm::CellList: {
                buildCellList();
                getCellNeighbors();
                break;
            }

            case Algorithm::OMPCellList: {
                buildOMPCellList();
                getOMPCellNeighbors();
                break;
            }

            case Algorithm::CudaCellList: {
                throw std::runtime_error("Cuda in WIP; not available");
                break;
            }

            default: throw std::runtime_error("Invalid algorithm");
                break;
        }
    }
    setAlgorithmSearch();
}

void Evaluation::reportBenchmarks() {
    
    for (auto it = algorithms.begin(); it != algorithms.end(); it++) {
        auto algorithm = *it;

        switch(algorithm) {
            case Algorithm::BruteForce: {
                printBruteForceTime();
                break;
            }

            case Algorithm::KDTree: {
                printTreeTime();
                break;
            }

            case Algorithm::CellList: {
                printCellTime();
                break;
            }

            case Algorithm::OMPCellList: {
                printOMPCellTime();
                break;
            }

            case Algorithm::CudaCellList: {
                throw std::runtime_error("Cuda in WIP; not available");
                break;
            }

            default: throw std::runtime_error("Invalid algorithm");
                break;
        }
    }
}

void Evaluation::validateAlgorithms() {

    if (algorithms.size() < 2) 
        throw std::runtime_error("Must have 2 or more algorithms to validate");

    const auto& baseline = *(algorithms.begin());

    std::vector<uint16_t> baseline_vector = algorithm_search[baseline].neighbors;

    if (baseline_vector.empty())
    throw std::runtime_error("Baseline algorithm was not run");

    for (auto next = std::next(algorithms.begin()); next != algorithms.end(); next++) {
        
        const auto& comparison = algorithm_search[*next].neighbors;

        if (!validate(baseline_vector, comparison)) {
            std::cerr << "Algorithm " << algorithm_search[*next].algorithm_name 
                << " mismatches baseline\n"; 
            throw std::runtime_error("Error: mismatch");
        }
    }
}

//Map between algorithms enum class and the relevant data
void Evaluation::setAlgorithmSearch() {

    //ensure that only valid algorithms become keys in algorithm_search
    for (auto it = algorithms.begin(); it != algorithms.end(); it++) {
        auto algorithm = *it;

        switch(algorithm) {
            case Algorithm::BruteForce: {

                brute_force.algorithm_name = "brute_force";
                brute_force.protein_name = protein_name;
                brute_force.num_atoms = atoms.x.size();

                algorithm_search[Algorithm::BruteForce] = brute_force;
                break;
            }

            case Algorithm::KDTree: {
                kd_tree.algorithm_name = "kd_tree";
                kd_tree.protein_name = protein_name;
                kd_tree.num_atoms = atoms.x.size();
    
                algorithm_search[Algorithm::KDTree] = kd_tree;
                break;
            }

            case Algorithm::CellList: {
                cell_list.algorithm_name = "cell_list";
                cell_list.protein_name = protein_name;
                cell_list.num_atoms = atoms.x.size();
    
                algorithm_search[Algorithm::CellList] = cell_list;
                break;
            }

            case Algorithm::OMPCellList: {
                OMP_cell_list.algorithm_name = "OMP_cell_list";
                OMP_cell_list.protein_name = protein_name;
                OMP_cell_list.num_atoms = atoms.x.size();
    
                algorithm_search[Algorithm::OMPCellList] = OMP_cell_list;
                break;
            }

            case Algorithm::CudaCellList: {
                throw std::runtime_error("Cuda in WIP; not available");
                break;
            }

            default: throw std::runtime_error("Invalid algorithm");
                break;
        }
    }
}

void Evaluation::exportJson() {

    namespace fs = std::filesystem;
    using json = nlohmann::ordered_json;

    //format json file
    if (algorithms.size() == 0)
        throw std::runtime_error("Must have at least one algorithm");
    
    json j;

    const auto& first = *(algorithms.begin());

    j["protein"] = protein_name;
    j["atoms"] = algorithm_search.at(first).num_atoms;

    for (auto it = algorithms.begin(); it != algorithms.end(); it++) {
        auto algorithm = *it;
        //checks for a valid algorithm_search entry
        const auto& entry = algorithm_search.at(algorithm);

        if (entry.num_atoms == 0)
            throw std::runtime_error("Algorithm was not benchmarked");

        j["algorithms"][entry.algorithm_name]["search"] = entry.search_ms;
        if (entry.build_ms != 0)
            j["algorithms"][entry.algorithm_name]["build"] = entry.build_ms;
        }

    //Now export
    fs::path directory{BENCHMARK_STATISTICS_DIR};

    fs::create_directories(directory);

    fs::path file_path = directory / (protein_name + ".json");

    std::ofstream file(file_path);

    if (!file)
    throw std::runtime_error("Could not open JSON output file");

    file << j.dump(4) << '\n';
}

void Evaluation::buildCellList() {

    atoms2 = atoms;
    cellList.reset();
    
    auto start = std::chrono::steady_clock::now();

    cellList.emplace(atoms2, cutoff);

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    cell_list.build_ms = static_cast<double>(ms);
}

//exactly the same as buildCellList; but needed for OMP_cell_list build_ms
void Evaluation::buildOMPCellList() {

    atoms2 = atoms;
    cellList.reset();

    auto start = std::chrono::steady_clock::now();

    cellList.emplace(atoms2, cutoff);

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    OMP_cell_list.build_ms = static_cast<double>(ms);
}

void Evaluation::buildTree() {
    
    auto start = std::chrono::steady_clock::now();

    kdTree.emplace(atoms);

    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    kd_tree.build_ms = static_cast<double>(ms);
}

void Evaluation::getBruteForceNeighbors() {
    
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

    brute_force.search_ms = static_cast<double>(ms);

    brute_force.neighbors = std::move(brute_force_neighbors);
}

void Evaluation::getTreeNeighbors() {

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

    kd_tree.search_ms = static_cast<double>(ms);

    kd_tree.neighbors = tree_neighbors;
}

void Evaluation::getCellNeighbors() {

    if (!cellList) {
        throw std::runtime_error("Cell list not initialized");
    }

    auto start = std::chrono::steady_clock::now();

    std::vector<uint16_t> neighbors = cellList->numNeighbors();

    auto end = std::chrono::steady_clock::now();

    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    cell_list.search_ms = ms;

    //purely to make pairwise comparison possible; not part of native cell list usage
    std::vector<uint16_t>reorderedNeighbors = cellList->degreesInOriginalOrder(neighbors);
    cell_list.neighbors = std::move(reorderedNeighbors);
}

void Evaluation::getOMPCellNeighbors() {
    if (!cellList) {
        throw std::runtime_error("Cell list not initialized");
    }

    auto start = std::chrono::steady_clock::now();

    std::vector<uint16_t> neighbors = cellList->parallelNumNeighbors();

    auto end = std::chrono::steady_clock::now();

    auto ms = std::chrono::duration<double, std::milli>(end-start).count();

    OMP_cell_list.search_ms = ms;

    //purely to make pairwise comparison possible; not part of native cell list usage
    std::vector<uint16_t>reorderedNeighbors = cellList->degreesInOriginalOrder(neighbors);
    OMP_cell_list.neighbors = std::move(reorderedNeighbors);
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
    std::cout << "Build: " << cell_list.build_ms << " milliseconds    ";
    std::cout << "Search: " << cell_list.search_ms << " milliseconds\n";
}

void Evaluation::printOMPCellTime() {
    std::cout << "OMP Cell list  ";
    std::cout << "Build: " << OMP_cell_list.build_ms << " milliseconds    ";
    std::cout << "Search: " << OMP_cell_list.search_ms << " milliseconds\n";
}

void Evaluation::printTreeTime() {
    std::cout << "KD-tree:       ";
    std::cout << "Build: " << kd_tree.build_ms << " milliseconds    ";
    std::cout << "Search: " << kd_tree.search_ms << " milliseconds\n";
}

void Evaluation::printBruteForceTime() {
    std::cout << "Brute force     ";
    std::cout << brute_force.search_ms << " milliseconds\n";
}