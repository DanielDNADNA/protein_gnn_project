#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <optional>
#include <initializer_list>
#include <map>

#include "atom.h"
#include "kd_tree.h"
#include "cell_list.h"
#include "algorithm.h"
#include "algorithm_data.h"

class Evaluation {

    public:

    Evaluation(AtomRecord& atoms, const std::string& fileName, 
        std::vector<Algorithm> algorithms);

    void validateAlgorithms();

    void reportBenchmarks();

    void exportJson();

    private: 

    AlgorithmDetails cell_list;
    AlgorithmDetails OMP_cell_list;
    AlgorithmDetails kd_tree;
    AlgorithmDetails brute_force;

    void buildCellList();
    void buildOMPCellList();
    void buildTree();

    void getBruteForceNeighbors();
    void getTreeNeighbors();
    void getCellNeighbors();
    void getOMPCellNeighbors();

    void printBruteForceTime();
    void printTreeTime();
    void printCellTime();
    void printOMPCellTime();

    bool validate(const std::vector<uint16_t>& base_neighbor, 
        const std::vector<uint16_t>& compare_neighbor); 
        
    void setAlgorithmSearch();

    std::map<Algorithm, AlgorithmDetails> algorithm_search;

    AtomRecord& atoms; 
    AtomRecord atoms2;

    std::optional<KDTree> kdTree;

    std::optional<CellList> cellList;

    std::string protein_name;

    std::vector<Algorithm> algorithms;

    float cutoff;
    float cutoff_sq;
};