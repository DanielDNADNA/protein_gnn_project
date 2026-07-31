#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <optional>
#include <initializer_list>

#include "atom.h"
#include "kd_tree.h"
#include "cell_list.h"


enum class Algorithm {
        BruteForce, 
        KDTree, 
        CellList, 
        CellListOMP, 
        CudaCellList
    };


class Evaluation {

    public:

    Evaluation(AtomRecord& atoms, const std::string& fileName);

    void evaluateAlgorithms(std::initializer_list<Algorithm> algorithms);

    void validateAlgorithms(std::initializer_list<Algorithm> algorithms);
    
    private: 

    float cutoff;
    float cutoff_sq;

    void buildCellList();
    void buildTree();

    std::vector<uint16_t> getBruteForceNeighbors();
    std::vector<uint16_t> getTreeNeighbors();
    std::vector<uint16_t> getCellNeighbors();
    std::vector<uint16_t> getOMPCellNeighbors();

    void printBruteForceTime();
    void printTreeTime();
    void printCellTime();
    void printOMPCellTime();

    bool validate(const std::vector<uint16_t>& base_neighbor, 
        const std::vector<uint16_t>& compare_neighbor);  
    
    std::vector<uint16_t> match_neighbor(Algorithm algorithm);

    AtomRecord& atoms; 
    AtomRecord atoms2;

    std::optional<KDTree> kdTree;
    double tree_build_time;
    double tree_search_time;
    std::vector<uint16_t> tree_neighbors;

    std::optional<CellList> cellList;
    double cell_list_build_time;
    double cell_list_search_time;
    std::vector<uint16_t> cell_neighbors;

    //build method is identical for OMPCellList and normal celllist
    double OMP_list_search_time;
    std::vector<uint16_t> OMP_cell_neighbors;

    double brute_force_time;
    std::vector<uint16_t> brute_force_neighbors;
};