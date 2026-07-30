#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <optional>

#include "atom.h"
#include "kd_tree.cpp"
#include "cell_list.cpp"

class Evaluation {

    public:

    Evaluation::Evaluation(AtomRecord& atoms, const std::string& fileName);

    void buildCellList();
    void buildTree();

    std::vector<uint16_t> getBruteForceNeighbors();
    std::vector<uint16_t> getTreeNeighbors();
    std::vector<uint16_t> getCellNeighbors();

    void printBruteForceTime();
    void printTreeTime();
    void printCellTime();

    bool validate(const std::vector<uint16_t>& base_neighbor, 
        const std::vector<uint16_t>& compare_neighbor);  

    AtomRecord& atoms; 
    AtomRecord atoms2;
    std::optional<KDTree> kdTree;
    double tree_build_time;
    double tree_search_time;

    std::optional<CellList> cellList;
    double cell_list_build_time;
    double cell_list_search_time;

    double brute_force_time;
    
    private: 

    float cutoff;
    float cutoff_sq;
};