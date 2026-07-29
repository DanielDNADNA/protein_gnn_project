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

    AtomRecord atoms2; //for use by the bruteforce search & KD-tree
    std::vector<uint16_t> reorderedNeighbors;
    std::vector<uint16_t> tree_neighbors;
    std::vector<uint16_t> brute_force_neighbors;
    
    std::optional<KDTree> kdTree;
    std::optional<CellList> cellList;

    Evaluation::Evaluation(AtomRecord& atoms, const std::string& fileName);

    void buildCellList(AtomRecord& atoms2);
    void buildTree(AtomRecord& atoms2);

    void getBruteForceNeighbors();
    void getTreeNeighbors();
    void getCellNeighbors();

    bool validate(const std::vector<uint16_t>& base_neighbor, const std::vector<uint16_t>& compare_neighbor);  
    
    private: 
    float cutoff;
    float cutoff_sq;
};