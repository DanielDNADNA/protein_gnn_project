#pragma once
#include <vector>

#include "atom.h"

struct Neighbor {
    int atom_index;
    float distance;
};

class KDTree {
    public: 
    
    std::vector<int> atom_index;
    std::vector<int> left; 
    std::vector<int> right;
    std::vector<int> indices;

    const AtomRecord& atoms;

    KDTree(const AtomRecord& atoms);

    bool less (int lhs, int rhs, int axis) const;

    size_t constructTree(int begin, int end, int depth);
    
    void radiusSearch(float qx, float qy, float qz, float radius, int depth, 
        std::vector<int>& neighbors, int index);
};

