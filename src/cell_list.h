#pragma once

#include <vector>
#include "atom.h"

class CellList {

    public: 

    //x, y, z coordinates
    CellList(const AtomRecord& atoms, float cutoff);

    size_t numNeighbors();

    private: 
    //edits: num_atoms, x, y, z max and min, calls set setNumBoxes
    void setAtoms(const AtomRecord& atoms);
    
    void setNumBoxes();

    //temporary; used to place atoms in cells
    std::vector<size_t> setCellCount(const AtomRecord& atoms);
    
    size_t getCellID(float x, float y, float z) const;

    void setBegin(const std::vector<size_t>& cell_count);

    void build(std::vector<size_t>& cell_count, const AtomRecord& atoms);

    //points to indices within atoms
    std::vector<size_t> flat_list;
    //points to indices within flat_list
    std::vector<size_t> begin;

    size_t num_atoms;

    float x_max, y_max, z_max;
    float x_min, y_min, z_min;

    //number of boxes in each direction
    int Nx, Ny, Nz; 
    //equals cutoff
    float cell_length_x, cell_length_y, cell_length_z;

}