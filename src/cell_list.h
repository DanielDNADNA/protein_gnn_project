#pragma once

#include <vector>
#include <array>

#include "atom.h"

class CellList {

    public: 

    //x, y, z coordinates
    CellList(AtomRecord& atoms, float cutoff);

    std::vector<uint16_t> numNeighbors();

    private: 
    //edits: num_atoms, x, y, z max and min, calls set setNumBoxes
    void setAtoms();
    
    void setNumBoxes();

    //temporary; used to place atoms in cells
    std::vector<uint32_t> setCellCount();
    
    size_t getCellID(float x, float y, float z) const;

    void setBegin(const std::vector<uint32_t>& cell_count);

    void build(std::vector<uint32_t>& cell_count, std::vector<uint32_t>& flat_list);

    void sort(std::vector<uint32_t>& flat_list);

    //Fills array with IDs of the 26 surrounding cell IDs
    void formOffsets ();
    
    //Adds neighbors to neighbors vector
    void checkCell(size_t cell_index, float x, float y, float z, uint16_t& num_neighbors);

    //will become implicit cell list
    AtomRecord& atoms;
    //points to indices within sorted atom_record
    std::vector<uint32_t> begin;

    uint32_t num_atoms;

    float x_max, y_max, z_max;
    float x_min, y_min, z_min;

    //number of boxes in each direction
    uint32_t Nx, Ny, Nz; 
    //equals cutoff
    float cell_length_x, cell_length_y, cell_length_z;

    float cutoff;

    float cutoff_sq;

    //gives indices of the surrounding 26 cells when added to the cell_index
    std::array<int, 27> offsets;

    std::vector<uint16_t> degree;
};