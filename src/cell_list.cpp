#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>
#include <array>
#include <iostream>
#include <cassert>

#include "cell_list.h"
#include "atom.h"

CellList::CellList(AtomRecord& atoms, float cutoff): atoms(atoms), cutoff(cutoff) {
    
    cutoff_sq = cutoff*cutoff;
    cell_length_x = cutoff;
    cell_length_y = cutoff;
    cell_length_z = cutoff;
    
    setAtoms();

    std::vector<uint32_t> flat_list(num_atoms);
    degree.resize(num_atoms);
    
    setNumBoxes();
    std::vector<uint32_t> cell_count = setCellCount();
    setBegin(cell_count);

    build(cell_count, flat_list);

    sort(flat_list);
}

void CellList::setAtoms() {

    num_atoms = atoms.x.size();

    x_max = atoms.x[0];
    x_min = atoms.x[0];
    y_max = atoms.y[0];
    y_min = atoms.y[0]; 
    z_max = atoms.z[0];
    z_min = atoms.z[0];

    for (size_t i = 1; i < atoms.x.size(); i++) {

        if (x_max < atoms.x[i]) x_max = atoms.x[i];
        if (x_min > atoms.x[i]) x_min = atoms.x[i];

        if (y_max < atoms.y[i]) y_max = atoms.y[i];
        if (y_min > atoms.y[i]) y_min = atoms.y[i];

        if (z_max < atoms.z[i]) z_max = atoms.z[i];
        if (z_min > atoms.z[i]) z_min = atoms.z[i];
    }
}

void CellList::setNumBoxes() {

    Nx = static_cast<size_t>(std::ceil((x_max-x_min)/cell_length_x));
    Ny = static_cast<size_t>(std::ceil((y_max-y_min)/cell_length_y));
    Nz = static_cast<size_t>(std::ceil((z_max-z_min)/cell_length_z));

    Nx = std::max(Nx, static_cast<uint32_t> (1));
    Ny = std::max(Ny, static_cast<uint32_t> (1));
    Nz = std::max(Nz, static_cast<uint32_t> (1));
}

std::vector<uint32_t> CellList::setCellCount() {

    std::vector<uint32_t> cell_count(Nx * Ny * Nz);

    for (size_t i = 0; i < atoms.x.size(); i++) {
        float x = atoms.x[i];
        float y = atoms.y[i];
        float z = atoms.z[i];

        cell_count[getCellID(x, y, z)]++;
    }

    return cell_count;
}

size_t CellList::getCellID(float x, float y, float z) const {

    size_t cell_x = static_cast<size_t>((x - x_min) / cell_length_x);
    size_t cell_y = static_cast<size_t>((y - y_min) / cell_length_y);
    size_t cell_z = static_cast<size_t>((z - z_min) / cell_length_z);

    cell_x = std::min(cell_x, static_cast<size_t>(Nx - 1));
    cell_y = std::min(cell_y, static_cast<size_t>(Ny - 1));
    cell_z = std::min(cell_z, static_cast<size_t>(Nz - 1));

    return cell_x + Nx * cell_y + Nx * Ny * cell_z;
}

void CellList::setBegin(const std::vector<uint32_t>& cell_count) {

    begin.resize(cell_count.size()+1);

    begin[0] = 0;

    for (size_t i = 0; i < cell_count.size(); i++) {
        begin[i+1] = cell_count[i] + begin[i];
    }
}

void CellList::build(std::vector<uint32_t>& cell_count, std::vector<uint32_t>& flat_list) {

    for (size_t i = 0; i < atoms.x.size(); i++) {
        size_t index = getCellID (atoms.x[i], atoms.y[i], atoms.z[i]);
        flat_list[begin[index] + --cell_count[index]] = i;
    }
}

//moves the atoms into the order given by flat_list; fewer cache misses during neighborhood checking
//for this application, neighborhood sort dominates num computations
void CellList::sort(std::vector<uint32_t>& flat_list) {
    AtomRecord sortedAtom;
    sortedAtom.resize(atoms.x.size());

    for (size_t i = 0; i < flat_list.size(); i++) {
        size_t index = flat_list[i];
        sortedAtom.x[i] = atoms.x[index];
        sortedAtom.y[i] = atoms.y[index];
        sortedAtom.z[i] = atoms.z[index];
        sortedAtom.atom_type[i] = atoms.atom_type[index];
        sortedAtom.locant[i] = atoms.locant[index];
        sortedAtom.residue_type[i] = atoms.residue_type[index];
        sortedAtom.atom_index[i] = atoms.atom_index[index];
        sortedAtom.b_factor[i] = atoms.b_factor[index];
    }

    atoms = std::move(sortedAtom);
}

//derive the x, y, z cell coords and then add 1/-1 to each to find the 
//adjacent and diagonal rows
uint16_t CellList::checkNeighboringCells (size_t cell_index, float x, float y, float z) {
    
    int cell_x = cell_index % Nx;
    int cell_y = (cell_index / Nx) % Ny;
    int cell_z = cell_index / (Nx * Ny);
    
    uint16_t num_neighbors = 0; 

    for (int dz = -1; dz <= 1; dz++) 
    for (int dy = -1; dy <= 1; dy++)
    for (int dx = -1; dx <= 1; dx++) {
        int nx = cell_x + dx;
        int ny = cell_y + dy; 
        int nz = cell_z + dz; 

        //filter invalid neighbors (boundary cases)
        if (nx < 0 || nx >= Nx) continue;
        if (ny < 0 || ny >= Ny) continue;
        if (nz < 0 || nz >= Nz) continue;

        size_t neighbor = nx + ny * Nx + Nx * Ny * nz;
        checkCell(neighbor, x, y, z, num_neighbors);
    }

    return num_neighbors;
 }

//note: for the purposes of this implemetation, self is allowed as a neighbor
void CellList::checkCell(size_t cell_index, float x, float y, float z, uint16_t& num_neighbors) {
    
    size_t start_atom = begin[cell_index];
    size_t end_atom = begin[cell_index+1];

    for (size_t i = start_atom; i < end_atom; i++) {
        float dx = atoms.x[i] - x;
        float dy = atoms.y[i] - y;
        float dz = atoms.z[i] - z;
        
        if (dx*dx + dy*dy + dz*dz <= cutoff_sq) {
            num_neighbors++;
        }
    }
}

std::vector<uint16_t> CellList::numNeighbors() {
    
    for (size_t cell_index = 0; cell_index < begin.size() - 1; cell_index++) {
        
        for (size_t atom = begin[cell_index]; atom < begin[cell_index+1]; atom++) {
            
            uint16_t num_neighbors = 0;
            float x = atoms.x[atom];
            float y = atoms.y[atom];
            float z = atoms.z[atom];

            num_neighbors = checkNeighboringCells(cell_index, x, y, z);
            degree[atom] = num_neighbors;
        }
    }
    
    return degree;
}

std::vector<uint16_t> CellList::degreesInOriginalOrder(const std::vector<uint16_t>& degree) 
    const {

    std::vector<uint16_t> reorderedNeighbors(num_atoms);

    for (size_t i = 0; i < atoms.atom_index.size(); i++) {
        uint32_t index = atoms.atom_index[i];
        reorderedNeighbors[index] = degree[i];
    }

    return reorderedNeighbors;
}