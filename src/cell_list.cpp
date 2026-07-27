#include <vector>
#include <algorithm>
#include <cmath>

#include "cell_list.h"
#include "atom.h"


CellList::CellList(const AtomRecord& atoms, float cutoff) {

    cell_length_x = cutoff;
    cell_length_y = cutoff;
    cell_length_z = cutoff;
    setAtoms(atoms);
    setNumBoxes();
    std::vector<size_t> cell_count = setCellCount(atoms);
    setBegin(cell_count);

    build(cell_count, atoms);

}

void CellList::setAtoms(const AtomRecord& atoms) {

    num_atoms = atoms.x.size();

    x_max = atoms.x[0];
    x_min = atoms.x[0];
    y_max = atoms.y[0];
    y_min = atoms.y[0]; 
    z_max = atoms.z[0];
    z_min = atoms.z[0];;

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

    Nx = std::max(Nx, 1);
    Ny = std::max(Ny, 1);
    Nz = std::max(Nz, 1);
}

std::vector<size_t> CellList::setCellCount(const AtomRecord& atoms) {

    std::vector<size_t> cell_count(Nx * Ny * Nz);

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

void CellList::setBegin(const std::vector<size_t>& cell_count) {

    begin.reserve(cell_count.size()+1);

    begin.push_back(0);

    for (size_t i = 0; i < cell_count.size(); i++) {
        begin.push_back(cell_count[i] + begin[i]);
    }
}

void CellList::build(std::vector<size_t>& cell_count, const AtomRecord& atoms) {
    
    flat_list.resize(num_atoms);

    for (size_t i = 0; i < atoms.x.size(); i++) {
        size_t index = getCellID (atoms.x[i], atoms.y[i], atoms.z[i]);
        flat_list[begin[index] + --cell_count[index]] = i;
    }
}


