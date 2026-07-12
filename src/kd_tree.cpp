
#include <vector>
#include <numeric>
#include <algorithm>

#include "kd_tree.h"
#include "atom.h"

    KDTree::KDTree(const AtomRecord& atoms)
    : atoms(atoms),
      indices(atoms.x.size())
    {
        left.reserve(atoms.x.size());
        right.reserve(atoms.x.size());
        atom_index.reserve(atoms.x.size());
        std::iota(indices.begin(), indices.end(), 0);
    }

    bool KDTree::less (int lhs, int rhs, int axis) const {
        
        switch (axis) {
            case 0: 
                if (atoms.x[lhs] != atoms.x[rhs]) {
                    return atoms.x[lhs] < atoms.x[rhs];
                }
            case 1: 
                if (atoms.y[lhs] != atoms.y[rhs]) {
                    return atoms.y[lhs] < atoms.y[rhs];
                }
            default: 
                if (atoms.z[lhs] != atoms.z[rhs]) {
                    return atoms.z[lhs] < atoms.z[rhs];
                }
        }
        //tie breaks
        if (atoms.x[lhs] != atoms.x[rhs]) {
                return atoms.x[lhs] < atoms.x[rhs];
            }
        else if (atoms.y[lhs] != atoms.y[rhs]) {
                return atoms.y[lhs] < atoms.y[rhs];
            }
        else if (atoms.z[lhs] != atoms.z[rhs]) {
                return atoms.z[lhs] < atoms.z[rhs];
            }
        //should never trigger as one of x, y, or z coords is guaranteed to be different
        else return rhs > lhs;
    }

    size_t KDTree::constructTree(int begin, int end, int depth) {

        if (begin >= end) return -1;

        //axis will cycle between 0-2; x, y, z
        int axis = depth % 3;

        auto median_it = indices.begin() + (begin + end)/2;

        //note: in run 1, grabs the median, but for future runs, guarantees >= and < property
        //for successive axes.
        std::nth_element(
                indices.begin() + begin,
                median_it,
                indices.begin() + end,
                [&](int a, int b)
                {
                    return less(a, b, axis);
                }
            );

        //necessary, as the bounds of the recursion is the median within the index vector; no relation 
        //to the atoms
        size_t median_position = median_it-indices.begin();
        //returns back the equivalent index; then searches the indices for the median value
        //node: median value = index of indices that corresponds to the median of atom.x, y, or z
        //depending on axis
        int median_index = indices[median_position];
        
        size_t node = atom_index.size(); 

        atom_index.push_back(median_index);
        left.push_back(-1);
        right.push_back(-1);

        //correctly assigns the right place for the nodes; .push_back() would fail due to recursion
        left[node] = constructTree(begin, median_position, depth+1);
        right[node] = constructTree(median_position+1, end, depth+1);

        return node;
    }
    
    void KDTree::radiusSearch(float qx, float qy, float qz, float radius, int depth, 
        std::vector<int>& neighbors, int index) {
        
        if (index >= atom_index.size() || index == -1) return;

        //points to a location within atoms
        int node = atom_index[index];

        if (node >= atoms.x.size()) return;

        int axis = depth % 3;

        //check if within range of the node in question
        float dx = atoms.x[node]-qx;
        float dy = atoms.y[node]-qy;
        float dz = atoms.z[node]-qz;

        int rightChild = right[index];
        int leftChild = left[index];

        //if a point is in range, it is a neighbor
        if (dx * dx + dy * dy + dz * dz <= radius * radius) {
            neighbors.push_back(node); 
        }

        if (right[index] == -1 || left[index] == -1) return;

        

        return;
    }

    int nodeDegree () {

    }


