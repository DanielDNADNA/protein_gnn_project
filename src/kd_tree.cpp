
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

        KDTree::constructTree(0, atoms.x.size(), 0);
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

    int KDTree::constructTree(int begin, int end, int depth) {

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
        
        int node = atom_index.size(); 

        atom_index.push_back(median_index);
        left.push_back(-1);
        right.push_back(-1);

        //correctly assigns the right place for the nodes; .push_back() would fail due to recursion
        left[node] = constructTree(begin, median_position, depth+1);
        right[node] = constructTree(median_position+1, end, depth+1);

        return node;
    }
    
    //given the coords of a point, finds all points within a given radius
    void KDTree::radiusSearch(float qx, float qy, float qz, float radius, int depth, 
        std::vector<int>& neighbors, int index) {
        //node: index corresponds to atom_index; points to a index in atoms object
        //so left[inx] = index -> atom_index[index] = node -> atoms.x[node] (the coord)
        //qx, qy, qz are coords of the point whose neighbors I am trying to find
        
        if (index >= atom_index.size() || index == -1) return;

        //points to a location within atoms
        int node = atom_index[index];

        if (node >= atoms.x.size() || node == -1) return;

        int axis = depth % 3;

        float x = atoms.x[node];
        float y = atoms.y[node];
        float z = atoms.z[node];
        //check if within range of the node in question
        float dx = qx - x;
        float dy = qy - y;
        float dz = qz - z;

        //if a point is in range, it is a neighbor
        if (dx * dx + dy * dy + dz * dz <= radius * radius) {
            neighbors.push_back(node); 
        }

        float distance;

        switch (axis) {
            case 0: distance = dx; 
                    break;
            case 1: distance = dy; 
                    break;
            default: distance = dz; 
                    break;
        }

        int rightChild = right[index];
        int leftChild = left[index];

        //the range of possible points is given by [-rad+qx, rad+qx]
        //if the original point is within that range, check both children

        if (std::abs(distance) <= radius) {
            if (rightChild != -1) {
                radiusSearch(qx, qy, qz, radius, depth+1, neighbors, rightChild);
            }
            if (leftChild != -1) {
                radiusSearch(qx, qy, qz, radius, depth+1, neighbors, leftChild);
            }
        }

        //if distance > 0, I am on the left side of the point, so move right
        else if (distance > 0 && rightChild != -1) {
            radiusSearch(qx, qy, qz, radius, depth+1, neighbors, rightChild);
        }
        
        else if (distance < 0 && leftChild != -1) {
            radiusSearch(qx, qy, qz, radius, depth+1, neighbors, leftChild);
        }

        return;
    }

    //.reserve will be done outside of the function
    size_t KDTree::getNodeDegree (float qx, float qy, float qz, float radius, std::vector<int>& neighbors) {
        neighbors.clear();
        radiusSearch(qx, qy, qz, radius, 0, neighbors, 0);

        //radius search includes self as neighbor; excluded by spatial degree
        if (neighbors.size() > 0) {
            return neighbors.size() - 1 ;
        }
        else return 0;
    }


