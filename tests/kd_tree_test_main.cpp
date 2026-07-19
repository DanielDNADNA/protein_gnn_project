#include "kd_tree.h"
#include "atom.h"
#include <chrono>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <random>

void bruteForceSearch(float qx, float qy, float qz, float radius, 
    const AtomRecord& atoms, std::vector<int>& neighbors);
void parseFile (const std::string& fileName, AtomRecord& coords1, float& radius);
void generateRandomAtoms(AtomRecord& atoms, size_t numAtoms);

int main() {
    AtomRecord coords1;
    int numAtoms = 20000;
    generateRandomAtoms(coords1, numAtoms);

    std::cout << "For " << numAtoms << " atoms: \n";

    float radius = 5.0f;

    auto start = std::chrono::high_resolution_clock::now();

    KDTree kdTree {coords1};

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "KD-tree construction: " << elapsed.count() << " ms\n";

    auto count = elapsed.count();

    std::mt19937 rng(123);

    std::vector<int> neighbors1;
    std::vector<int> neighbors2;
    neighbors1.reserve(numAtoms/10);
    neighbors2.reserve(numAtoms/10);

    float totalNeighbors = 0;

    /*start = std::chrono::high_resolution_clock::now();
    for (int q = 0; q < coords1.x.size(); q++) {

        float qx = coords1.x[q];
        float qy = coords1.y[q];
        float qz = coords1.z[q];

        neighbors1.clear();

        bruteForceSearch(
            qx, qy, qz,
            radius,
            coords1,
            neighbors1
        );
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Brute force search: " << elapsed.count() << " ms\n";

    */
    start = std::chrono::high_resolution_clock::now();
    for (int q = 0; q < coords1.x.size(); q++) {

        float qx = coords1.x[q];
        float qy = coords1.y[q];
        float qz = coords1.z[q];

        neighbors2.clear();

        kdTree.radiusSearch(
            qx, qy, qz,
            radius,
            0,
            neighbors2,
            0
        );

        totalNeighbors += neighbors2.size();
    }

    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "KD-tree search: " << elapsed.count() << " ms\n";
    std::cout << "Total KD-tree time: " << count + elapsed.count() << " ms \n";

    std::cout << "Average degree: " << static_cast<double>(totalNeighbors) / coords1.x.size() << '\n';

    /*std::sort(neighbors1.begin(), neighbors1.end());
    std::sort(neighbors2.begin(), neighbors2.end());

    if (neighbors1 == neighbors2) {
        std::cout << "Validated \n";
    }

    else std::cout << "Error \n";
    */
    return 0;
}

void bruteForceSearch(float qx, float qy, float qz, float radius, 
    const AtomRecord& atoms, std::vector<int>& neighbors) {
    
    neighbors.clear();
    
    for (size_t j = 0; j < atoms.x.size(); j++) {
        float dx = qx - atoms.x[j];
        float dy = qy - atoms.y[j];
        float dz = qz - atoms.z[j];

        if (dx * dx + dy * dy + dz * dz <= radius * radius) {
            neighbors.push_back(j);
        }
    }

    return;
}

void parseFile (const std::string& fileName, AtomRecord& coords1, float& radius) {
    
    std::fstream inputFile (fileName);

    if (!inputFile) {
    throw std::runtime_error("Could not open file: " + fileName);
    }

    int numAtoms = 0; 
    inputFile >> numAtoms >> radius; 

    coords1.reserve(numAtoms);

    float x, y, z;

    while (inputFile >> x >> y >> z) {
        coords1.x.push_back(x);
        coords1.y.push_back(y);
        coords1.z.push_back(z);
    }
    
    return;
}

void generateRandomAtoms(AtomRecord& atoms, size_t numAtoms) {
    std::mt19937 rng(42); // fixed seed = reproducible

    std::normal_distribution<float> dist(0.0f, 35.0f);

    atoms.reserve(numAtoms);

    for (size_t i = 0; i < numAtoms; i++) {
        atoms.x.push_back(dist(rng));
        atoms.y.push_back(dist(rng));
        atoms.z.push_back(dist(rng));
    }
}

/*
100,000 atoms
---------------
KD-tree build      : 72 ms
KD-tree queries    : 482 ms
Total              : 554 ms

Brute force        : 51,355 ms

Average degree     : 28.3
Speedup (queries)  : 106.5×
Speedup (overall)  : 92.7×
*/
/*
int main() {

    std::string fileName = "sample.txt";
    AtomRecord coords1;
    float radius;

    parseFile (fileName, coords1, radius);
    std::vector<int> neighbors1;
    
    std::vector<int> neighbors2;
    KDTree kdTree {coords1};

for (size_t i = 0; i < coords1.x.size(); i++) {
    float qx = coords1.x[i];
    float qy = coords1.y[i];
    float qz = coords1.z[i];

    neighbors2.clear();

    bruteForceSearch(qx, qy, qz, radius, coords1, neighbors1);
    kdTree.radiusSearch(qx, qy, qz, radius, 0, neighbors2, 0);

    std::sort(neighbors1.begin(), neighbors1.end());
    std::sort(neighbors2.begin(), neighbors2.end());

    if (neighbors1.size() != neighbors2.size()) {
        std::cout << "Mismatch size" << '\n';
        std::cout << "Brute force: " << neighbors1.size() << '\n';
        std::cout << "KD tree: " << neighbors2.size() << '\n';
        return 1;
    }
    
    for (size_t j = 0; j < neighbors1.size(); j++) {
        if (neighbors1[j] != neighbors2[j]) {
            std::cout << "Mismatch at index " << j << '\n';
            return 1; 
        }
    }
}

    std::cout << "Validation: KD-tree matches brute force" << '\n';
    return 0;
}*/