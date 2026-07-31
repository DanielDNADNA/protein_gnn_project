#include "protein_gnn_project.h"
#include "atom.h"
#include "evaluation.h"

#include <vector>
#include <iostream>
#include <string>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <cassert>

namespace fs = std::filesystem;

int main() {

    std::vector<fs::path> buckets;

    for (const auto& entry : fs::directory_iterator("benchmark_bin")) {
        if (entry.is_directory()) {
            buckets.push_back(entry.path());
        }
    }

    //sorts the files least to greatest
    std::sort(buckets.begin(), buckets.end(), 
        [](const fs::path& a, const fs::path& b) {
            auto parse_size = [](const fs::path& p) {
                std::string s = p.filename().string();

                if (s.back() == ('k')) return std::stoi(s.substr(0, s.size()-1)) * 1000;

                if (s.back() == ('m')) return std::stoi(s.substr(0, s.size()-1)) *1000000;

                return 0;
            };   
            return parse_size(a) < parse_size(b);
        });

    for (const auto& bucket : buckets) {
        std::cout << "===== " << bucket.filename().filename() << " =====\n";

        for (const auto& file : fs::directory_iterator(bucket)) {
            if (file.path().extension() != ".bin")
                continue;

            AtomRecord atoms;
            atoms = atoms.load(file.path().string());

            Evaluation evaluate(atoms, file.path().stem().string());
            
            evaluate.buildCellList();
            std::vector<uint16_t> cell_neighbors = evaluate.getCellNeighbors();
            
            evaluate.buildTree();
            std::vector<uint16_t> tree_neighbors = evaluate.getTreeNeighbors();
            
            if (atoms.x.size() <= 5000) {
                std::vector<uint16_t> brute_force_neighbors = evaluate.getBruteForceNeighbors();

                std::cout << "Evaluating tree_neighbors using brute force: \n";

                if (!evaluate.validate(tree_neighbors, brute_force_neighbors)) {
                    throw std::runtime_error("KD tree and brute force mismatch");
                }
        
                std::cout << "KD tree validated \n";

                std::cout << "Evaluating cell_neighbors using brute force: \n";
            
                if (!evaluate.validate(brute_force_neighbors, cell_neighbors)) {
                    throw std::runtime_error("Brute force and cell list mismatch");
                }
                
                std::cout << "Cell list validated \n";

                std::cout << "Evaluating cell_neighbors with kd tree: \n";
                if (!evaluate.validate(tree_neighbors, cell_neighbors)) {
                    throw std::runtime_error("KD tree and cell list mismatch");
                }

                evaluate.printBruteForceTime();
            }

            else {
                std::cout << "Comparing cell list to kd tree \n";

                if (!evaluate.validate(tree_neighbors, cell_neighbors)) {
                    throw std::runtime_error("KD tree and cell list mismatch");
                }

                std::cout << "Cell list matches kd tree \n";
            }

            evaluate.printCellTime();
            evaluate.printTreeTime();
            std::cout << "\n";
        }
    }

    return 0;
}
