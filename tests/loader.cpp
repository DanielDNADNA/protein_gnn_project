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

    //Add algorithms here to validate/evaluate them: 
    //BruteForce, KDTree, CellList, CellListOMP, CudaCellList
    auto algorithms = {
            Algorithm::CellList, 
            Algorithm::OMPCellList
        };

    //stops validating/evaluating if max_num_atoms passes this mark
    size_t max_num_atoms = 10000000;

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

            size_t num_atoms = atoms.x.size();
            
            if (num_atoms <= max_num_atoms) {

                Evaluation evaluate(atoms, file.path().stem().string());
                evaluate.benchmarkAlgorithms(algorithms);

                //validation uses the first algorithm as the base case
                if (algorithms.size() >= 2) {
                    evaluate.validateAlgorithms(algorithms);
                }
            }
            else break;
        }
    }

    return 0;
}
