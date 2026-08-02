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
#include <random>

namespace fs = std::filesystem;

int main() {

    //Add algorithms here to validate/evaluate them: 
    //BruteForce, KDTree, CellList, CellListOMP, CudaCellList
    std::vector<Algorithm> algorithms = {
            Algorithm::CellList, 
            Algorithm::OMPCellList,
            Algorithm::KDTree    
        };
    
    //to randomly shuffle algorithms
    std::mt19937 rng(42);

    //stops validating/evaluating if max_atoms passes this mark
    size_t max_atoms = 10000000;
    size_t min_atoms = 25000;

    std::vector<fs::path> buckets;

    for (const auto& entry : fs::directory_iterator("benchmark_bin")) {
        if (entry.is_directory()) {
            buckets.push_back(entry.path());
        }
    }

    //sorts the files least atoms to most
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
            
            std::shuffle(algorithms.begin(), algorithms.end(), rng);

            AtomRecord atoms;
            atoms = atoms.load(file.path().string());

            size_t num_atoms = atoms.x.size();
            
            if (num_atoms <= max_atoms && num_atoms >= min_atoms) {

                Evaluation evaluate(atoms, file.path().stem().string(), algorithms);
                evaluate.reportBenchmarks();

                //validation uses the first algorithm as the base case
                if (algorithms.size() >= 2) {
                    evaluate.validateAlgorithms();
                }

                evaluate.exportJson();
            }
            
            if (num_atoms > max_atoms) break;
            if (num_atoms < min_atoms) continue;
        }
    }

    return 0;
}
