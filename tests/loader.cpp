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
    
    auto algorithms = {
                Algorithm::CellList, 
                Algorithm::CellListOMP};

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
            
            evaluate.evaluateAlgorithms(algorithms);
            evaluate.validateAlgorithms(algorithms);
        }
    }

    return 0;
}
