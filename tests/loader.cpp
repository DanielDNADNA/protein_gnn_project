#include "protein_gnn_project.h"
#include "atom.h"
#include "evaluation.h"

#include <vector>
#include <iostream>
#include <filesystem>
#include <cmath>
#include <cassert>

namespace fs = std::filesystem;

int main() {

    for (const auto& bucket : fs::directory_iterator("benchmark_bin")) {
        if (!bucket.is_directory())
            continue;
        std::cout << "===== " << bucket.path().filename() << " =====\n";

        for (const auto& file : fs::directory_iterator(bucket.path())) {
            if (file.path().extension() != ".bin")
                continue;

            AtomRecord atoms;
            atoms = atoms.load(file.path().string());

            Evaluation evaluate(atoms, bucket.path().filename());
            evaluate.buildCellList();
            std::vector<uint16_t> cell_neighbors = evaluate.getCellNeighbors();
            evaluate.buildTree();
            std::vector<uint16_t> tree_neighbors = evaluate.getTreeNeighbors();
            std::vector<uint16_t> brute_force_neighbors = evaluate.getBruteForceNeighbors();

            std::cout << "Evaluating tree_neighbors using brute force: \n";
            assert(evaluate.validate(tree_neighbors, brute_force_neighbors));
            std::cout << "KD tree validated \n";
            std::cout << "Evaluating cell_neighbors using brute force: \n";
            assert(evaluate.validate(brute_force_neighbors, cell_neighbors));
            std::cout << "Cell list validated \n";

            evaluate.printCellTime();
            evaluate.printTreeTime();
            evaluate.printBruteForceTime();
        }
    }

    return 0;
}
