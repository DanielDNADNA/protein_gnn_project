#include <vector>
#include <iostream>
#include <filesystem>

#include "atom.h"
#include "protein_gnn_project.h"

namespace fs = std::filesystem;

int main() {

    for (const auto& entry : fs::directory_iterator(BENCHMARK_RAW_DIR)) {

        if (entry.path().extension() != ".cif")
            continue;

        auto pdb_id = entry.path().stem().string();

        std::cout << "Processing: " << pdb_id << "\n";

        Parse_mmCIF parse(entry.path().string());

        const AtomRecord& atoms = parse.atoms;

        size_t num_atoms = atoms.x.size();

        std::string folder;

        if (num_atoms < 1000) {
            folder = "1k";
        }
        
        else if (num_atoms >= 1000 && num_atoms < 5000) {
            folder = "5k";
        }

        else if (num_atoms >= 5000 && num_atoms < 10000) {
            folder = "10k";
        }

        else if (num_atoms >= 10000 && num_atoms < 25000) {
            folder = "25k";
        }

        else if (num_atoms >= 25000 && num_atoms < 50000) {
            folder = "50k";
        }

        else if (num_atoms >= 50000 && num_atoms < 100000) {
            folder = "100k";
        }
        
        else if (num_atoms >= 100000 && num_atoms < 250000) {
            folder = "250k";
        }

        else if (num_atoms >= 250000 && num_atoms < 500000) {
            folder = "500k";
        }

        else if (num_atoms >= 500000 && num_atoms < 1000000) {
            folder = "1m";
        }

        else if (num_atoms >= 1000000 && num_atoms < 2000000) {
            folder = "2m";
        }

        else if (num_atoms >= 2000000 && num_atoms < 3000000) {
            folder = "3m";
        }

        fs::path bin_dir = fs::path(BENCHMARK_BIN_DIR) / folder;

        fs::create_directories(bin_dir);

        fs::path output_file = bin_dir / (pdb_id + ".bin");

        atoms.save(output_file);
    }

    return 0;
}