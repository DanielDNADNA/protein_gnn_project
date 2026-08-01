#pragma once

#include <map>
#include <vector>
#include <string>
#include <cstdint>

#include "algorithm.h"

struct AlgorithmDetails {

    double build_ms = 0.0;
    double search_ms = 0.0;

    std::vector<uint16_t> neighbors;

    std::string algorithm_name;

    size_t num_atoms = 0;

    std::string protein_name;

};