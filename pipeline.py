import sys
import os

import spatial_engine
import numpy as np

from pathlib import Path
from dataset_statistics import Statistics
from normalizer import Normalizer

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "tests" / "data"

protein_files = sorted(DATA_DIR.glob("*.cif.gz*"))

def main():

    stat_features = Statistics() 
    stat_target = Statistics()
    stat_edge_attr = Statistics()

    x_stats = Statistics()
    y_stats = Statistics()
    z_stats = Statistics()

    for file in protein_files:
        parser = spatial_engine.ParsemmCIF(str(file))
        coords, features, target, edge_index, edge_attr = parser.exportTensors()

        x_stats.WelfordAlgorithm(coords[:,0]),
        y_stats.WelfordAlgorithm(coords[:,1]),
        z_stats.WelfordAlgorithm(coords[:,2])
        
        coords_mean = np.array([
            x_stats.mean, 
            y_stats.mean, 
            z_stats.mean
        ])

        stat_features.WelfordAlgorithm(features[:,15])

        stat_target.WelfordAlgorithm(target[:,0])
        
        stat_edge_attr.WelfordAlgorithm(edge_attr[:,4])
        

    for file in protein_files:
        parser = spatial_engine.ParsemmCIF(str(file))
        coords, features, target, edge_index, edge_attr = parser.exportTensors()
        tensor_tuple = (coords, features, target, edge_index, edge_attr)

        features_std = stat_features.CalculatePopulationStd()
        target_std = stat_target.CalculatePopulationStd()
        edge_attr_std = stat_edge_attr.CalculatePopulationStd()

        stats_tuple = (coords_mean, stat_target.mean, target_std, stat_features.mean, features_std, stat_edge_attr.mean, edge_attr_std)
        normalizer = Normalizer(tensor_tuple, stats_tuple)

if __name__ == "__main__":
    main()