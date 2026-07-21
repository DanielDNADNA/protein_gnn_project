import sys
import os

import spatial_engine
import numpy as np

from pathlib import Path
from dataset_statistics import Statistics

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "tests" / "training"

protein_files = sorted(DATA_DIR.glob("*.cif.gz*"))

def main():

    stat_features = Statistics() 
    stat_target = Statistics()
    stat_edge_attr = Statistics()

    x_stats = Statistics()
    y_stats = Statistics()
    z_stats = Statistics()

    for file in protein_files:
        parser = spatial_engine.Parse_mmCIF(str(file))
        coords, features, target, edge_index, edge_attr = parser.exportTensors()

        x_stats.WelfordAlgorithm(coords[:,0]),
        y_stats.WelfordAlgorithm(coords[:,1]),
        z_stats.WelfordAlgorithm(coords[:,2])

        stat_features.WelfordAlgorithm(features[:,14])

        stat_target.WelfordAlgorithm(target)
        
        stat_edge_attr.WelfordAlgorithm(edge_attr[:,4])

        print ("Preprocessing step")

        print ("First row of coords: ", coords[0])

        print ("First row of features: ", features[0])

        print ("First row of target: ", target[0])
        print ("Second row of target:", target[1])

        print ("First row of edge_attr: ", edge_attr[0])

    for file in protein_files:
        parser = spatial_engine.Parse_mmCIF(str(file))
        coords, features, target, edge_index, edge_attr = parser.exportTensors()

        stat_features.Standardize(features[:,14]) 
        stat_target.Standardize(target)
        stat_edge_attr.Standardize(edge_attr[:,4])

        x_stats.Standardize(coords[:,0])
        y_stats.Standardize(coords[:,1])
        z_stats.Standardize(coords[:,2])

        features_std = stat_features.std
        target_std = stat_target.std
        edge_attr_std = stat_edge_attr.std

        print ("Processing statistics: ")
        print ("x_mean: ", x_stats.mean)
        print ("x_std: ", x_stats.std)
        print ("y_mean: ", y_stats.mean)
        print ("y_std: ", y_stats.std)
        print ("z_mean: ", z_stats.mean)
        print ("z_std: ", z_stats.std)

        print ("mean of target: ", stat_target.mean)
        print ("std of target: ", target_std)
        print ("mean of edge_attr", stat_edge_attr.mean)
        print ("std of edge_attr", edge_attr_std)
        print ("mean of features", stat_features.mean)
        print ("std of features", features_std)

        print ("Post-processing:" )
        print ("First row of coords: ", coords[0])
        print(coords.shape)

        print ("First row of features: ", features[0])
        print(features.shape)

        print ("First row of target: ", target[0])
        print ("Second row of target: ", target[1])
        print(target.shape)

        print ("First row of edge_attr: ", edge_attr[0])
        print (edge_attr.shape)

if __name__ == "__main__":
    main()