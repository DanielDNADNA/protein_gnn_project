import sys
import os

import spatial_engine
import numpy as np

from pathlib import Path
from ml.dataset_statistics import Statistics

from torch_geometric.data import Data
import torch

ROOT = Path(__file__).parent
DATA_DIR = ROOT / "tests" / "training"

PROCESSED_DIR = ROOT / "processed"
PROCESSED_DIR.mkdir(exist_ok=True)

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

    for file in protein_files:
        parser = spatial_engine.Parse_mmCIF(str(file))
        coords, features, target, edge_index, edge_attr = parser.exportTensors()

        stat_features.Standardize(features[:,14]) 
        stat_target.Standardize(target)
        stat_edge_attr.Standardize(edge_attr[:,4])

        x_stats.Standardize(coords[:,0])
        y_stats.Standardize(coords[:,1])
        z_stats.Standardize(coords[:,2])

        data = Data(x=features, edge_index=edge_index, edge_attr=edge_attr, y=target, pos=coords)

        output = PROCESSED_DIR / (file.stem.replace(".cif", "") + ".pt")

        torch.save(data, output)

if __name__ == "__main__":
    main()