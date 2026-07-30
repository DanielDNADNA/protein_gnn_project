import sys
import os

import spatial_engine
import numpy as np

from pathlib import Path

ROOT = Path(__file__).parent

protein_file = ROOT / "tests" / "data" / "1DSB.cif.gz"


def main():
    parser = spatial_engine.Parse_mmCIF(str(protein_file))

    coords, features, target, edge_index, edge_attr = parser.exportTensors()
    print(coords.shape, coords.dtype, sep = ", ")
    print(features.shape, features.dtype, sep = ", ")
    print(target.shape, target.dtype, sep = ", ")
    print(edge_index.shape, edge_index.dtype, sep = ", ")
    print(edge_attr.shape, edge_attr.dtype, sep = ", ")
    
    print(features[0])

    print(edge_attr[29])

    assert edge_index.max() < coords.shape[0]
    assert edge_index.min() >= 0
    assert np.all(edge_attr[:,4] > 0)

    bond_types = edge_attr[:,:4]
    assert np.all(bond_types.sum(axis=1) == 1)

    edges = edge_index.T

    edge_set = set(map(tuple, edges))

    for u, v in edges:
        assert (v, u) in edge_set

    assert not np.any(edge_index[0] == edge_index[1])

    for i in range(edge_index.shape[1]):
        a = edge_index[0,i]
        b = edge_index[1,i]

        calculated = np.linalg.norm(coords[a]-coords[b])

        assert np.isclose(
            calculated,
            edge_attr[i,4],
            atol=1e-5
        )

    avg_degree = 2*edge_index.shape[1]/coords.shape[0]
    print(avg_degree)

    assert np.all(features[:,0:6].sum(axis=1)==1)

    assert np.all(features[:,6:14].sum(axis=1)==1)

    assert np.all(features[:,15:].sum(axis=1)==1)

    assert np.all(edge_attr[:,4] < 2.0)

if __name__ == "__main__":
    main()