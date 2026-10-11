# Extraction of Molecular Features and Spatial Data Structures Benchmarking

## Overview 
On the C++ side, this project extracts molecular features from mmCIF files and then calculates the distances between 
each molecule within a 5.0 Angstrom radius. The molecular features are packaged for use by a Pytorch GNN (unfinished), 
while the distance calculations were benchmarked for manually-implemented bruteforce, KD-tree and cell list approaches.

## Features
- **Selection and Extraction:**
The RCSB Protein Data Bank API was used to select and download mmCIF files of molecules of varying sizes.
- **Structure Parsing:**
Gemmi was then used to extract atomic coordinates, residue IDs, and other relevant features. These were then organized into tensors.
- **Neighborhood search:**
The number of atoms within 5.0 Angstroms of each atom was calculated using a cell list, which was benchmarked against
bruteforce and the KD-tree.
- **Benchmarking:**
The cell list and KD-tree were validated using brute force for molecules under 25,000 atoms. Beyond that, due to scaling,
the cell list and KD-tree were validated against each other. The time required to build the data structures
and make the search were recorded and compared as well.
- **GNN:**
Battaglia-style Meta GNN made in accordance to Sarparast et al's 2024 Arxiv preprint paper. It uses molecular features to
predict the B-factor of each atom in a protein, or the ability of the atom to deviate from its initial position. 
(Unfinished)

## Implementation

### Parsing and Feature Extraction
This project used Gemmi, a C++ based parser of the PDB/mmCIF files that are stored on the PDB databank. It provides superior 
performance to Biopython, another popular library, as well as being in the author's preferred language.

The features and the GNN architecture were exactly those used in a 2024 Arxiv preprint by Sarparast et al., as this project
began as a replication of the Meta GNN used within the paper. 

As explained within the paper, there are four matrices representing each matrix (section 5.2): 
node features, 3D coordinates, edge features, target node feature. The number of node and edge features was not explicitly given
and were assumed by the author. 

#### The matrix of node features. 

A [num atoms, 36]  tensor with the following features: 
- **Atom type:** Indices 0-5: Carbon, nitrogen, oxygen, sulfur, phosphorous, and other encoded via one hot
- **Atom locant:** Indices 6-13: Alpha, beta, gamma, ..., eta, and other assigned respectively.
- **Atom degree:** Index 14: Number of neighbors, normalized. 
- **Residue type:** Indices 15-36: the residue the atom belongs to (20 choices) and another for none/unknown.

#### The matrix of 3D coordinates.

A [num atoms, 3] tensor with x, y, z coordinates of that atom respectively, normalized. 

#### The matrix of edge features

A [num edges, 4] tensor with the following features: 
- **Bond type:** Indices 0-3: Single, double, triple, and delocalized respectively; default of single.
- **Bond distance:** Index 4: The distance between two bonded atoms, normalized.


