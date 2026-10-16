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

### KD-tree and Cell List 

The number of neighbors of each atom can be found efficiently by the kd tree and the cell lists, which 
algorithmically prune areas that are outside of the radius. 

#### KD-tree

The kd-tree is a special case of binary tree, which instead divides space based on the axes. For example, if the axis is 1, 
which corresponds to the x-axis, all of the right children will have a greater or equal x coordinate than all of the left children. The next recursion does the same for the y-axis, and the next for the z-axis, before looping back to x. 

This implementation uses an implicit kd-tree; where the indices of the atoms are stored in a separate vector, and the right and left vectors store the positions of the next right and left child respectively. 

The search method uses the size guarantees to prune, on average, half of all atoms each recursion. For example, if a point is (1,0,0), the radius is 5.0 Angstroms, and the current point is (7,0,0), the distance would be calculated 1-7 = -6. The x value would be too large, so the algorithm would check only the left children. Since on average the tree is roughly balanced, that eliminates about half of all atoms. 

Because of this, the search algorithm has a $O(\log n)$ time complexity per neighbor. 

#### Cell list

The cell list uses spatial hashing to narrow down the atoms to a smaller box, and then uses the brute-force method to search for neighbors within that bounded range. Each atom is assigned to a box of size 5.0f via an ID. 

The search algorithm converts the atom's coordinates into its ID, and then checks its box and the 9 surrounding ones for points within the radius of the atom. 

While the theoretical time complexity for a brute force approach is $O(N^2)$, the number of elements within a box is relatively fixed for an object of constant density, which proteins roughly are, yielding an average search time of $O(1)$ per atom. 

This cell list was implemented using a CSR format, as well as an implicit cell list structure, by ordering the AtomRecord struct such that atoms in the same box are adjacent to each other. 

In addition, a basic OMP cell list was implemented, by making the search within the 9 boxes parallel. 

## Benchmarking

While the cell list has a superior theoretical time complexity than the kd-tree, $\log n$ is often nearly constant over 
a few orders of magnitude, thus constant factors may enable it to be a faster algorithm. In general, the cell list does better in structures of roughly regular density, while doing worse in ones without (as some areas may have an extremely 
large amount of atoms to check). By contrast, the kd-tree does best in structures with roughly the same number of atoms on each 'split' of the axes, as otherwise less than half of the points would be pruned each recursion, so is not dependent on density. 

The build and search steps for the kd-tree and cell list were timed using the chrono library from STL.
For each trial, a choice was made between using the cell list, an OMP cell list, and a kd-tree first. This prevented any one data structure from having the advantage of going first, although this did introduce some variability. 

The data was then fit to a quadratic, linear, or n log n model for that algorithm, and the $r^2$ and rmse of each were reported. 

### Specs
- **Hardware:** Intel(R) Core(TM) Ultra 7 265 
- **Compiler:** gcc 15.2.0, no flags

### Results (in miliseconds)
#### Brute force:
$f(n) = 4.49 \cdot 10^{-7}n^2 + 7.18+10^{-5}n + 0.256$, $r^2$ = 0.9997, rmse = 0.92

The KD-tree and the cell list's times were divided between search and build. 

| Left Aligned | Centered | Right Aligned |
| :---         | :---:    | ---:          |
|          | Text     | Text          |
#### KD-tree: 

Search:
Unexpectedly, the linear case performed best, with an $r^2$ score of 0.996 and an rmse of 29.00, 
$f(n) = 9.56 \cdot 10^{-4}n -0.77$

While the n log n form, the theoretical time complexity, had an $r^2$ score of 0.995 and an rmse of 37, along with a very large constant term and a small n ln n term, which would normally suggest that it was not very relevant. 
$f(n) = 6.64 \cdot 10^{-5}n \log n + 11.054$

#### Cell List:

Search: 

