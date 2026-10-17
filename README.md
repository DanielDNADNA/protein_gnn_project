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
bruteforce and the kd-tree.
- **Benchmarking:**
The cell list and kd-tree were validated using brute force for molecules under 25,000 atoms. Beyond that, due to scaling,
the cell list and kd-tree were validated against each other. The time required to build the data structures
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
node features, 3D coordinates, edge features, target node feature. The number of node and edge features was not explicitly given and were assumed by the author. 

#### The matrix of node features. 

A [num atoms, 36]  tensor with the following features: 
- **Atom type:** Indices 0-5: Carbon, nitrogen, oxygen, sulfur, phosphorous, and other encoded via one hot
- **Atom locant:** Indices 6-13: Alpha, beta, gamma, ..., eta, and other assigned respectively.
- **Atom degree:** Index 14: Number of neighbors, normalized. 
- **Residue type:** Indices 15-36: the residue the atom belongs to (20 choices) and another for none/unknown.

#### The matrix of 3D coordinates.

A [num atoms, 3] tensor with x, y, z coordinates of that atom respectively, normalized. 

#### The matrix of edge features

A [num edges, 5] tensor with the following features: 
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
large amount of atoms to check). By contrast, the kd-tree does best in structures with roughly the same number of atoms on each 'split' of the axes, as otherwise less than half of the points would be pruned each recursion, so is not dependent on density. It is also possible that the kd-tree may be faster over some size ranges than the cell list and vice versa. The benchmarking aimed to resolve these questions for specifically proteins. 

20 atoms were chosen from the range of (1000, 5000), (5000, 10000, (10000, 50000) ... (1 * 10^6, 3 * 10^6) to measure the performance of the algorithms over a wide range of atom numbers. Because of this, however, the data below should not be considered representative of the algorithms' performance with protein data in general. As stated before, the brute force algorithm was tested only up to 25,000 atoms due to scaling concerns. 

The build and search steps for the kd-tree and cell list were timed using the chrono library from STL.
For each trial, a choice was made between using the cell list, an OMP cell list, and a kd-tree first. This prevented any one data structure from having the advantage of going first, although this did introduce some variability. 

The data was then fit to a quadratic, linear, or n log n model for that algorithm, and the $r^2$ and rmse of each were reported. The $r^2$ term measures the amount of time the algorithm takes that is explained by the number of atoms; generally a higher number is best. In general practice, an $r^2$ of positive or negative 0.90 is considered best, although in this case a value closer to 0.99 is desired. The number of atoms is the dominant factor in the runtime, so lower values are much more likely indicative of a poorly fitted model. By contrast, a low rmse is desired, as it measures the mean error the model has with the actual data. There is no universal standard for a good rmse, and it should be used to compare competing models for the same algorithm. 

### Specs
- **Hardware:** Intel(R) Core(TM) Ultra 7 265 
- **Compiler:** gcc 15.2.0, no flags
- **Environment:** Visual Studios, WSL 3.0.1

### Results (in miliseconds)

The kd-tree and the cell list's times were divided between search and build. 

Functions of Best Fit for Brute Force, K-d Tree, and Cell List; time in milliseconds, n in number of atoms, excluding water
| Brute Force | Linear | Quadratic | n Ln n |
| :---:       | :---:  | :---:     | :---:            |
| Search      | $f(n) = 9.36 \cdot 10^{-3}n - 27.98$ | $f(n) = 4.49 \cdot 10^{-7}n^2 + 7.18 \cdot 10^{-5}n + 0.256$ | $f(n) = 9.42 \cdot 10^{-4} n \ln n - 22.95$ |
| rmse       | 18.81    | 0.921     | 16.62   | 
| r^2        |  0.892   | 0.9997    | 0.916   |

| KD-tree | Linear | Quadratic | n Ln n |
| :---:       | :---:  | :---:     | :---:            |
| Search      | $f(n) = 9.5 \cdot 10^{-4}n - 0.768$ | $f(n) = -4.11 \cdot 10^{-11}n^2 + 1.03 \cdot 10^{-3}n - 6.320$ | $f(n) = 6.64 \cdot 10^{-5} n \ln n + 11.05$ |
| rmse       | 28.96    | 25.30     | 37.10   | 
| r^2        |  0.997   | 0.998    | 0.995   |
| Build      | $f(n) = 1.81 \cdot 10^{-4}n - 1.54$ | $f(n) = 4.19 \cdot 10^{-12}n^2 + 1.73 \cdot 10^{-4}n - 0.97$ | $f(n) = 1.26 \cdot 10^{-5} n \ln n + 0.64$ |
| rmse       | 3.67    | 3.37     | 3.40   | 
| r^2        |  0.999   | 0.999    | 0.999  |

| Cell List | Linear | Quadratic | n Ln n |
| :---:       | :---:  | :---:     | :---:            |
| Search      | $f(n) = 3.70 \cdot 10^{-4}n + 0.385$ | $f(n) = -5.39 \cdot 10^{-12}n^2 + 3.79 \cdot 10^{-4}n - 0.344$ | $f(n) = 2.57 \cdot 10^{-5} n \ln n + 4.92$ |
| rmse       | 8.63    | 8.42 | 11.38   | 
| r^2        |  0.998   | 0.998 | 0.997   |
| Build      | $f(n) = 2.72 \cdot 10^{-5}n - 0.32$ | $f(n) = 3.38 \cdot 10^{-12}n^2 + 2.13+10^{-5}n - 0.13$ | $f(n) = 1.90 \cdot 10^{-6} n \ln n - 0.007$ |
| rmse       | 4.78    | 4.65     | 4.72   | 
| r^2        |  0.904   | 0.910    | 0.907  |

### Discussion
The only algorithm with a clear fit is brute force, with a nearly perfect fit to the quadratic form and a much poorer fit to the other ones. This is not too surprising; an extra quadratic term enables the equation fit around any noise that may exist i.e. overfitting. Thus, in general, when a simpler model has nearly identical rmse and $r^2$ values, it is preferable. 

With this in mind, the cell list does match its linear theoretical time complexity. In both the build and search models, the rmse and $r^2$ values are quite similar. Moreover, the quadratic terms in the quadratic models are miniscule, thus likely overfitting. The main deviation from expected is the kd-tree; which has linear behavior over its build and search functions. 
This is perhaps not surprising, as from 5000 atoms to 2.5 million atoms, the ln term only increases by 1.7 times; perhaps not large enough to be distinguished from noise. The addage that the ln term is essentially constant is not contradicted here. 

As for the speed comparison, the speed may be estimated by summing the polynomials for the building and search algorithms.  

Cell list: $f(n) = 3.70 \cdot 10^{-4}n + 0.385 + f(n) = 2.72 \cdot 10^{-5}n - 0.32$.  

KD-tree: $f(n) = 6.64 \cdot 10^{-5} n \ln n + 11.05 + 1.26 \cdot 10^{-5} n \ln n + 0.64$.  

Comparing values from 5000-1,000,000, the cell list outperforms the kd-tree by roughly 3-7 times.

## Usage
A few steps: clone the repository, then download all the required libraries. The desired order to run the files is: get_pdb_files.py (choosing the ranges as desired), download_pdb_files.py, save_pdb_files.cpp, loader.cpp, then last download_pdb_files.py.

For loader.cpp, if brute force is used, add that to the enum (Algorithm::BruteForce) and then restrict the max_atoms to 25,000. Once complete, remove bruteforce from the enum, then have (min_atoms, max_atoms) be (25,000, 10^7). 

## Limitations and Future Work
The results from the benchmark should neither be considered typical for proteins or for the data structures in general; there are other cases where kd-tree does outperform cell list. Moreover, a more realistic usage of these data structures would involve moving, not static proteins, although this adds technical complexity and a scenario out of scope of this project. 

The validation is not fully complete either; as it only scales to 25,000 atoms. There is a slim possibility that there exists a coding error that produces the exact same error in the kd-tree and the cell list beyond 25,000 atoms. A more rigorous validation process would involve using an external library for those two data structures, as they would be fast enough across the whole atom range. 
In addition, a more analysis of time and space complexity of the time and space usage of the kd-tree and cell list would be preferable to the rough estimation made using the polynomials. 

The GNN is technically complete, with the graph object and data pipeline mostly complete as well; however, there is no script for downloading the necessary pdb files, and the model has not be evaluated or tested for bugs. Thus, this side of the project remains incomplete. 

## Credit
Sarparast, S., Zaimi, A., Ebert, M., Goldsmith, M.R. (2024). _Advanced Atom-level Representations For Protein
Flexibility Prediction Utilizing Graph Neural Networks_. Arxiv. https://doi.org/10.48550/arXiv.2408.12519  
This paper was what this project initially attempted to replicate (the GNN and data pipeline). The idea and structure of the GNN, tensors, and data parsing came from this paper.  

Thanks to Gemmi for processing the .pdb files. https://gemmi.readthedocs.io/en/latest/
Of course thanks to the pdb databank and pdb rcsb for hosting the files publicly and enabling convenient downloads respectively. https://www.rcsb.org/ and https://data.rcsb.org/.  

Thanks to all the other libraries, which made processing, fitting, and displaying the data possible. 
