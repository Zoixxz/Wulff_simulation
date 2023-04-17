# Wulff_simulation
Code for the simulating the Wulff shape on the $\mathbb{Z}^2$ lattice.
The library used for the simulation is iGraph which is implemented in C. For a performance comparison between different libraries see https://graph-tool.skewed.de/performance.

How does the algorithm work?

The main idead is to approximate the measure which we are interested in, namely $P_p( \cdot  | |C(0)| = m)$, using a Markov chain which converges to this distribution. To achieve this we use Metropolis transition probabilities.

What is in the folders?

- The folder "reversible_MCs_simulation" contains the files used to simulate large finite clusters for supercitical bernoulli bond percolation on the square lattice.
- The folder "time_scaling_analysis" contains adapted files to facilitate the analysis of the time scaling behaviour of the code.
- The folder "nx_wulff_simulation" is NOT MAINTAINED and contains a preliminary version of the simulation using the NetworkX library which is significantly slower than iGraph.

"gif_creator.py" is a short pyhton code which can be used to create GIFs using the clusters generated via simulation.
