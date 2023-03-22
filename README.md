# Wulff_simulation
Code for the simulating the Wulff shape on the $\mathbb{Z}^2$ lattice.
There are currently two version of the same code using a different graph library. Once the library NetworkX is used which is a pure pyhton implementation and therefore does not perform too well for large graphs. The other library is iGraph which in generally is much more efficient. For a comparison see https://graph-tool.skewed.de/performance.
How does the algorithm work?
The main idead is to approximate the measure which we are interested in, namely $P_p( \cdot  | |C(0)| = m)$, using a Markov chain which converges to this distribution. To achieve this we use Metropolis transition probabilities.
