# Wulff_simulation
Code for simulating large but finite clusters for Bernoulli bond percolation on the $\mathbb{Z}^2$ lattice.
The library used for the simulation is iGraph which is implemented in C. For a performance comparison between different libraries see https://graph-tool.skewed.de/performance.

How does the algorithm work?

The main idead is to approximate the measure which we are interested in using a Markov chain which converges to this distribution. To achieve this we use Metropolis transition probabilities.

What is in the folders?

- The folder "reversible_MCs_simulation" contains the files used to simulate large finite clusters for supercitical bernoulli bond percolation on the square lattice. The file also contains the transistion probabilites for another Markov chain converging to the desired distribution.
- The folder "time_scaling_analysis" contains adapted files to facilitate the analysis of the time scaling behaviour of the code.
- The folder "nx_wulff_simulation" is NOT MAINTAINED and contains a preliminary version of the simulation using the NetworkX library which is significantly slower than iGraph. Furthermore it does NOT compute the correct transition probabilities but an approximate version.

"gif_creator.py" is a short pyhton code which can be used to create GIFs using the clusters generated via simulation.

How are the files in "reversible_MCs_simulatio" interconnected?

- "paralell_data_evalutation.py" is the file containing the main function where the parameters of the simulation can be specified.
- "ig_percoSimAux.py" contains functions needed for the initalisation of a cluster using a norm. The 1-norm and max-norm are implemented.
- "ig_percoSimIteration.py" contains the function handling one iteration step.
- "ig_percoSimIterationAux.py" contains functions mainly used in the iteration step, such as computing the acceptance probability.
- "ig_python_simulation_wulff.py" handles how often the iteration is performed and the data collection and plotting

