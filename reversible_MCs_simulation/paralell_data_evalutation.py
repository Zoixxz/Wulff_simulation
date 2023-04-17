import sys
import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
from matplotlib import cm
import ig_python_simulation_wulff as sim
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm

"""
the following print statements give infro on the versions of
python used and the igraph library
"""
print("Python version")
print(sys.version)
print("Version info.")
print(sys.version_info)
print("igraph version: " + ig.__version__)


N = 24  # we build a box from -N to N
m = 15  # m is size of inside box, must be strictly smaller than N
p = 0.6  # percolation parameter
iterations = 100000  # number of iterations to be performed
drawst = 20000  # steps after which to draw the graph
do_plts = True  # if False no graphs will be drawn

number_of_simulations = 1  # number of simulations to be performed
size_simulation_box = np.full(number_of_simulations, N, dtype=int)      # N
size_simulation_cluster = np.full(number_of_simulations, m, dtype=int)  # m
percolation_parameter = np.full(number_of_simulations, p, dtype=float)  # p
iterations_zip = np.full(number_of_simulations, iterations, dtype=int)  # iter
draw_step_size = np.full(number_of_simulations, drawst, dtype=int)      # dst
do_plots = np.full(number_of_simulations, do_plts, dtype=bool)          # plts
labeling_folders = np.arange(number_of_simulations)                     # j

# create iterator for joblib Parallel
iterator = zip(size_simulation_box, size_simulation_cluster,
               percolation_parameter, iterations_zip, draw_step_size,
               do_plots, labeling_folders)

num_cores = multiprocessing.cpu_count()
print('number of cores ', num_cores)
# note: tqdm only tracks queueing progress as mentioned here
# https://stackoverflow.com/questions/24983493
results = Parallel(n_jobs=num_cores)(delayed(sim.performSimulation)(v)
                                     for v in tqdm(iterator))

# create array with data
data = np.zeros((iterations, number_of_simulations))
for i in range(number_of_simulations):
    data[:, i] = results[i]

print('length results = ', len(results))

"""
plot the precentiles of the boudnary size
"""
colormap = cm.Blues  # set the colomap for the quantiles
# percentiles have to always include 0.5 and contain an odd number
percentileQuantiles = [0.1, 0.25, 0.4, 0.5, 0.6,
                       0.75, 0.9]  # percentiles for computation

n = len(percentileQuantiles)
half = int((n - 1) / 2)

folderName = ('G_figures/' + 'rev_G_simulation_N' + str(N) + 'm' + str(m)
              + 'iterations' + str(iterations) + 'p' + str(p) + '/'
              + 'percentiles')
sim.createFolder(folderName)

percData = np.zeros((iterations, len(percentileQuantiles)))
for i in range(len(percentileQuantiles)):
    for t in range(iterations):
        percData[t, i] = np.percentile(
            data[t, :], percentileQuantiles[i] * 100)

"""
here we compute the mean of the boundary size for the performed
simulations
"""
mean = np.zeros(iterations)
for t in range(iterations):
    mean[t] = np.mean(percData[t, :])
"""
plot the boundary size analysis
"""
fig, (ax1) = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(8, 4))
x_axis = np.arange(0, iterations, 1)
ax1.plot(x_axis, mean, color='red', label='mean')
ax1.plot(x_axis, percData[:, half], color='k',
         label='median')
for i in range(half):
    value = percentileQuantiles[i]
    label = str(int(value * 100)) + 'th perc.'
    ax1.fill_between(x_axis, percData[:, i], percData[:, - (i + 1)],
                     color=colormap(i * 7 / (8 * half) + 1 / 8),
                     label=label)
ax1.set_title(
    f"boundary size of {number_of_simulations} realisations", fontsize=15)
ax1.tick_params(labelsize=11.5)
ax1.set_xlabel('steps', fontsize=14)
ax1.set_ylabel('boundary size', fontsize=14)
ax1.legend()
fig.tight_layout()
plt.savefig(folderName + '/' + 'percentiles' + '.pdf')
plt.clf()
plt.close()
