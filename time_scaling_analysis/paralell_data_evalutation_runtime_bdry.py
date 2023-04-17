import os
import sys
import numpy as np
from scipy.optimize import curve_fit
import igraph as ig
import matplotlib.pyplot as plt
from matplotlib import cm
import ig_python_simulation_wulff_runtime as sim
from joblib import Parallel, delayed
import multiprocessing
from datetime import datetime
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

"""
this file is dedicated to analysing the scaling of the time steps of the code
"""
p = 0.6  # parameter for percolation

number_of_simulations = 15  # number of simulations to be performed

size_simulation_box = np.array([20, 26, 36, 45, 60, 75, 90,
                                105, 120, 135, 150, 170, 190, 210, 230])  # N
size_simulation_cluster = np.array([13, 18, 24, 30, 40, 50, 60,
                                    70, 80, 90, 100, 115, 130, 150, 170])  # m
percolation_parameter = np.full(number_of_simulations, p, dtype=float)
iterations_zip = np.array([1000, 1000, 1000, 1000, 100, 100, 100,
                           50, 50, 50, 50, 50, 50, 50, 20])  # sample size

# create folder for plots
folderName = ('runtime_evaluation/')
sim.createFolder(folderName)

# check the number of available cores for parallelisation
num_cores = multiprocessing.cpu_count()

# build iterator for joblib Parallel
iterator = zip(size_simulation_box, size_simulation_cluster,
               percolation_parameter, iterations_zip)

print('number of cores ', num_cores)
"""
note: tqdm only tracks queueing progress as mentioned here
https://stackoverflow.com/questions/24983493
"""

results = Parallel(n_jobs=num_cores)(delayed(sim.performSimulation)(v) for v in tqdm(iterator))

"""
this way of computing is here for testing purposes
results = [sim.performSimulation(v) for v in iterator]
"""
timeData = []
for i in range(number_of_simulations):
    timeData.append(results[i])

print(timeData)
means = np.zeros(number_of_simulations)
sterr = np.zeros(number_of_simulations)

for i in range(number_of_simulations):
    means[i] = np.mean(timeData[i])
    sterr[i] = np.std(timeData[i]) / np.sqrt(len(timeData[i]))
# typically the standard error will be too small to show on plot

"""
plot the time step analysis
"""
x_axis = [(2 * (v[0] - 2)) * ((2 * v[1] + 1) **2)
          for v in zip(size_simulation_cluster, size_simulation_box)]
"""
# fitting polynomial
coeffs1 = np.polyfit(x_axis, means, deg=1)
poly1 = np.poly1d(coeffs1)
"""
"""
fit a linear function going through 0
"""


def fit_fun(x, a):
    return a * x


params = curve_fit(fit_fun, x_axis, means)
[a] = params[0]
y_axis = [fit_fun(x, a) for x in x_axis]
print('coeff ', a)
print('x_axis ', x_axis)


fig, (ax1) = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(8, 4))
ax1.plot(x_axis, y_axis, color='r', linestyle='--', label=f'linear order')
ax1.scatter(x_axis, means, color='k', marker='o', label='time per iteration')
ax1.legend()

ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.tick_params(labelsize=11.5)
ax1.set_title(r'scaling of timestep w. r. to $\beta \cdot |V|$', fontsize=14)
ax1.set_xlabel(r'$\beta\cdot |V|$', fontsize=14)
ax1.set_ylabel('average time per iteration [s]', fontsize=14)
fig.tight_layout()
plt.savefig(folderName + '/' + 'runtime_combined' + '.pdf')
plt.clf()
plt.close()
