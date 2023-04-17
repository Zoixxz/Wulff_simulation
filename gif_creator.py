import os
import sys
import numpy as np
import imageio
print(sys.version)

file_str = os.getcwd()


G_dir_str = '/G_figures/'
dir_str = '_G_simulation_N50m35iterations4000001p0.7EXP'  # name of directory
stri = file_str + G_dir_str + dir_str
print(stri)
frames = []

names = np.arange(0, 4000001, 100000)
for t in names:
    image = imageio.imread(stri + '/' + str(t) + '.png')
    frames.append(image)

imageio.mimsave(stri + '/sim_4x10e6_p07.gif',  # output gif
                frames,
                fps=1)         # optional: frames per second
