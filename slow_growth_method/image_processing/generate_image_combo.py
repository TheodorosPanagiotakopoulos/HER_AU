import data
import os
import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
from ase.io import read
import json as js
import numpy as np
from figure import get_figure_2

os.system("rm -rf frame_*")
loc = os.getcwd() + "/"
dirs = data.get_dirs()
dirs = [loc + x for x in dirs]
count = -1
step = 250
start = -1

idx1, idx2, idx3 = data.get_standarized_ICONST_data(loc)
print(idx1, idx2, idx3)

system = read(loc + '/RUN1/POSCAR')

# Define the rotation angles
rotations = [ "180z, -60x, 30y", "180z, -60x, 40y", "180z, -60x, 50y", "180z, -60x, 60y", "180z, -60x, 70y",  "180z, -70x, 30y", "180z, -70x, 40y", "180z, -70x, 50y", "180z, -70x, 60y", "180z, -70x, 70y", "180z, -80x, 30y", "180z, -80x, 40y", "180z, -80x, 50y", "180z, -80x, 60y", "180z, -80x, 70y", "180z, -60x, -30y", "180z, -60x, -40y", "180z, -60x, -50y", "180z, -60x, -60y", "180z, -60x, -70y",  "180z, -70x, -30y", "180z, -70x, -40y", "180z, -70x, -50y", "180z, -70x, -60y", "180z, -70x, -70y", "180z, -80x, -30y", "180z, -80x, -40y", "180z, -80x, -50y", "180z, -80x, -60y", "180z, -80x, -70y" ]

for dir in dirs:
    with open(dir + '/DATA.js', 'r') as f:
        data = js.load(f)
    npoints = len(data.keys())
    traj = Trajectory(dir + '/MOVIE.traj', 'r')
    for i in range(len(traj)):
        count += 1
        if count > start and count % step == 0:
            print(count)
            img = traj[i]
            ref = img[41].position * 1.
            img[184].symbol = 'He'
            img[185].symbol = 'He'
            img[186].symbol = 'He'
            img[187].symbol = 'He'

            img[idx1].symbol = 'He'
            if system[idx2].symbol != "N":
                img[idx2].symbol = 'S'

            if not np.isnan(idx3):
                img[idx3].symbol = "He"

            center = (img.cell[0, :] + img.cell[1, :]) / 2.
            center[2] = 0.
            img.positions[:, 0] -= ref[0]
            img.positions[:, 1] -= ref[1]
            img.positions[:, 0] += center[0]
            img.positions[:, 1] += center[1]
            img.wrap()
            prefix = str(count)

            for rot in rotations:
                rot_filename = rot.replace(", ", "_").replace(" ", "_")
                fout = f'frame_{prefix}_{rot_filename}'
                if not os.path.isfile(fout + '.png'):
                    get_figure_2(img, fout, rot=rot, w=13, h=12)
            os.system('rm -f *.ini *.pov')
