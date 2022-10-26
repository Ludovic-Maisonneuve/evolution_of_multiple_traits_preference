import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib
import os

from functions import *

SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 25

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

eps = 0.01
res = 250

L_log_a = np.linspace(-1, 1, res)
L_c = np.linspace(0, 0.5 * eps, res)

R = np.zeros((res, res))

thisdir = 'results/a'

for r, d, f in os.walk(thisdir):  # r=root, d=directories, f = files
    for file in f:
        if file != '.DS_Store':
            f = file.split('_')
            i = int(f[0])
            j = int(f[1][:-4])
            result = open(thisdir + '/' + file, 'r')
            list_results = result.read().split('\n')
            #    os.remove(thisdir + '/' + file)
            count = int(list_results[0])
            if count == 0:
                print('count == 0')
                os.remove(thisdir + '/' + file)

            L_g = list_results[1][1:-1].split(',')
            L_gamma_ES = []
            for g in L_g:
                L_gamma_ES.append(float(g))

            L_n = list_results[2][1:-1].split(',')
            L_nature = []
            for i_n, txt in enumerate(L_n):

                if i_n == 0:
                    L_nature.append(txt[1:-1])
                else:
                    L_nature.append(txt[2:-1])

            if count == 1:
                if L_nature[0] == 'x':
                    R[i,j] = L_gamma_ES[0]

            #elif count == 2:
                #print('issue: count == 2')
                #print(i,j)
                #os.remove(thisdir + '/' + file)

            elif count == 3:
                i_o = L_nature.index('o')
                gamma_o = L_gamma_ES[i_o]
                L_x = []
                for i_x in range(3):
                    if i_x != i_o:
                        L_x.append(L_gamma_ES[i_x])

                gamma_x1 = L_x[0]
                gamma_x2 = L_x[1]
                if (gamma_x1 - 1 / 2) * (gamma_x2 - 1 / 2) < 0:
                    if gamma_o > 1/2:
                        R[i, j] = np.min([gamma_x1, gamma_x2])
                    else:
                        R[i, j] = np.max([gamma_x1, gamma_x2])
                else:
                    if gamma_x1 > 1/2:
                        R[i, j] = np.min([gamma_x1, gamma_x2])
                    else:
                        R[i, j] = np.max([gamma_x1, gamma_x2])

PiYG = matplotlib.cm.get_cmap('PiYG')

L = []
for t in np.linspace(0, 1, 101):
    L.append(PiYG(t))

cm = LinearSegmentedColormap.from_list(
    'cmap', L, N=101)

cm = LinearSegmentedColormap.from_list(
    'cmap', L, N=101)

fig, ax = plt.subplots(1)
plt.imshow(R.T, extent=[-1, 1, 0, 0.5 * eps], aspect= 2 / (0.5 * eps), vmin = 0, vmax=1, origin='lower', cmap=cm)
plt.xticks([-1, 0, 1])
plt.yticks([0, 0.0025, 0.005])
#CS = ax.contour(L_log_a, L_c, R.T, levels=[0, 1], colors=[np.array([0.001,0.634, 1, 1])], linewidths=4, linestyles= 'dashed')
plt.tight_layout()