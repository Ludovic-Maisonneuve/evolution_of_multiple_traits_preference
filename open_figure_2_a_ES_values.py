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
                if 'x' in txt:
                    L_nature.append('x')
                else:
                    L_nature.append('o')

            if count == 1:
                if L_nature[0] == 'x':
                    R[i,j] = L_gamma_ES[0]

            #elif count == 2:
                #print('issue: count == 2')
                #print(i,j)
                #os.remove(thisdir + '/' + file)

            elif count == 3:
                if 0 in L_gamma_ES and 1 in L_gamma_ES:
                    R[i, j] = -1
                else:
                    #print(count, L_gamma_ES, L_nature)
                    R[i, j] = -2
                    #print(i,j)
                    #os.remove(thisdir + '/' + file)
            else:
                #print('count > 3', i, j)
                #print(count, L_gamma_ES, L_nature)
                R[i, j] = -2

PiYG = matplotlib.cm.get_cmap('PiYG')

L = []
for t in np.linspace(-2, 1, 301):
    if t < -1:
        L.append(np.array([30/255,144/255,1,1]))
    elif t < 0:
        L.append(np.array([0, 0, 0, 1]))
    else:
        L.append(PiYG(t))

cm = LinearSegmentedColormap.from_list(
    'cmap', L, N=301)

plt.figure()
plt.imshow(R.T, extent=[-1, 1, 0, 0.5 * eps], aspect= 2 / (0.5 * eps), vmin = -2, vmax=1, origin='lower', cmap=cm)
plt.xticks([-1, 0, 1])
plt.yticks([0, 0.0025, 0.005])
plt.tight_layout()