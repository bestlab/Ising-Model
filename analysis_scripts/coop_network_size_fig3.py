from __future__ import unicode_literals
import mdtraj as md

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# for creating cooperative network
from collections import defaultdict
from copy import copy

# for creating native contact map
from matplotlib import colors
import matplotlib.patches as mpatches

# for analyzing frequencies
from scipy import stats

numResid = 39

# loops through all DESRES trajectories
lst_sigma = []
lst_network_size = []
f_in = 'network_size_stable_Ntl9.dat'
with open(f_in,'r') as f:
    for line in f:
        sigma = float(line.split(', ')[0].split('(')[1])
        size = float(line.split(', ')[1].split(')')[0])
        lst_sigma.append(sigma)
        lst_network_size.append(size)

plt.rcParams['xtick.labelsize'] = 40
plt.rcParams['ytick.labelsize'] = 40
# create plot
fig,ax1 = plt.subplots(figsize=(14,10))
plt.plot(lst_sigma,lst_network_size,linewidth=6,c='blue')
plt.axhline(y=60,color='r',linewidth=6)
plt.ylabel('Largest\n' + 'Cooperative Network',linespacing=1.5,fontsize=40)
ax1.yaxis.set_label_coords(-.13, .5)
plt.xlabel('Threshold Cooperativity (Z)',fontsize=40)
ax1.xaxis.set_label_coords(.50, -.10)
plt.subplots_adjust(bottom = 0.13, top=0.95, left=0.20,right = 0.97,wspace=0.25,hspace=0.25)
#title = "NTL9 Cooperative Network Threshold"
#plt.title(title,fontsize=44)
out_f = "network_size_stable_NTL9_fig3.png"
plt.savefig(out_f)
plt.close()
