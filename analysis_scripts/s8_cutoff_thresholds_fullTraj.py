from __future__ import unicode_literals
import mdtraj as md

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# for creating native contact map
from matplotlib import colors
import matplotlib.patches as mpatches

# for analyzing frequencies
from scipy import stats

# pass in data structures for counting
def network_size(f_param):
    lst_sigma = []
    lst_network_size = []
    with open(f_param,'r') as f:
        for line in f:
            sigma = float(line.split(', ')[0].split('(')[1])
            size = float(line.split(', ')[1].split(')')[0])
            lst_sigma.append(sigma)
            lst_network_size.append(size)
    return lst_sigma, lst_network_size

# initialize variables
f_params = [
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/Ising_Model_all/network_size_stable_2JOF_all.dat',
    'BBA/Combined_1FME/Ising_Model_all/network_size_stable_BBA_all.dat',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/Ising_Model_all/network_size_stable_villin_all.dat',
    'WW domain/Combined_GTT/Ising_Model_all/network_size_stable_ww_domain_all.dat',
    'Ntl9/Combined_NTL9/Ising_Model_all/network_size_stable_Ntl9_all.dat',
    'Protein G/Combined_NuG2/Ising_Model_all/network_size_stable_NuG2_all.dat',
    'A3D/Combined_A3D/Ising_Model_all/network_size_stable_A3D_all.dat',
    'Ubiquitin/Combined_1UBQ/Ising_Model_10ns_filter_all/network_size_stable_ubq.dat',
    'l-repressor/Combined_lambda/Ising_Model_all/network_size_stable_lambda_all.dat']
lst_expected_resid = [12,32,50,50,60,103,146,154,164]
lst_names = ['Trp-cage','BBA','Villin','WW Domain','NTL9','NuG2',r'$\alpha$3D','Ubiquitin',r'$\lambda$-repressor']

# initialize lists
lst_sigmas = []
lst_network_sizes = []
for i in range(len(f_params)):
    sig, size = network_size(f_params[i])
    lst_sigmas.append(sig)
    lst_network_sizes.append(size)
    print(i)

# reset tick label size
plt.rcParams['xtick.labelsize']=40
plt.rcParams['ytick.labelsize']=40
# create subplot of all the figures
fig, ax = plt.subplots(3,3,figsize=(35,30))
# graph each figure individually
for i in range(9):
    ax[int(i/3)][i%3].plot(lst_sigmas[i],lst_network_sizes[i],linewidth=6,c='blue')
    ax[int(i/3)][i%3].axhline(y=lst_expected_resid[i],color='r',linewidth=6)
    ax[int(i/3)][i%3].tick_params(labelsize=50)
    ax[int(i/3)][i%3].set_title(lst_names[i],fontsize=56)

# create legend
fig.supylabel('Largest Cooperative Network',fontsize=56)
fig.supxlabel('Standard Deviations above the Mean',x=0.52,fontsize=56)
plt.subplots_adjust(bottom = 0.09, top=0.94, left=0.11,right = 0.97,wspace=0.25,hspace=0.25)
out_f = "s8_network_size_cutoffs_fullTraj.png"
plt.savefig(out_f)
plt.close()
