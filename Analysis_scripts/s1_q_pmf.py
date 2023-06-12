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

# function to extract q values
def q(f_param_q):
    f_in = f_param_q
    q = []
    time = []
    count = 0.0
    with open(f_in) as f:
        for line in f:
            lstLine = line.split(' ')
            newLine = []
            for i in range(len(lstLine)):
                if (lstLine[i] != ''):
                    newLine.append(lstLine[i])
            q.append(float(newLine[1]))
            count += 0.0002
            time.append(count)
    return q,time

# function to extract unfolded state boundaries
def dual_cutoff(f_in):
    unfolded = []
    folded = []
    with open(f_in) as f:
        for line in f:
            unfolded.append(int(line.split(':')[0]))
            folded.append(int(line.split(':')[1].split('\n')[0]))
    return unfolded, folded    

# initialize variables
f_q = [
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Q/2JOF-0-protein_q.dat',
    'BBA/Combined_1FME/Q/q_all_best.dat',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Q/2F4K-0-protein_q.dat',
    'WW domain/Combined_GTT/Q/q_all_best.dat',
    'Ntl9/Combined_NTL9/Q/q_all_best.dat',
    'Protein G/Combined_NuG2/Q/q_all_best.dat',
    'A3D/Combined_A3D/Q/q_all_best.dat',
    'Ubiquitin/Q/q_all_best.dat',
    'l-repressor/Combined_lambda/Q/q_all_best.dat',
]
f_boundaries = [
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Q/unfolded_intervals.dat',
    'BBA/Combined_1FME/Q/unfolded_intervals.dat',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Q/unfolded_intervals.dat',
    'WW domain/Combined_GTT/Q/unfolded_intervals.dat',
    'Ntl9/Combined_NTL9/Q/unfolded_intervals.dat',
    'Protein G/Combined_NuG2/Q/unfolded_intervals.dat',
    'A3D/Combined_A3D/Q/unfolded_intervals.dat',
    'Ubiquitin/Q/unfolded_intervals.dat',
    'l-repressor/Combined_lambda/Q/unfolded_intervals.dat',
]
lst_names = ['Trp-cage','BBA','Villin','WW Domain','NTL9','NuG2',r'$\alpha$3D','Ubiquitin',r'$\lambda$-repressor']
fname = ['2JOF','bba','villin','GTT','ntl9','nug2','a3d','ubq','lambda']
unfolded_q = [0.10,0.05,0.28,0.05,0.12,0.30,0.19,0.10,0.40]
folded_q = [0.95,0.77,0.94,0.92,0.96,0.92,0.70,0.95,0.90]

# get the q values for each protein
lst_q = []
lst_t = []
for num in range(len(f_q)):
    q_prot, t_prot = q(f_q[num])
    lst_q.append(q_prot)
    lst_t.append(t_prot)
    print(num)

# calculate and plot PMF
lst_freq = []
for num in range(len(lst_q)):
    data = plt.hist(lst_q[num],bins=np.linspace(0,1,51),density='True')
    freq = data[0]
    for i in range(len(freq)):
        freq[i] = -1* np.log(freq[i])
    lst_freq.append(freq)

qVal = np.linspace(0.01,0.99,50)

# reset tick label size
plt.rcParams['xtick.labelsize']=40
plt.rcParams['ytick.labelsize']=40
# create subplot of all the figures
fig, ax = plt.subplots(3,3,figsize=(35,30))
# graph each figure individually
for i in range(9):
    ax[int(i/3)][i%3].plot(qVal,lst_freq[i],linewidth=6,c='blue')
    ax[int(i/3)][i%3].fill_between(qVal,lst_freq[i],min(lst_freq[i]))
    ax[int(i/3)][i%3].axvline(x=unfolded_q[i],color='r',linewidth=6,linestyle='dashed')
    ax[int(i/3)][i%3].axvline(x=folded_q[i],color='r',linewidth=6,linestyle='dashed')
    ax[int(i/3)][i%3].tick_params(labelsize=50)
    ax[int(i/3)][i%3].set_title(lst_names[i],fontsize=56)

# create legend
fig.supxlabel('Q',x=0.515,fontsize=60)
plt.subplots_adjust(bottom = 0.09, top=0.94, left=0.06,right = 0.97,wspace=0.25,hspace=0.25)
out_f = "s1_pmf.png"
plt.savefig(out_f)
plt.close()

# plot all Q values 
# get the unfolded and folded boundaries for each protein
lst_unfolded = []
lst_folded = []
for num in range(len(f_boundaries)):
    u,f = dual_cutoff(f_boundaries[num])
    lst_unfolded.append(u)
    lst_folded.append(f)

out_folder = '00_Figures/s1_q/'

for num in range(9):
    # plot average MJ Contact
    fig, ax1 = plt.subplots(figsize=(32,14))
    plt.plot(lst_t[num],lst_q[num],linewidth=1.5,color='blue')
    for i in range(len(lst_unfolded[num])):
        u = lst_unfolded[num][i]
        f = lst_folded[num][i]
        plt.plot(lst_t[num][u:f],lst_q[num][u:f],linewidth=1.5,color='orange')
    blue_patch = mpatches.Patch(color='blue', label='Folded')
    orange_patch = mpatches.Patch(color='orange',label='Unfolded')
    plt.legend(handles=[blue_patch,orange_patch],loc='upper left',fontsize=40,borderaxespad=0.10)
    plt.xlabel(r'Time ($\mu$s)', fontsize=60)
    ax1.xaxis.set_label_coords(.5, -.12)
    plt.ylabel('Q', fontsize=60)
    ax1.yaxis.set_label_coords(-.07, .5)
    plt.subplots_adjust(bottom = 0.17, top=0.90, left=0.1,right = 0.97,wspace=0.25,hspace=0.25)
    plt.title(lst_names[num], fontsize=60)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=60)
    out_f = out_folder + fname[num] + '_q.png'
    plt.savefig(out_f)
    plt.close('all')
