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

# pass in data structures for counting
def num_coop_interactions(f_param,f_template,numResid):
    lower = 1
    upper = 100
    # function for sorting networks size
    def getValue(e):
        return e['value']
    # make initial list of contacts
    lstPotContact = []
    for i in range(0,numResid):
        for j in range(i+1, numResid):
            if (abs(i-j) > 3):
                lstPotContact.append([i,j])
    nameContacts = []
    for j in range(len(lstPotContact)):
        nameContacts.append(lstPotContact[j][0]*numResid + lstPotContact[j][1])
    name2pair = {}
    for j in range(len(lstPotContact)):
        name = lstPotContact[j][0]*numResid + lstPotContact[j][1]
        nameContacts.append(name)
        name2pair[name] = (lstPotContact[j][0],lstPotContact[j][1])
    # calculate the native contacts 
    traj = md.load(f_template)
    prot = traj.atom_slice(traj.topology.select("protein"))
    # calculate distance between contacts (distance in nanometers)
    dist, index_contact = md.compute_contacts(prot,lstPotContact,scheme='closest-heavy',periodic=False)
    dist = dist[0]
    # filter by distance 
    NATIVE_CUTOFF = 0.45
    lstNoContact = []
    lstContact = []
    for i in range(0,len(dist)):
        if (dist[i] < NATIVE_CUTOFF):
            lstContact.append(index_contact[i])
        else:
            lstNoContact.append(index_contact[i])
    # load name of native contacts
    lstNameNative = []
    for i in range(len(lstContact)):
        lstNameNative.append(lstContact[i][0]*numResid + lstContact[i][1])
    lstNameNonnative = []
    for i in range(len(lstNoContact)):
        lstNameNonnative.append(lstNoContact[i][0]*numResid + lstNoContact[i][1])
    # make list of J matrix values and get cutoff values
    lst_j = []
    f_in = f_param
    with open(f_in,'r') as f:
        for line in f:
            if line[0] == 'J':
                value = float(line.split(' ')[5].split('\n')[0])
                lst_j.append(value)
    CUT1 = np.average(lst_j) + lower*np.std(lst_j)
    if upper == 100:
        CUT2 = max(lst_j) + 1
    else:
        CUT2 = np.average(lst_j) + upper*np.std(lst_j)
    # make list of high contacts
    lstCoop = []
    with open(f_in,'r') as f:
        for line in f:
            if line[0] == 'J':
                c1 = int(line.split(' ')[1])
                c2 = int(line.split(' ')[2])
                value = float(line.split(' ')[5].split('\n')[0])
                # restricts for cutoff value and replicates
                if (value >= CUT1) and (value < CUT2) and (c1 != c2):
                    lstCoop.append((nameContacts[c1],nameContacts[c2],value))
                    lstCoop.append((nameContacts[c2],nameContacts[c1],value))
    # create adjacency list from the list of highly cooperative contacts
    adjList = defaultdict(list)
    for i in range(len(lstCoop)):
        adjList[lstCoop[i][0]].append(lstCoop[i][1])
    # Create the map of number of edges for each residue
    contact_coop = np.zeros((numResid,numResid))
    for contact in adjList.keys():
        resid1 = name2pair[contact][0]
        resid2 = name2pair[contact][1]
        value = len(adjList[contact])
        contact_coop[resid1][resid2] = value
    # Find minimum and maximum of contact graphs
    min_coopSum = np.min(contact_coop)
    max_coopSum = np.max(contact_coop)
    # load in the native contacts on the other side
    for contact in lstNameNative:
        resid1 = name2pair[contact][0]
        resid2 = name2pair[contact][1]
        contact_coop[resid2][resid1]=contact_coop[resid1][resid2]
    mask = np.zeros((numResid,numResid))
    for i in range(numResid):
        for j in range(numResid):
            if (contact_coop[i][j] == 0):
                mask[i][j] = int(1)
    values_fe = np.ma.masked_array(contact_coop/max_coopSum,mask)
    min_vertices = np.min(contact_coop)
    max_vertices = np.max(contact_coop)
    return values_fe, min_vertices, max_vertices

# initialize variables
f_params = [
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/Ising_Model/trp-cage_stable.bm_param',
    'BBA/Combined_1FME/Ising_Model/bba_stable.bm_param',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/Ising_Model/villin_stable.bm_param',
    'WW domain/Combined_GTT/Ising_Model/ww_domain_stable.bm_param',
    'Ntl9/Combined_NTL9/Ising_Model/ntl9_stable.bm_param',
    'Protein G/Combined_NuG2/Ising_Model/nug2_stable.bm_param',
    'A3D/Combined_A3D/Ising_Model/a3d_stable.bm_param',
    'Ubiquitin/Combined_1UBQ/Ising_Model/ubiquitin_stable_cont.bm_param',
    'l-repressor/Combined_lambda/Ising_Model/lambda_stable.bm_param']
f_templates = [
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/Ising_Model/2JOF_native.pdb',
    'BBA/Combined_1FME/Ising_Model/1FME_native.pdb',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/Ising_Model/template_2F4K.pdb',
    'WW domain/Combined_GTT/Ising_Model/template_GTT.pdb',
    'Ntl9/Combined_NTL9/Ising_Model/NTL9_native.pdb',
    'Protein G/Combined_NuG2/Ising_Model/template_NuG2.pdb',
    'A3D/Combined_A3D/Ising_Model/A3D_native.pdb',
    'Ubiquitin/Combined_1UBQ/Ising_Model/template_1ubq.pdb',
    'l-repressor/Combined_lambda/Ising_Model/lambda_native.pdb']
lst_nresid = [20,28,35,35,39,56,73,76,80]
lst_coop_cutoff = [7.6,7.4,6.4,10.6,8.5,13.5,15.9,17.8,11.0]
lst_alpha = [.7,.5,.2,.2,.1,.1,.1,.03,.05]
lst_names = ['Trp-cage','BBA','Villin','WW Domain','NTL9','NuG2',r'$\alpha$3D','Ubiquitin',r'$\lambda$-repressor']
fname = ['2JOF','bba','villin','GTT','ntl9','nug2','a3d','ubq','lambda']
# initialize lists
# get number contacts matrix for each protein
lst_num_coop = []
lst_min_num_coop = []
lst_max_num_coop = []
for num in range(len(f_params)):
    contact_matrix, min_num_coop, max_num_coop = num_coop_interactions(f_params[num],f_templates[num],lst_nresid[num])
    lst_num_coop.append(contact_matrix)
    lst_min_num_coop.append(min_num_coop)
    lst_max_num_coop.append(max_num_coop)
    print(num)

out_folder = '00_Figures/s9_num_coop_interactions/'

# reset tick label size
plt.rcParams['xtick.labelsize']=40
plt.rcParams['ytick.labelsize']=40
# create subplot of all the figures
fig, ax = plt.subplots(3,3,figsize=(28,25))
# graph each figure individually
for num in range(9):
    numResid = lst_nresid[num]
    # Output contact potential distribution to a 2D discrete plot
    # create graph and colorbar
    fig,ax1 = plt.subplots(figsize=(10,9))
    bounds = np.linspace(lst_min_num_coop[num],lst_max_num_coop[num],14)
    pos = ax1.imshow(lst_num_coop[num],cmap='coolwarm',extent=(0,numResid*1000,0,numResid*1000))
    cbar = fig.colorbar(pos,ax=ax1)
    cbar.ax.set_xlabel("Frames",fontsize=36)
    cbar.ax.tick_params(labelsize = 36)
    # set ylabel parameters
    if numResid < 50:
        ax1.set_yticks(np.arange(numResid*1000-4500,499,-5000,dtype='int'))
        ax1.set_yticklabels(np.arange(5,numResid+1,5,dtype='int'),fontsize=40)
        # set xlabel parameters
        ax1.set_xticks(np.arange(4500,numResid*1000+500,5000,dtype='int'))
        ax1.set_xticklabels(np.arange(5,numResid+1,5,dtype='int'),fontsize=40)
    else:
        ax1.set_yticks(np.arange(numResid*1000-9500,499,-10000,dtype='int'))
        ax1.set_yticklabels(np.arange(10,numResid+1,10,dtype='int'),fontsize=40)
        # set xlabel parameters
        ax1.set_xticks(np.arange(9500,numResid*1000+500,10000,dtype='int'))
        ax1.set_xticklabels(np.arange(10,numResid+1,10,dtype='int'),fontsize=40)
    # plot all the connections
    # set ylabel parameters
    plt.ylabel('Residue',fontsize=40)
    ax1.yaxis.set_label_coords(-.17, .5)
    # set xlabel parameters
    plt.xlabel('Residue',fontsize=40)
    ax1.xaxis.set_label_coords(.5, -.14)
    # set title
    plt.title(lst_names[num],fontsize=40)
    plt.subplots_adjust(bottom = 0.12, top=0.94, left=0.18,right = 0.96,wspace=0.25,hspace=0.25)
    out_f = out_folder + fname[num] + '_num_coop.png'
    plt.savefig(out_f)
    plt.close()
