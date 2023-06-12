# Notes: Use "conda activate md" instead of "conda activate figures" or colorbar axes labels will not appear
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

def num_contacts_tp(f_in,f_template,numResid):
    # make initial list of contacts
    lstPotContact = []
    for i in range(0,numResid):
        for j in range(i+1, numResid):
            if (abs(i-j) > 3):
                lstPotContact.append([i,j])
    nameContacts = []
    name2pair = {}
    for j in range(len(lstPotContact)):
        name = lstPotContact[j][0]*numResid + lstPotContact[j][1]
        nameContacts.append(name)
        name2pair[name] = (lstPotContact[j][0],lstPotContact[j][1])
    # Find native contacts and native contact names
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
    # initialize number contacts values 
    numContacts = {}
    lst_num = []
    with open(f_in,'r') as f:
        for line in f:
            c1 = int(line.split('-')[0].split(' ')[1])
            c2 = int(line.split('-')[1].split(' ')[0])
            key = (c1,c2)
            value = int(float(line.split(' ')[2].split('\n')[0]))
            numContacts[key] = value
            if value > 0:
                lst_num.append(value)
    # Create cooperativity map
    freqMap = np.zeros((numResid,numResid))
    for key in numContacts.keys():
        freqMap[key[1]][key[0]] = numContacts[key]
    # Fill in native contacts
    for contact in lstNameNative:
        resid1 = name2pair[contact][0]
        resid2 = name2pair[contact][1]
        freqMap[resid1][resid2]=freqMap[resid2][resid1]
    # Create mask
    mask = np.zeros((numResid,numResid))
    for i in range(numResid):
        for j in range(numResid):
            if (freqMap[i][j] == 0):
                mask[i][j] = int(1)
    # find the median number of contacts
    min_freq = np.min(lst_num)
    # create a log of the frequency map
    fe_map = np.zeros((numResid,numResid))
    for i in range(len(freqMap)):
        for j in range(len(freqMap[0])):
            if freqMap[i][j] != 0:
                fe_map[i][j] = freqMap[i][j]/1000
            else:
                fe_map[i][j] = 0
    values_fe = np.ma.masked_array(fe_map,mask)
    return values_fe

def num_contacts_tp_native(f_in,f_template,numResid):
    # make initial list of contacts
    lstPotContact = []
    for i in range(0,numResid):
        for j in range(i+1, numResid):
            if (abs(i-j) > 3):
                lstPotContact.append([i,j])
    nameContacts = []
    name2pair = {}
    for j in range(len(lstPotContact)):
        name = lstPotContact[j][0]*numResid + lstPotContact[j][1]
        nameContacts.append(name)
        name2pair[name] = (lstPotContact[j][0],lstPotContact[j][1])
    # Find native contacts and native contact names
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
    # initialize number contacts values 
    numContacts = {}
    lst_num = []
    with open(f_in,'r') as f:
        for line in f:
            c1 = int(line.split('-')[0].split(' ')[1])
            c2 = int(line.split('-')[1].split(' ')[0])
            key = (c1,c2)
            value = int(float(line.split(' ')[2].split('\n')[0]))
            numContacts[key] = value
            c = key[0]*numResid+key[1]
            if (value > 0) and (c in lstNameNative):
                lst_num.append(value)
    # find the median number of contacts
    max_freq = np.max(lst_num)
    min_freq = np.min(lst_num)
    # Create cooperativity map
    freqMap = np.zeros((numResid,numResid))
    for key in numContacts.keys():
        c = key[0]*numResid+key[1]
        if c in lstNameNative:
            freqMap[key[1]][key[0]] = numContacts[key]
    # Fill in native contacts
    for contact in lstNameNative:
        resid1 = name2pair[contact][0]
        resid2 = name2pair[contact][1]
        freqMap[resid1][resid2]=max_freq
    # Create mask
    mask = np.zeros((numResid,numResid))
    for i in range(numResid):
        for j in range(numResid):
            if (freqMap[i][j] == 0):
                mask[i][j] = int(1)
    # create a log of the frequency map
    fe_map = np.zeros((numResid,numResid))
    for i in range(len(freqMap)):
        for j in range(len(freqMap[0])):
            if freqMap[i][j] != 0:
                fe_map[i][j] = freqMap[i][j]/1000
            else:
                fe_map[i][j] = 0
    values_fe = np.ma.masked_array(fe_map,mask)
    return values_fe

f_num_contacts = [
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/tp_num_contacts.dat',
    'BBA/Combined_1FME/tp_num_contacts.dat',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/tp_num_contacts.dat',
    'WW domain/Combined_GTT/tp_num_contacts.dat',
    'Ntl9/Combined_NTL9/tp_num_contacts.dat',
    'Protein G/Combined_NuG2/tp_num_contacts.dat',
    'A3D/Combined_A3D/tp_num_contacts.dat',
    'Ubiquitin/Combined_1UBQ/tp_num_contacts.dat',
    'l-repressor/Combined_lambda/tp_num_contacts.dat',
]
f_templates = [
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/Ising_Model/2JOF_native.pdb',
    'BBA/Combined_1FME/Ising_Model/1FME_native.pdb',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/Ising_Model/template_2F4K.pdb',
    'WW domain/Combined_GTT/Ising_Model/template_GTT.pdb',
    'Ntl9/Combined_NTL9/Ising_Model/NTL9_native.pdb',
    'Protein G/Combined_NuG2/Ising_Model/template_NuG2.pdb',
    'A3D/Combined_A3D/Ising_Model/A3D_native.pdb',
    'Ubiquitin/Combined_1UBQ/Ising_Model/template_1ubq.pdb',
    'l-repressor/Combined_lambda/Ising_Model/lambda_native.pdb'
]
f_mj = 'Ntl9/Combined_NTL9/mj.txt'
lst_nresid = [20,28,35,35,39,56,73,76,80]
lst_names = ['Trp-cage','BBA','Villin','WW Domain','NTL9','NuG2',r'$\alpha$3D','Ubiquitin',r'$\lambda$-repressor']
fname = ['2JOF','bba','villin','GTT','ntl9','nug2','a3d','ubq','lambda']

# get number contacts in transition path for each protein
lst_numContacts = []
for num in range(len(f_num_contacts)):
    contact_matrix = num_contacts_tp(f_num_contacts[num],f_templates[num],lst_nresid[num])
    lst_numContacts.append(contact_matrix)
    print(num)

# get number native contacts in transition path for each protein
lst_numContacts_native = []
for num in range(len(f_num_contacts)):
    contact_matrix = num_contacts_tp_native(f_num_contacts[num],f_templates[num],lst_nresid[num])
    lst_numContacts_native.append(contact_matrix)
    print(num)

out_folder = '00_Figures/s8_tp_fe/'

plt.rcParams['xtick.labelsize']=40
plt.rcParams['ytick.labelsize']=40
# plot all numContacts FE transition path plots
for i in range(9):
    numResid = lst_nresid[i]
    # Output contact potential distribution to a 2D discrete plot
    # create graph and colorbar
    fig,ax1 = plt.subplots(figsize=(12,10))
    pos = ax1.imshow(lst_numContacts[i],cmap='coolwarm',extent=(0,numResid*1000,0,numResid*1000))
    cbar = fig.colorbar(pos,ax=ax1,pad=0.1)
    cbar.ax.set_xlabel("Frames" + "\n" + r"(10^3)",fontsize=36)
    cbar.ax.tick_params(labelsize = 40)
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
    # set ylabel parameters
    plt.ylabel('Residue',fontsize=40)
    ax1.yaxis.set_label_coords(-.14, .5)
    # set xlabel parameters
    plt.xlabel('Residue',fontsize=40)
    ax1.xaxis.set_label_coords(.53, -.11)
    # set title
    plt.title(lst_names[i] + " Transition Path",fontsize=40)
    plt.subplots_adjust(bottom = 0.12, top=0.94, left=0.14,right = 0.98,wspace=0.25,hspace=0.25)
    out_f = out_folder + fname[i] + '_tp.png'
    plt.savefig(out_f)
    plt.close()

# plot all numContacts FE transition path native plots
for i in range(9):
    numResid = lst_nresid[i]
    # Output contact potential distribution to a 2D discrete plot
    # create graph and colorbar
    fig,ax1 = plt.subplots(figsize=(12,10))
    pos = ax1.imshow(lst_numContacts_native[i],cmap='coolwarm',extent=(0,numResid*1000,0,numResid*1000))
    cbar = fig.colorbar(pos,ax=ax1,pad=0.1)
    cbar.ax.set_xlabel("Frames" + "\n" + r"(10^3)",fontsize=36)
    cbar.ax.tick_params(labelsize = 40)
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
    # set ylabel parameters
    plt.ylabel('Residue',fontsize=40)
    ax1.yaxis.set_label_coords(-.14, .5)
    # set xlabel parameters
    plt.xlabel('Residue',fontsize=40)
    ax1.xaxis.set_label_coords(.53, -.11)
    # set title
    plt.title(lst_names[i],fontsize=40)
    plt.subplots_adjust(bottom = 0.12, top=0.94, left=0.14,right = 0.98,wspace=0.25,hspace=0.25)
    out_f = out_folder + fname[i] + '_tp_native.png'
    plt.savefig(out_f)
    plt.close()
