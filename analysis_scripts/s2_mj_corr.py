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

# initialize maps for reading amino acids
amiAciToPot = {
    'ALA': 0,
    'CYS': 1,
    'ASP': 2,
    'GLU': 3,
    'PHE': 4,
    'GLY': 5,
    'HIS': 6,
    'HIP': 6,
    'ILE': 7,
    'LYS': 8,
    'LEU': 9,
    'MET': 10,
    'NLE': 10,
    'ASN': 11,
    'PRO': 12,
    'GLN': 13,
    'ARG': 14,
    'SER': 15,
    'THR': 16,
    'VAL': 17,
    'TRP': 18,
    'TYR': 19}

aaToPot = {
    'A': 0,
    'C': 1,
    'D': 2,
    'E': 3,
    'F': 4,
    'G': 5,
    'H': 6,
    'I': 7,
    'K': 8,
    'L': 9,
    'M': 10,
    'N': 11,
    'P': 12,
    'Q': 13,
    'R': 14,
    'S': 15,
    'T': 16,
    'V': 17,
    'W': 18,
    'Y': 19}

def num_contacts(f_in,f_template,numResid):
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
    values = np.ma.masked_array(freqMap,mask)
    # find the median number of contacts
    min_freq = np.min(lst_num)
    # create a log of the frequency map
    fe_map = np.zeros((numResid,numResid))
    for i in range(len(freqMap)):
        for j in range(len(freqMap[0])):
            if freqMap[i][j] != 0:
                fe_map[i][j] = -1*np.log(freqMap[i][j]/min_freq)
            else:
                fe_map[i][j] = 0
    values_fe = np.ma.masked_array(fe_map,mask)
    return values_fe

def mj(mj_pot,f_template,numResid):
    traj = md.load(f_template)
    prot = traj.atom_slice(traj.topology.select("protein"))
    # make initial list of contacts
    lstPotContact = []
    for i in range(0,numResid):
        for j in range(i+1, numResid):
            if (abs(i-j) > 3):
                lstPotContact.append([i,j])
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
    # Find contact potentials
    contactPot = np.zeros((numResid,numResid))
    # Find potential for all residue pairs
    for contact in index_contact:
        resid1 = contact[0]
        resid2 = contact[1]
        i_pot = amiAciToPot[str(prot.topology.residue(resid1))[0:3]]
        j_pot = amiAciToPot[str(prot.topology.residue(resid2))[0:3]]
        contactPot[resid2][resid1] = mj_pot[i_pot,j_pot]
    # Plot potential for native contacts
    for i in range(0,len(lstContact)):
        resid1 = lstContact[i][0]
        resid2 = lstContact[i][1]
        i_pot = amiAciToPot[str(prot.topology.residue(resid1))[0:3]]
        j_pot = amiAciToPot[str(prot.topology.residue(resid2))[0:3]]
        contactPot[resid1][resid2] = mj_pot[i_pot,j_pot]   
    # Create mask
    mask = np.zeros((numResid,numResid))
    for i in range(numResid):
        for j in range(numResid):
            if (contactPot[i][j] == 0):
                mask[i][j] = int(1)
    values = np.ma.masked_array(contactPot,mask)
    return values

def numContacts_mj_corr(f_in,mj_pot,f_template,numResid):
    traj = md.load(f_template)
    prot = traj.atom_slice(traj.topology.select("protein"))
    # make initial list of contacts
    lstPotContact = []
    for i in range(0,numResid):
        for j in range(i+1, numResid):
            if (abs(i-j) > 5):
                lstPotContact.append([i,j])
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
    contactPot = np.zeros((numResid,numResid))
    # Find potential for residue pairs labeled as native contacts
    for i in range(0,len(lstContact)):
        resid1 = lstContact[i][0]
        resid2 = lstContact[i][1]
        i_pot = amiAciToPot[str(prot.topology.residue(resid1))[0:3]]
        j_pot = amiAciToPot[str(prot.topology.residue(resid2))[0:3]]
        if (resid1 < resid2):
            contactPot[resid1][resid2] = mj_pot[i_pot,j_pot]   
        else:
            contactPot[resid2][resid1] = mj_pot[i_pot,j_pot]
    # Find potentials of residues not in contact
    for i in range(0,len(lstPotContact)):
        resid1 = lstPotContact[i][0]
        resid2 = lstPotContact[i][1]
        i_pot = amiAciToPot[str(prot.topology.residue(resid1))[0:3]]
        j_pot = amiAciToPot[str(prot.topology.residue(resid2))[0:3]]
        if (resid1 < resid2):
            contactPot[resid2][resid1] = mj_pot[i_pot,j_pot]   
        else:
            contactPot[resid1][resid2] = mj_pot[i_pot,j_pot]
    nameContacts = []
    name2pair = {}
    for j in range(len(lstPotContact)):
        name = lstPotContact[j][0]*numResid + lstPotContact[j][1]
        nameContacts.append(name)
        name2pair[name] = (lstPotContact[j][0],lstPotContact[j][1])
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
    for contact in lstNameNative:
        resid1 = name2pair[contact][0]
        resid2 = name2pair[contact][1]
        if (resid1 < resid2):
            freqMap[resid1][resid2]=freqMap[resid2][resid1]
        else:
            freqMap[resid1][resid2]=freqMap[resid1][resid2]
    # find the median number of contacts
    min_freq = np.min(lst_num)
    # create a log of the frequency map
    fe_map = np.zeros((numResid,numResid))
    for i in range(len(freqMap)):
        for j in range(len(freqMap[0])):
            if freqMap[i][j] != 0:
                fe_map[i][j] = -1*np.log(freqMap[i][j]/min_freq)
            else:
                fe_map[i][j] = 0
    lst_mj = []
    for i in range(len(mj_pot)):
        for j in range(len(mj_pot[0])):
            lst_mj.append(mj_pot[i][j])
    lst_numContacts = []
    for value in numContacts.values():
        lst_numContacts.append(value)
    # Calculate the minimum number of contacts and minimum number of 
    min_numContacts = min(numContacts.values())
    min_mj = np.max(mj_pot)
    avg_numContacts = np.average(lst_numContacts)
    avg_mj = np.average(lst_mj)
    lst_Freq = []
    lst_mjContact = []
    pairs = []
    for i in range(len(freqMap)):
        for j in range(len(freqMap[0])):
            if (i>j) and (contactPot[i][j] != 0):
                pairs.append((i,j))
                lst_Freq.append(fe_map[i][j])
                lst_mjContact.append(contactPot[i][j])
    lst_nativeFreq = []
    lst_nativeMj = []
    for contact in lstNameNative:
        resid1 = name2pair[contact][0]
        resid2 = name2pair[contact][1]
        lst_nativeFreq.append(fe_map[resid1][resid2])
        lst_nativeMj.append(contactPot[resid1][resid2])
    r,p = stats.pearsonr(lst_Freq,lst_mjContact)
    return lst_Freq,lst_mjContact,lst_nativeFreq,lst_nativeMj,r,p

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
f_num_contacts = [
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/numContacts.dat',
    'BBA/Combined_1FME/numContacts.dat',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/numContacts.dat',
    'WW domain/Combined_GTT/numContacts.dat',
    'Ntl9/Combined_NTL9/numContacts.dat',
    'Protein G/Combined_NuG2/numContacts.dat',
    'A3D/Combined_A3D/numContacts.dat',
    'Ubiquitin/Combined_1UBQ/numContacts.dat',
    'l-repressor/Combined_lambda/numContacts.dat',
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

# get number contacts matrix for each protein
lst_numContacts = []
for num in range(len(f_num_contacts)):
    contact_matrix = num_contacts(f_num_contacts[num],f_templates[num],lst_nresid[num])
    lst_numContacts.append(contact_matrix)
    print(num)

# Create MJ potential from text file
mj_params = np.zeros((20,20))
with open(f_mj,'r') as f:
    for line in f:
        mj_params[amiAciToPot[line.split(' ')[0]], amiAciToPot[line.split(' ')[1]]] = float(line.split(' ')[2].split('\n')[0])
        mj_params[amiAciToPot[line.split(' ')[1]], amiAciToPot[line.split(' ')[0]]] = float(line.split(' ')[2].split('\n')[0])

# get the mj values for each protein
lst_mj_contacts = []
for num in range(len(f_templates)):
    contact_matrix = mj(mj_params,f_templates[num],lst_nresid[num])
    lst_mj_contacts.append(contact_matrix)
    print(num)

# get the mj value, numContact correlations for each protein
lst_freqs = []
lst_mj = []
lst_nativefreqs = []
lst_nativemj = []
lst_r = []
for num in range(len(f_templates)):
    freq,mj,native_freq,native_mj,r,p = numContacts_mj_corr(f_num_contacts[num],mj_params,f_templates[num],lst_nresid[num])
    lst_freqs.append(freq)
    lst_mj.append(mj)
    lst_nativefreqs.append(native_freq)
    lst_nativemj.append(native_mj)
    lst_r.append((r,p))
    print(num)

out_folder = '00_Figures/s2_mj_corr/'

plt.rcParams['xtick.labelsize']=40
plt.rcParams['ytick.labelsize']=40
# plot all numContacts FE plots
for i in range(9):
    numResid = lst_nresid[i]
    # Output contact potential distribution to a 2D discrete plot
    # create graph and colorbar
    fig,ax1 = plt.subplots(figsize=(11,9))
    pos = ax1.imshow(lst_numContacts[i],cmap='coolwarm_r',extent=(0,numResid*1000,0,numResid*1000))
    if i == 0:
        cbar = fig.colorbar(pos,ax=ax1,ticks = [0,-3,-6,-9,-12,-15],pad=0.08)
    else:
        cbar = fig.colorbar(pos,ax=ax1,pad=0.08)
    cbar.ax.set_xlabel(r"$\Delta$G" + "\n" + r"($k_{B}T$)",fontsize=36)
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
    # set ylabel parameters
    plt.ylabel('Residue',fontsize=40)
    ax1.yaxis.set_label_coords(-.17, .5)
    # set xlabel parameters
    plt.xlabel('Residue',fontsize=40)
    ax1.xaxis.set_label_coords(.5, -.13)
    # set title
    plt.title(lst_names[i],fontsize=40)
    plt.subplots_adjust(bottom = 0.16, top=0.91, left=0.18,right = 0.96,wspace=0.25,hspace=0.25)
    out_f = out_folder + fname[i] + '_FE.png'
    plt.savefig(out_f)
    plt.close()

# plot all MJ plots
for i in range(9):
    numResid = lst_nresid[i]
    # Output contact potential distribution to a 2D discrete plot
    # create graph and colorbar
    fig,ax1 = plt.subplots(figsize=(11,9))
    pos = ax1.imshow(lst_mj_contacts[i],cmap='coolwarm_r',extent=(0,numResid*1000,0,numResid*1000))
    cbar = fig.colorbar(pos,ax=ax1,pad=0.08)
    cbar.ax.set_xlabel("MJ",fontsize=36)
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
    # set ylabel parameters
    plt.ylabel('Residue',fontsize=40)
    ax1.yaxis.set_label_coords(-.17, .5)
    # set xlabel parameters
    plt.xlabel('Residue',fontsize=40)
    ax1.xaxis.set_label_coords(.5, -.13)
    # set title
    plt.title(lst_names[i],fontsize=40)
    plt.subplots_adjust(bottom = 0.16, top=0.91, left=0.18,right = 0.96,wspace=0.25,hspace=0.25)
    out_f = out_folder + fname[i] + '_MJ.png'
    plt.savefig(out_f)
    plt.close()

# plot all free energy mj correlation plots
for i in range(9):
# reset tick label size
    numResid = lst_nresid[i]
    fig,ax1 = plt.subplots(figsize=(13,9))
    plt.scatter(lst_freqs[i],lst_mj[i],c = 'blue',s=200)
    plt.scatter(lst_nativefreqs[i],lst_nativemj[i],c='red',s=200)
    plt.xlabel(r'MD Free Energy ($k_BT$)',fontsize=40)
    ax1.xaxis.set_label_coords(.555, -.12)
    plt.yticks(np.arange(0,-8,-2,dtype='int'))
    plt.ylabel('MJ Energy',fontsize=40)
    ax1.yaxis.set_label_coords(-.12, .5)
    plt.title(lst_names[i],fontsize=40)
    #Create legend
    red_patch = mpatches.Patch(color='red', label='Native')
    blue_patch = mpatches.Patch(color='blue', label= 'Nonnative')
    plt.legend(handles=[red_patch,blue_patch], bbox_to_anchor=[0.99, 0.01], loc='lower right', borderaxespad=0., fontsize=30)
    plt.subplots_adjust(bottom = 0.17, top=0.91, left=0.15,right = 0.96,wspace=0.25,hspace=0.25)
    out_f = out_folder + fname[i] + '_corr.png'
    plt.savefig(out_f)
    plt.close()
