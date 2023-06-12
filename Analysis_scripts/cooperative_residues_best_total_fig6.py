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
def count_aa(f_in,path,f_template):
    # initialize dictionary for counting residues
    aa_count = {
        'ALA': 0,
        'VAL': 0,
        'ILE': 0,
        'LEU': 0,
        'MET': 0,
        'PHE': 0,
        'TYR': 0,
        'TRP': 0,
        'HIS': 0,
        'PRO': 0,
        'ARG': 0,
        'LYS': 0,
        'ASP': 0,
        'GLU': 0,
        'SER': 0,
        'THR': 0,
        'ASN': 0,
        'GLN': 0,
        'GLY': 0,
    }
    aa_instances = {
        'ALA': 0,
        'VAL': 0,
        'ILE': 0,
        'LEU': 0,
        'MET': 0,
        'PHE': 0,
        'TYR': 0,
        'TRP': 0,
        'HIS': 0,
        'PRO': 0,
        'ARG': 0,
        'LYS': 0,
        'ASP': 0,
        'GLU': 0,
        'SER': 0,
        'THR': 0,
        'ASN': 0,
        'GLN': 0,
        'GLY': 0,
    }
    aa_index = {
        'ALA': 0,
        'VAL': 1,
        'ILE': 2,
        'LEU': 3,
        'MET': 4,
        'PHE': 5,
        'TYR': 6,
        'TRP': 7,
        'HIS': 8,
        'PRO': 9,
        'ARG': 10,
        'LYS': 11,
        'ASP': 12,
        'GLU': 13,
        'SER': 14,
        'THR': 15,
        'ASN': 16,
        'GLN': 17,
        'GLY': 18,
    }
    # open the specific sigma file for Ntl9
    numContacts = []
    lst_residues = []
    with open(f_in,'r') as f:
        for line in f:
            num = int(line.split('\t')[0].split(' ')[1])
            res1 = line.split('\t')[1].split(',')[0]
            res2 = line.split('\t')[1].split(',')[1]
            numContacts.append(num)
            lst_residues.append((res1,res2))
    # count all amino acids
    for i in range(len(numContacts)):
        res1 = lst_residues[i][0][0:3]
        res2 = lst_residues[i][1][0:3]
        if res1 in aa_index.keys() and res2 in aa_index.keys():
            aa_count[res1] += numContacts[i]
            aa_count[res2] += numContacts[i]
            aa_instances[res1] += 1
            aa_instances[res2] += 1
    lst_aa = []
    lst_aa_count = []
    for aa in aa_count.keys():
        lst_aa.append(aa)
        lst_aa_count.append(aa_count[aa])
    #Calculate contact potential
    aa_pair_count = np.zeros((19,19))
    aa_pair_instances = np.zeros((19,19))
    for i in range(len(numContacts)):
        res1 = lst_residues[i][0][0:3]
        res2 = lst_residues[i][1][0:3]
        if res1 in aa_index.keys() and res2 in aa_index.keys():
            index1 = aa_index[res1]
            index2 = aa_index[res2]
            aa_pair_count[index1][index2] += numContacts[i]
            aa_pair_count[index2][index1] += numContacts[i]
            aa_pair_instances[index1][index2] += 1
            aa_pair_instances[index2][index1] += 1
    return aa_count,aa_instances,aa_pair_count,aa_pair_instances

# main function which calls up all other functions
if __name__ == "__main__":
    # initialize variables
    f_ins = [
    'A3D/Combined_A3D/Ising_Model/J_network_vertices/A3D_vertices_list_1sigma.dat',
    'BBA/Combined_1FME/Ising_Model/J_network_vertices/BBA_vertices_list_1sigma.dat',
    'l-repressor/Combined_lambda/Ising_Model/J_network_vertices/lambda_vertices_list_1sigma.dat',
    'Ntl9/Combined_NTL9/Ising_Model/J_network_vertices/NTL9_vertices_list_1sigma.dat',
    'Protein G/Combined_NuG2/Ising_Model/J_network_vertices/NuG2_vertices_list_1sigma.dat',
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/Ising_Model/J_network_vertices/2JOF_vertices_list_1sigma.dat',
    'Ubiquitin/Combined_1UBQ/Ising_Model/J_network_vertices/ubq_vertices_list_1sigma.dat',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/Ising_Model/J_network_vertices/villin_vertices_list_1sigma.dat',
    'WW domain/Combined_GTT/Ising_Model/J_network_vertices/ww_domain_vertices_list_1sigma.dat']
    paths = [
    'A3D/Combined_A3D/Ising_Model/J_network_vertices/',
    'BBA/Combined_1FME/Ising_Model/J_network_vertices/',
    'l-repressor/Combined_lambda/Ising_Model/J_network_vertices/',
    'Ntl9/Combined_NTL9/Ising_Model/J_network_vertices/',
    'Protein G/Combined_NuG2/Ising_Model/J_network_vertices/',
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/Ising_Model/J_network_vertices/',
    'Ubiquitin/Combined_1UBQ/Ising_Model/J_network_vertices/',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/Ising_Model/J_network_vertices/',
    'WW domain/Combined_GTT/Ising_Model/J_network_vertices/']
    f_templates = [
    'A3D/Combined_A3D/A3D_native.pdb',
    'BBA/Combined_1FME/1FME_native.pdb',
    'l-repressor/Combined_lambda/lambda_native.pdb',
    'Ntl9/Combined_NTL9/NTL9_native.pdb',
    'Protein G/Combined_NuG2/template_NuG2.pdb',
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/2JOF_native.pdb',
    'Ubiquitin/Combined_1UBQ/template_1ubq.pdb',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/template_2F4K.pdb',
    'WW domain/Combined_GTT/template_GTT.pdb']
    aa_count_total = {
        'ALA': 0.0,
        'VAL': 0.0,
        'ILE': 0.0,
        'LEU': 0.0,
        'MET': 0.0,
        'PHE': 0.0,
        'TYR': 0.0,
        'TRP': 0.0,
        'HIS': 0.0,
        'PRO': 0.0,
        'ARG': 0.0,
        'LYS': 0.0,
        'ASP': 0.0,
        'GLU': 0.0,
        'SER': 0.0,
        'THR': 0.0,
        'ASN': 0.0,
        'GLN': 0.0,
        'GLY': 0.0,
    }
    aa_instances_total = {
        'ALA': 0.0,
        'VAL': 0.0,
        'ILE': 0.0,
        'LEU': 0.0,
        'MET': 0.0,
        'PHE': 0.0,
        'TYR': 0.0,
        'TRP': 0.0,
        'HIS': 0.0,
        'PRO': 0.0,
        'ARG': 0.0,
        'LYS': 0.0,
        'ASP': 0.0,
        'GLU': 0.0,
        'SER': 0.0,
        'THR': 0.0,
        'ASN': 0.0,
        'GLN': 0.0,
        'GLY': 0.0,
    }
    aa_avg = {
        'ALA': 0.0,
        'VAL': 0.0,
        'ILE': 0.0,
        'LEU': 0.0,
        'MET': 0.0,
        'PHE': 0.0,
        'TYR': 0.0,
        'TRP': 0.0,
        'HIS': 0.0,
        'PRO': 0.0,
        'ARG': 0.0,
        'LYS': 0.0,
        'ASP': 0.0,
        'GLU': 0.0,
        'SER': 0.0,
        'THR': 0.0,
        'ASN': 0.0,
        'GLN': 0.0,
        'GLY': 0.0,
    }
    #Calculate contact potential
    aa_pair_count_total = np.zeros((19,19))
    aa_pair_instances_total = np.zeros((19,19))
    # loops through all 9 proteins
    for i in range(len(f_ins)):
        count_dict, count_instances_dict, count_matrix, count_instances_matrix = count_aa(f_ins[i],paths[i],f_templates[i])
        for key in count_dict.keys():
            aa_count_total[key] += count_dict[key]
            aa_instances_total[key] += count_instances_dict[key]
        dim1 = count_matrix.shape[0]
        dim2 = count_matrix.shape[1]
        for j in range(dim1):
            for k in range(dim2):
                aa_pair_count_total[j][k] += count_matrix[j][k]
                aa_pair_instances_total[j][k] += count_instances_matrix[j][k]
        print(np.max(aa_pair_count_total))
    # create bar graph of total count
    lst_aa = []
    lst_aa_count_total = []
    for aa in aa_count_total.keys():
        lst_aa.append(aa)
        lst_aa_count_total.append(aa_count_total[aa])
    # create bar graph of total count normalized
    for aa in aa_avg.keys():
        if aa_instances_total[aa] != 0:
            aa_avg[aa] = float(aa_count_total[aa]/aa_instances_total[aa])
    lst_aa_avg = []
    for aa in aa_avg.keys():
        lst_aa_avg.append(aa_avg[aa])
    # Graph total contact pairs
    # normalize total contact pairs
    for i in range(19):
        for j in range(19):
            if aa_pair_instances_total[i][j] <= 2:
                aa_pair_count_total[i][j] = 0
            else:
                value = aa_pair_count_total[i][j]/aa_pair_instances_total[i][j]
                aa_pair_count_total[i][j] = value
    # create mask for contact potential
    mask = np.zeros((19,19))
    for i in range(19):
        for j in range(19):
            if (aa_pair_count_total[i][j] == 0):
                mask[i][j] = int(1)
    values_fe = np.ma.masked_array(aa_pair_count_total,mask)
    max_count = np.max(aa_pair_count_total)
    min_count = np.min(aa_pair_count_total)
    bounds = np.linspace(min_count,max_count,14)
    fig,ax1 = plt.subplots(figsize=(19,19))
    pos = ax1.imshow(values_fe,cmap='Blues',extent=(0,19000,0,19000))
    cbar = fig.colorbar(pos,ax=ax1)
    cbar.ax.set_xlabel("Avg." + '\n' + "Num. Coop." +  '\n' + "Interactions",fontsize=32)
    cbar.ax.tick_params(labelsize = 32)
    # plot hydrophobic boundaries
    #ax1.axhline(y=10000,xmin=0/19000,xmax=9000/19000,color='r',linewidth=6)
    #ax1.axvline(x=9000,ymin=10000/19000,ymax=1,color='r',linewidth=6)
    # set ylabel parameters
    ax1.set_yticks(np.arange(19*1000-500,0,-1000,dtype='int'))
    ax1.set_yticklabels(lst_aa,fontsize=36)
    ax1.yaxis.set_label_coords(-.14, .5)
    plt.ylabel('Amino Acid Type',fontsize=40)
    # set xlabel parameters
    ax1.set_xticks(np.arange(500,19*1000+500,1000,dtype='int'))
    ax1.set_xticklabels(lst_aa,fontsize=36)
    ax1.xaxis.set_label_coords(.5, -.14)
    plt.xticks(rotation=90)
    plt.xlabel('Amino Acid Type',fontsize=40)
    #plt.title("Cooperative Amino Acid Pairs",fontsize=44)
    plt.subplots_adjust(bottom = 0.16, top=0.91, left=0.18,right = 0.96,wspace=0.25,hspace=0.25)
    plt.savefig("aa_pair_avg_best_total.png")
    plt.close()
