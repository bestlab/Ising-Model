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
def coop_hist(f_param,f_template,numResid):
    # find the native contacts
    traj = md.load(f_template)
    prot = traj.atom_slice(traj.topology.select("protein"))
    numResid = numResid
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
    # make list of high contacts
    lstValues = []
    # sort contacts into native, native given nonnative, nonnative given native, and nonnative
    lstNative = []
    lstMixed = []
    lstNonnative = []
    with open(f_param,'r') as f:
        for line in f:
            if line[0] == 'J':
                c1 = int(line.split(' ')[1])
                c2 = int(line.split(' ')[2])
                value = float(line.split(' ')[5].split('\n')[0])
                lstValues.append(value)
                if nameContacts[c1] in lstNameNative and nameContacts[c2] in lstNameNative:
                    lstNative.append(value)
                elif nameContacts[c1] in lstNameNonnative and nameContacts[c2] in lstNameNonnative:
                    lstNonnative.append(value)
                else:
                    lstMixed.append(value)
    return lstNative,lstMixed,lstNonnative

# main function which calls up all other functions
if __name__ == "__main__":
    # initialize variables
    f_params = [
    'A3D/Combined_A3D/Ising_Model/a3d_stable.bm_param',
    'BBA/Combined_1FME/Ising_Model/bba_stable.bm_param',
    'l-repressor/Combined_lambda/Ising_Model/lambda_stable.bm_param',
    'Ntl9/Combined_NTL9/Ising_Model/ntl9_stable.bm_param',
    'Protein G/Combined_NuG2/Ising_Model/nug2_stable.bm_param',
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/Ising_Model/trp-cage_stable.bm_param',
    'Ubiquitin/Combined_1UBQ/Ising_Model_10ns_filter/ubiquitin_stable.bm_param',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/Ising_Model/villin_stable.bm_param',
    'WW domain/Combined_GTT/Ising_Model/ww_domain_stable.bm_param']
    f_templates = [
    'A3D/Combined_A3D/Ising_Model/A3D_native.pdb',
    'BBA/Combined_1FME/Ising_Model/1FME_native.pdb',
    'l-repressor/Combined_lambda/Ising_Model/lambda_native.pdb',
    'Ntl9/Combined_NTL9/Ising_Model/NTL9_native.pdb',
    'Protein G/Combined_NuG2/Ising_Model/template_NuG2.pdb',
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/Ising_Model/2JOF_native.pdb',
    'Ubiquitin/Combined_1UBQ/Ising_Model_10ns_filter/template_1ubq.pdb',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/Ising_Model/template_2F4K.pdb',
    'WW domain/Combined_GTT/Ising_Model/template_GTT.pdb']
    lst_nresid = [73,28,80,39,56,20,76,35,35]
    lst_names = [r'$\alpha$3D','BBA',r'$\lambda$-repressor','NTL9','NuG2','Trp-cage','Ubiquitin','Villin','WW Domain']
    # initialize lists
    lst_native_all = []
    lst_nonnative_all = []
    lst_mixed_all = []
    for i in range(len(f_params)):
        native,mixed,nonnative= coop_hist(f_params[i],f_templates[i],lst_nresid[i])
        lst_native_all.append(native)
        lst_mixed_all.append(mixed)
        lst_nonnative_all.append(nonnative)
        print(i)
    # reset tick label size
    plt.rcParams['xtick.labelsize']=20
    plt.rcParams['ytick.labelsize']=20
    # create subplot of all the figures
    fig, ax = plt.subplots(3,3)
    # graph each figure individually
    for i in range(len(ax)):
        native_max = np.max(lst_native_all[i])
        native_min = np.min(lst_native_all[i])
        mixed_max = np.max(lst_mixed_all[i])
        mixed_min = np.min(lst_mixed_all[i])
        nonnative_max = np.max(lst_nonnative_all[i])
        nonnative_min = np.min(lst_nonnative_all[i])
        ax[int(i/3)][i%3].hist(lst_native_all[i],bins=np.linspace(native_min,native_max,100),alpha=0.5,label='native',color='darkblue')
        ax[int(i/3)][i%3].hist(lst_nonnative_all[i],bins=np.linspace(nonnative_min,nonnative_max,100),alpha=0.5,label='nonnative',color='red')
        ax[int(i/3)][i%3].hist(lst_mixed_all[i],bins=np.linspace(mixed_min,mixed_max,100),alpha=0.5,label='mixed',color='yellow')
        ax.title = lst_names[i]
        ax.yscale("log")
    # create legend
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center')
    fig.supxlabel('J Value',fontsize=20)
    fig.title('Cooperativity Distributions',fontsize=24)
    out_f = "J_Hist_all_log.png"
    plt.savefig(out_f)
    plt.close()

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
lst_names = ['Trp-cage','BBA','Villin','WW Domain','NTL9','NuG2',r'$\alpha$3D','Ubiquitin',r'$\lambda$-repressor']
# initialize lists
lst_native_all = []
lst_nonnative_all = []
lst_mixed_all = []
for i in range(len(f_params)):
    native,mixed,nonnative= coop_hist(f_params[i],f_templates[i],lst_nresid[i])
    lst_native_all.append(native)
    lst_mixed_all.append(mixed)
    lst_nonnative_all.append(nonnative)
    print(i)

# reset tick label size
plt.rcParams['xtick.labelsize']=40
plt.rcParams['ytick.labelsize']=40
# create subplot of all the figures
fig, ax = plt.subplots(3,3,figsize=(30,30))
# graph each figure individually
for i in range(9):
    native_max = np.max(lst_native_all[i])
    native_min = np.min(lst_native_all[i])
    mixed_max = np.max(lst_mixed_all[i])
    mixed_min = np.min(lst_mixed_all[i])
    nonnative_max = np.max(lst_nonnative_all[i])
    nonnative_min = np.min(lst_nonnative_all[i])
    ax[int(i/3)][i%3].hist(lst_native_all[i],bins=np.linspace(native_min,native_max,100),alpha=0.5,label='native',color='darkblue')
    ax[int(i/3)][i%3].hist(lst_nonnative_all[i],bins=np.linspace(nonnative_min,nonnative_max,100),alpha=0.5,label='nonnative',color='red')
    ax[int(i/3)][i%3].hist(lst_mixed_all[i],bins=np.linspace(mixed_min,mixed_max,100),alpha=0.5,label='mixed',color='yellow')
    ax[int(i/3)][i%3].tick_params(labelsize=50)
    ax[int(i/3)][i%3].set_title(lst_names[i],fontsize=50)
    ax[int(i/3)][i%3].set_yscale("log")

# create legend
blue_patch = mpatches.Patch(color='darkblue', label='Native-Native')
red_patch = mpatches.Patch(color='red',label='Nonnative-Nonnative')
yellow_patch = mpatches.Patch(color='yellow',label='Native-Nonnative')
fig.legend(handles=[red_patch,yellow_patch,blue_patch],loc='upper right',fontsize=40,borderaxespad=0.10,bbox_to_anchor=[0.97,1.0])
fig.supxlabel(r'$\varepsilon_{ij}$ Value',x=0.52,fontsize=60)
plt.subplots_adjust(bottom = 0.08, top=0.87, left=0.07,right = 0.97,wspace=0.25,hspace=0.25)
out_f = "J_Hist_all_log_alt2.png"
plt.savefig(out_f)
plt.close()
