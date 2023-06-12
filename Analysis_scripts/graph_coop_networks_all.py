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
def get_largest_coop_network(f_param,numResid,cutoff):
    lower = cutoff
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
    # conduct depth first search using list of cooperative contacts
    lstNetworks = []
    vTot = []
    replaceableAdj = copy(adjList)
    for contact in list(adjList):
        if not (contact in vTot):
            network = []
            network.append(contact)
            queue = []
            queue.append(contact)
            visited = []
            visited.append(contact)
            while queue:
                s = queue.pop(0)
                for coopContact in adjList[s]:
                    if not(coopContact in visited):
                        network.append(coopContact)
                        queue.append(coopContact)
                        visited.append(coopContact)
            lstNetworks.append(network)
            for i in range(len(visited)):
                vTot.append(visited[i])
    # get network size
    network_size = []
    for i in range(len(lstNetworks)):
        network_size.append({'index':i,'value':len(lstNetworks[i])})
    # sort networks by size
    network_size.sort(key=getValue,reverse=True)
    return lstNetworks, adjList

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
    'Ubiquitin/Combined_1UBQ/Ising_Model/ubiquitin_stable_cont.bm_param',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/Ising_Model/villin_stable.bm_param',
    'WW domain/Combined_GTT/Ising_Model/ww_domain_stable.bm_param']
    f_templates = [
    'A3D/Combined_A3D/Ising_Model/A3D_native.pdb',
    'BBA/Combined_1FME/Ising_Model/1FME_native.pdb',
    'l-repressor/Combined_lambda/Ising_Model/lambda_native.pdb',
    'Ntl9/Combined_NTL9/Ising_Model/NTL9_native.pdb',
    'Protein G/Combined_NuG2/Ising_Model/template_NuG2.pdb',
    'Trp-cage/DESRES-Trajectory_2JOF-0-protein/Contacts/Ising_Model/2JOF_native.pdb',
    'Ubiquitin/Combined_1UBQ/Ising_Model/template_1ubq.pdb',
    'Villin/DESRES-Trajectory_2F4K-0-protein/Contacts/Ising_Model/template_2F4K.pdb',
    'WW domain/Combined_GTT/Ising_Model/template_GTT.pdb']
    lst_nresid = [73,28,80,39,56,20,76,35,35]
    lst_coop_cutoff = [15.9,7.4,11.1,8.5,13.4,7.6,12.2,6.4,2.3]
    lst_names = [r'$\alpha$3D','BBA',r'$\lambda$-repressor','NTL9','NuG2','Trp-cage','Ubiquitin','Villin','WW Domain']
    # initialize lists
    lst_networks = []
    lst_adjLists = []
    for i in range(len(f_params)):
        network, adjList = get_largest_coop_network(f_params[i],lst_nresid[i],lst_coop_cutoff[i])
        lst_networks.append(network)
        lst_adjLists.append(adjList)
        print(i,len(network))
    # reset tick label size
    plt.rcParams['xtick.labelsize']=20
    plt.rcParams['ytick.labelsize']=20
    # create subplot of all the figures
    fig, ax = plt.subplots(3,3,figsize=(40,40))
    # graph each figure individually
    for num in range(9):
        numResid = lst_nresid[num]
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
        # Find native contacts and native contact names
        traj = md.load(f_templates[num])
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
        # Calculate number of networks
        num_network = len(lst_networks[num])
        edge_cmap = plt.cm.Purples_r(np.linspace(0,1,num_network))
        # graph potential on individual network
        coopMap = np.zeros((numResid,numResid))
        # Find potentials of residues in contact
        for contact in lst_networks[num]:
            resid1 = name2pair[contact][0]
            resid2 = name2pair[contact][1]
            if (resid1 < resid2):
                if contact in lstNameNative:
                    coopMap[resid1][resid2] = 2
                else:
                    coopMap[resid1][resid2] = 1
            else:
                if contact in lstNameNative:
                    coopMap[resid2][resid1] = 2
                else:
                    coopMap[resid2][resid1] = 1
        for contact in lstNameNative:
            resid2 = name2pair[contact][0]
            resid1 = name2pair[contact][1]
            coopMap[resid1][resid2]=2
        # Output contact potential distribution to a 2D discrete plot
        # create discrete colormap for labeling
        cmap = colors.ListedColormap(['white','#4169E1','#F32C2C'])
        bounds = [-0.5,0.5,1.5,2.5]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        # create graph and colorbar
        pos = ax[int(num/3)][num%3].imshow(coopMap,cmap=cmap,norm=norm,extent=(0,numResid*1000,0,numResid*1000))
        # set ylabel parameters
        ax[int(num/3)][num%3].set_yticks(np.arange(numResid*1000-500,500,-4000,dtype='int'))
        ax[int(num/3)][num%3].set_yticklabels(np.arange(1,numResid+1,4,dtype='int'),fontsize=28)
        # set xlabel parameters
        ax[int(num/3)][num%3].set_xticks(np.arange(500,numResid*1000+500,4000,dtype='int'))
        ax[int(num/3)][num%3].set_xticklabels(np.arange(1,numResid+1,4,dtype='int'),fontsize=28)
        ax[int(num/3)][num%3].set(xlabel='Residue',ylabel='Residue')
        # plot all the connections
        for contact1 in lst_adjLists[num].keys():
            if contact1 in lst_networks[num]:
                for contact2 in lst_adjLists[num][contact1]:
                    x1 = 1000*name2pair[contact1][1]+500
                    x2 = 1000*name2pair[contact2][1]+500
                    y1 = (numResid*1000-500) - (1000*name2pair[contact1][0])
                    y2 = (numResid*1000-500) - (1000*name2pair[contact2][0])
                    ax[int(num/3)][num%3].plot((x1,x2),(y1,y2),marker = 'o',color=edge_cmap[0],markersize=0.5,linewidth=0.6,alpha=0.1)
    # create legend
    fig.suptitle('Largest Cooperative Network',fontsize=24)
    out_f = "coop_networks.png"
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
lst_coop_cutoff = [7.6,7.4,6.4,10.6,8.5,13.5,15.9,17.8,11.0]
lst_alpha = [.7,.5,.2,.2,.1,.1,.1,.03,.05]
lst_names = ['Trp-cage','BBA','Villin','WW Domain','NTL9','NuG2',r'$\alpha$3D','Ubiquitin',r'$\lambda$-repressor']
# initialize lists
lst_networks = []
lst_adjLists = []
for i in range(len(f_params)):
    networks, adjList = get_largest_coop_network(f_params[i],lst_nresid[i],lst_coop_cutoff[i])
    lst_networks.append(networks)
    lst_adjLists.append(adjList)
    print(i,len(networks))

# reset tick label size
plt.rcParams['xtick.labelsize']=40
plt.rcParams['ytick.labelsize']=40
# create subplot of all the figures
fig, ax = plt.subplots(3,3,figsize=(25,28))
# graph each figure individually
for num in range(9):
    numResid = lst_nresid[num]
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
    # Find native contacts and native contact names
    traj = md.load(f_templates[num])
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
    # Calculate number of networks
    num_network = len(lst_networks[num])
    edge_cmap = plt.cm.Purples_r(np.linspace(0,1,num_network))
    # graph the total contact network
    coopMapTot = np.zeros((numResid,numResid))
    # iterate through the list of networks
    for i in range(0,len(lst_networks[num])):
        # Find potentials of residues in contact
        for contact in lst_networks[num][i]:
            resid1 = name2pair[contact][0]
            resid2 = name2pair[contact][1]
            if (resid1 < resid2):
                if contact in lstNameNative:
                    coopMapTot[resid1][resid2] = 2
                else:
                    coopMapTot[resid1][resid2] = 1
            else:
                if contact in lstNameNative:
                    coopMapTot[resid2][resid1] = 2
                else:
                    coopMapTot[resid2][resid1] = 1
        for contact in lstNameNative:
            resid2 = name2pair[contact][0]
            resid1 = name2pair[contact][1]
            coopMapTot[resid1][resid2]=2
    # Output contact potential distribution to a 2D discrete plot
    # create discrete colormap for labeling
    cmap = colors.ListedColormap(['white','#4169E1','#F32C2C'])
    bounds = [-0.5,0.5,1.5,2.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
    # create graph and colorbar
    pos = ax[int(num/3)][num%3].imshow(coopMapTot,cmap=cmap,norm=norm,extent=(0,numResid*1000,0,numResid*1000))
    # set ylabel parameters
    if numResid < 50:
        ax[int(num/3)][num%3].set_yticks(np.arange(numResid*1000-4500,499,-5000,dtype='int'))
        ax[int(num/3)][num%3].set_yticklabels(np.arange(5,numResid+1,5,dtype='int'),fontsize=40)
        # set xlabel parameters
        ax[int(num/3)][num%3].set_xticks(np.arange(4500,numResid*1000+500,5000,dtype='int'))
        ax[int(num/3)][num%3].set_xticklabels(np.arange(5,numResid+1,5,dtype='int'),fontsize=40)
    else:
        ax[int(num/3)][num%3].set_yticks(np.arange(numResid*1000-9500,499,-10000,dtype='int'))
        ax[int(num/3)][num%3].set_yticklabels(np.arange(10,numResid+1,10,dtype='int'),fontsize=40)
        # set xlabel parameters
        ax[int(num/3)][num%3].set_xticks(np.arange(9500,numResid*1000+500,10000,dtype='int'))
        ax[int(num/3)][num%3].set_xticklabels(np.arange(10,numResid+1,10,dtype='int'),fontsize=40)
    ax[int(num/3)][num%3].set_title(lst_names[num],fontsize=50)
    # plot all the connections by iterating through the list of networks
    for i in range(0,len(lst_networks[num])):
        for contact1 in lst_adjLists[num].keys():
            if contact1 in lst_networks[num][i]:
                for contact2 in lst_adjLists[num][contact1]:
                    x1 = 1000*name2pair[contact1][1]+500
                    x2 = 1000*name2pair[contact2][1]+500
                    y1 = (numResid*1000-500) - (1000*name2pair[contact1][0])
                    y2 = (numResid*1000-500) - (1000*name2pair[contact2][0])
                    ax[int(num/3)][num%3].plot((x1,x2),(y1,y2),marker = 'o',color=edge_cmap[0],markersize=0.5,linewidth=0.6,alpha=lst_alpha[num])

# create legend
red_patch = mpatches.Patch(color='#F32C2C', label='Native')
blue_patch = mpatches.Patch(color='#4169E1',label='Nonnative')
fig.legend(handles=[red_patch,blue_patch],loc='upper right',fontsize=40,borderaxespad=0.1,bbox_to_anchor=[0.97,1.0])
fig.supxlabel('Residue',x=0.535,fontsize=56)
fig.supylabel('Residue',fontsize=56)
plt.subplots_adjust(bottom = 0.07, top = 0.90, left=0.1,right = 0.97,wspace=0.25,hspace=0.25)
out_f = "coop_networks_all.png"
plt.savefig(out_f)
plt.close()
