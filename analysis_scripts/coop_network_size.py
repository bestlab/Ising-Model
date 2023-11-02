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

def largest_subgraph(lower_bound,upper_bound):
    ### CONSTANTS
    lower = lower_bound
    upper = upper_bound
    # Load all potential contacts
    # Make initial list of contacts and get file names
    lstPotContact = []
    for j in range(0,numResid):
        for k in range(j+1, numResid):
            if (abs(j-k)>3):
                lstPotContact.append([j,k])
    nameContacts = []
    for j in range(len(lstPotContact)):
        nameContacts.append(lstPotContact[j][0]*numResid + lstPotContact[j][1])
    name2pair = {}
    for j in range(len(lstPotContact)):
        name = lstPotContact[j][0]*numResid + lstPotContact[j][1]
        nameContacts.append(name)
        name2pair[name] = (lstPotContact[j][0],lstPotContact[j][1])
    # calculate the native contacts 
    traj = md.load('NTL9_native.pdb')
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
    f_in = 'ntl9_stable_freq.bm_param'
    with open(f_in,'r') as f:
        for line in f:
            if line[0] == 'J':
                value = float(line.split(' ')[3].split('\n')[0])
                lst_j.append(value)
    if lower == -100:
        CUT1 = min(lst_j) - 1
    else:
        CUT1 = np.average(lst_j) + lower*np.std(lst_j)
    CUT2 = np.average(lst_j) + upper*np.std(lst_j)
    # make list of high contacts
    lstCoop = []
    f_in = 'ntl9_stable_freq.bm_param'
    with open(f_in,'r') as f:
        for line in f:
            if line[0] == 'J':
                c1 = int(line.split(' ')[1])
                c2 = int(line.split(' ')[2])
                value = float(line.split(' ')[3].split('\n')[0])
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
    network_size = []
    for i in range(len(lstNetworks)):
        network_size.append(len(lstNetworks[i]))
    return max(network_size)

# main function which calls up all other functions
if __name__ == "__main__":
    # loops through all DESRES trajectories
    lst_sigma = []
    lst_network_size = []
    outStr = ""
    for sigma in range(-101,0,1):
        lst_sigma.append(sigma/10)
        value = largest_subgraph(-100,sigma/10)
        lst_network_size.append(value)
        print((sigma/10,value))
        outStr += str((sigma/10,value)) + '\n'
    out_f = "network_size_stable_Ntl9_freq.dat"
    with open(out_f,'w') as f:
        f.write(outStr)
    plt.rcParams['xtick.labelsize'] = 20
    plt.rcParams['ytick.labelsize'] = 20
    # create plot
    plt.figure(figsize=(14,10))
    plt.plot(lst_sigma,lst_network_size,linewidth=3,c='blue')
    plt.axhline(y=60,color='r',linewidth=2)
    plt.ylabel('Largest Connected Subgraph',fontsize=20)
    plt.xlabel('Standard Deviations above the Mean',fontsize=20)
    title = "NTL9"
    plt.title(title,fontsize=24)
    out_f = "network_size_stable_NTL9_freq.png"
    plt.savefig(out_f)
    plt.close()
