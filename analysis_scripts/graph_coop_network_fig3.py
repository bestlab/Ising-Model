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

lower = 8.5
upper = 100
numResid = 39

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

f_in = 'ntl9_stable.bm_param'
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

f_in = 'ntl9_stable.bm_param'
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

# Calculate number of networks
num_network = len(lstNetworks)
edge_cmap = plt.cm.Purples_r(np.linspace(0,1,num_network))

# Find native contacts and native contact names
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

# graph the individual contact networks
# iterate through the list of networks
for i in range(0,len(lstNetworks)):
    # graph potential on individual network
    coopMap = np.zeros((numResid,numResid))
    # Find potentials of residues in contact
    for contact in lstNetworks[network_size[i]['index']]:
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
    fig,ax1 = plt.subplots(figsize=(10,10))
    pos = ax1.imshow(coopMap,cmap=cmap,norm=norm,extent=(0,numResid*1000,0,numResid*1000))
    # set ylabel parameters
    ax1.set_yticks(np.arange(numResid*1000-4500,499,-5000,dtype='int'))
    ax1.set_yticklabels(np.arange(5,numResid+1,5,dtype='int'),fontsize=36)
    ax1.yaxis.set_label_coords(-.12, .5)
    plt.ylabel('Residue',fontsize=36)
    # set xlabel parameters
    ax1.set_xticks(np.arange(4500,numResid*1000+500,5000,dtype='int'))
    ax1.set_xticklabels(np.arange(5,numResid+1,5,dtype='int'),fontsize=36)
    ax1.xaxis.set_label_coords(.50, -.10)
    plt.xlabel('Residue',fontsize=36)
    # plot all the connections
    for contact1 in adjList.keys():
        if contact1 in lstNetworks[network_size[i]['index']]:
            for contact2 in adjList[contact1]:
                x1 = 1000*name2pair[contact1][1]+500
                x2 = 1000*name2pair[contact2][1]+500
                y1 = (numResid*1000-500) - (1000*name2pair[contact1][0])
                y2 = (numResid*1000-500) - (1000*name2pair[contact2][0])
                plt.plot((x1,x2),(y1,y2),marker = 'o',color=edge_cmap[0],markersize=0.5,linewidth=0.6,alpha=0.2)
    plt.subplots_adjust(bottom = 0.115, top=0.97, left=0.17,right = 0.97,wspace=0.25,hspace=0.25)
    # set title
    #if i == 0:
        #title = 'NTL9 Largest Cooperative Network'.format(num = i)
        #plt.title(title,fontsize=32)
    #else:
        #title = 'NTL9 Cooperative Network {num}'.format(num = i)
        #plt.title(title,fontsize=32)
    f_out = 'NTL9_stable_coopNetwork_{num}_fig3.png'.format(num = i)
    plt.savefig(f_out)
    plt.close()

# graph the total contact network
coopMapTot = np.zeros((numResid,numResid))
# iterate through the list of networks
for i in range(0,len(lstNetworks)):
    # Find potentials of residues in contact
    for contact in lstNetworks[network_size[i]['index']]:
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
fig,ax1 = plt.subplots(figsize=(10,10))
pos = ax1.imshow(coopMapTot,cmap=cmap,norm=norm,extent=(0,numResid*1000,0,numResid*1000))
# set ylabel parameters
ax1.set_yticks(np.arange(numResid*1000-4500,499,-5000,dtype='int'))
ax1.set_yticklabels(np.arange(5,numResid+1,5,dtype='int'),fontsize=36)
ax1.yaxis.set_label_coords(-.12, .5)
plt.ylabel('Residue',fontsize=36)
# set xlabel parameters
ax1.set_xticks(np.arange(4500,numResid*1000+500,5000,dtype='int'))
ax1.set_xticklabels(np.arange(5,numResid+1,5,dtype='int'),fontsize=36)
ax1.xaxis.set_label_coords(.50, -.10)
plt.xlabel('Residue',fontsize=36)
plt.subplots_adjust(bottom = 0.115, top=0.97, left=0.17,right = 0.97,wspace=0.25,hspace=0.25)
# plot all the connections
for i in range(0,len(lstNetworks)):
    for contact1 in adjList.keys():
        if contact1 in lstNetworks[network_size[i]['index']]:
            for contact2 in adjList[contact1]:
                x1 = 1000*name2pair[contact1][1]+500
                x2 = 1000*name2pair[contact2][1]+500
                y1 = (numResid*1000-500) - (1000*name2pair[contact1][0])
                y2 = (numResid*1000-500) - (1000*name2pair[contact2][0])
                plt.plot((x1,x2),(y1,y2),marker = 'o',color=edge_cmap[0],markersize=0.5,linewidth=0.6,alpha=0.2)

# set title
#title = 'NTL9 All Cooperative Networks'.format(num = i)
#plt.title(title,fontsize=32)
f_out = 'NTL9_stable_coopNetworkAll_fig3.png'.format(num = i)
plt.savefig(f_out)
plt.close()
