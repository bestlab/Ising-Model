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

# function for sorting networks size
def getValue(e):
    return e['value']

def largest_subgraph(lower_bound,upper_bound,f_template,f_param,numResid):
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
    with open(f_param,'r') as f:
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
    with open(f_param,'r') as f:
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

def get_network_size(lower, upper, interval, f_template, f_param, numResid, write, plot):
    # loops through all DESRES trajectories
    lst_sigma = []
    lst_network_size = []
    outStr = ""
    for sigma in range(lower,upper,interval):
        lst_sigma.append(sigma/10)
        value = largest_subgraph(-100,sigma/10, f_template, f_param, numResid)
        lst_network_size.append(value)
        print((sigma/10,value))
        outStr += str((sigma/10,value)) + '\n'
    if write == True:
        out_f = "network_size.dat"
        with open(out_f,'w') as f:
            f.write(outStr)
    if plot == True:
        plt.rcParams['xtick.labelsize'] = 20
        plt.rcParams['ytick.labelsize'] = 20
        # create plot
        plt.figure(figsize=(14,10))
        plt.plot(lst_sigma,lst_network_size,linewidth=3,c='blue')
        plt.axhline(y=60,color='r',linewidth=2)
        plt.ylabel('Largest Connected Subgraph',fontsize=20)
        plt.xlabel('Threshold Cooperativity (Z)',fontsize=20)
        out_f = "network_size.png"
        plt.savefig(out_f)
        plt.close()
    return lst_sigma, lst_network_size

def calculate_threshold(lst_sigma,lst_network,expected_resid):
    critical_index = -1
    for i in range(len(lst_network)):
        if lst_network[i] == expected_resid:
            critical_index = i
        elif (lst_network[i] < expected_resid and lst_network[i-1] > expected_resid):
            critical_index = i
    if critical_index != -1:
        return lst_sigma[critical_index]
    else:
        return "Threshold cooperativity not found"

def get_coop_network(f_param,f_template,numResid,threshold,plot,plot_total):
    lower = -100
    upper = threshold
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
    with open(f_param,'r') as f:
        for line in f:
            if line[0] == 'J':
                value = float(line.split(' ')[3].split('\n')[0])
                lst_j.append(value)
    CUT2 = np.average(lst_j) + upper*np.std(lst_j)
    if lower == -100:
        CUT1 = min(lst_j) - 1
    else:
        CUT1 = np.average(lst_j) + lower*np.std(lst_j)
    # make list of high contacts
    lstCoop = []
    with open(f_param,'r') as f:
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
    if plot == True:
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
            fig,ax1 = plt.subplots(figsize=(20,20))
            pos = ax1.imshow(coopMap,cmap=cmap,norm=norm,extent=(0,numResid*1000,0,numResid*1000))
            # set ylabel parameters
            ax1.set_yticks(np.arange(numResid*1000-500,500,-4000,dtype='int'))
            ax1.set_yticklabels(np.arange(1,numResid+1,4,dtype='int'),fontsize=40)
            plt.ylabel('Residue',fontsize=40)
            # set xlabel parameters
            ax1.set_xticks(np.arange(500,numResid*1000+500,4000,dtype='int'))
            ax1.set_xticklabels(np.arange(1,numResid+1,4,dtype='int'),fontsize=40)
            plt.xlabel('Residue',fontsize=40)
            # plot all the connections
            for contact1 in adjList.keys():
                if contact1 in lstNetworks[network_size[i]['index']]:
                    for contact2 in adjList[contact1]:
                        x1 = 1000*name2pair[contact1][1]+500
                        x2 = 1000*name2pair[contact2][1]+500
                        y1 = (numResid*1000-500) - (1000*name2pair[contact1][0])
                        y2 = (numResid*1000-500) - (1000*name2pair[contact2][0])
                        plt.plot((x1,x2),(y1,y2),marker = 'o',color=edge_cmap[0],markersize=0.5,linewidth=0.6,alpha=0.2)
            # set title
            title = 'Cooperative Network {num}'.format(num = i)
            plt.title(title,fontsize=36)
            f_out = 'coopNetwork_{num}.png'.format(num = i)
            plt.savefig(f_out)
            plt.close()
    if plot_total == True:
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
        fig,ax1 = plt.subplots(figsize=(20,20))
        pos = ax1.imshow(coopMapTot,cmap=cmap,norm=norm,extent=(0,numResid*1000,0,numResid*1000))
        # set ylabel parameters
        ax1.set_yticks(np.arange(numResid*1000-500,500,-4000,dtype='int'))
        ax1.set_yticklabels(np.arange(1,numResid+1,4,dtype='int'),fontsize=40)
        plt.ylabel('Residue',fontsize=40)
        # set xlabel parameters
        ax1.set_xticks(np.arange(500,numResid*1000+500,4000,dtype='int'))
        ax1.set_xticklabels(np.arange(1,numResid+1,4,dtype='int'),fontsize=40)
        plt.xlabel('Residue',fontsize=40)
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
        title = 'Every Cooperative Network'
        plt.title(title,fontsize=36)
        f_out = 'coopNetworkAll.png'.format(num = i)
        plt.savefig(f_out)
        plt.close()
    return lstNetworks, adjList

# main function which calls up all other functions
if __name__ == "__main__":
    # User input parameters
    f_param = "ntl9_stable.bm_param"
    f_template = "NTL9_native.pdb"
    numResid = 39
    lower = -100
    upper = 1
    interval = 1
    expected_resid = 60
    # calculate size of networks over an interval
    # the lower, upper, intervals all get divided by a factor of 10
    sigma, network_size = get_network_size(lower,upper,interval,f_template,f_param,numResid,True,True)
    # calculate threshold cooperativity (Z)
    threshold = calculate_threshold(sigma,network_size,expected_resid)
    # calculate all contacts which are in network
    get_coop_network(f_param,f_template,numResid,threshold,True,True)
