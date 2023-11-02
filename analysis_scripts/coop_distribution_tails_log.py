import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mdtraj as md

# for creating native contact map
from matplotlib import colors
import matplotlib.patches as mpatches

# inputs
numResid = 39
f_in = 'NTL9_native.pdb'
f_out = "J_dist_tails_NTL9_stable_log.png"

# find the native contacts
traj = md.load('NTL9_native.pdb')
prot = traj.atom_slice(traj.topology.select("protein"))

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

f_in = 'ntl9_stable.bm_param'
with open(f_in,'r') as f:
    for line in f:
    	if line[0] == 'J':
            c1 = int(line.split(' ')[1])
            c2 = int(line.split(' ')[2])
            value = float(line.split(' ')[3].split('\n')[0])
            lstValues.append(value)
            if nameContacts[c1] in lstNameNative and nameContacts[c2] in lstNameNative:
                lstNative.append(value)
            elif nameContacts[c1] in lstNameNonnative and nameContacts[c2] in lstNameNonnative:
                lstNonnative.append(value)
            else:
                lstMixed.append(value)

coop_max = np.max(lstValues)
coop_min = np.min(lstValues)
coop_avg = np.average(lstValues)
coop_std = np.std(lstValues)

# sort the lists
lstNative.sort()
lstMixed.sort()
lstNonnative.sort()

# compute the fraction of contacts that are above a certain fraction
e = np.linspace(coop_min,coop_avg-coop_std*1.95,1000)

native_tail = []
mixed_tail = []
nonnative_tail = []

for i in range(len(e)):
    native_count = 0
    nonnative_count = 0
    mixed_count = 0
    for j in range(len(lstNative)):
        if lstNative[j] < e[i]:
            native_count += 1
        else:
            break
    native_tail.append(native_count)
    for j in range(len(lstMixed)):
        if lstMixed[j] < e[i]:
            mixed_count += 1
        else:
            break
    mixed_tail.append(mixed_count)
    for j in range(len(lstNonnative)):
        if lstNonnative[j] < e[i]:
            nonnative_count += 1
        else:
            break
    nonnative_tail.append(nonnative_count)

# reset tick label size
plt.rcParams['xtick.labelsize']=40
plt.rcParams['ytick.labelsize']=40

# create histogram of density of all the contacts together
plt.figure(figsize=(14,10))
plt.plot(e,native_tail,linewidth=3,color='red')
plt.plot(e,nonnative_tail,linewidth=3,color='darkblue')
plt.plot(e,mixed_tail,linewidth=3,color='purple')
plt.yscale("log")
# create legend
blue_patch = mpatches.Patch(color='red', label='Native-Native')
red_patch = mpatches.Patch(color='darkblue',label='Nonnative-Nonnative')
yellow_patch = mpatches.Patch(color='purple',label='Native-Nonnative')
plt.legend(handles=[red_patch,yellow_patch,blue_patch],loc='upper right',fontsize=35,borderaxespad=0.10,bbox_to_anchor=[0.58,1.0])
plt.subplots_adjust(bottom = 0.13, top=0.95, left=0.10,right = 0.97,wspace=0.25,hspace=0.25)
plt.xlabel(r'$\varepsilon_{ij}$ Value',x=0.50,fontsize=40)
plt.savefig(f_out)
plt.close()
