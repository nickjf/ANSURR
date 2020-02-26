# args pdb decomp_list

import sys, os
from collections import Counter

pdb = open(sys.argv[1],'r')
decomp = open(sys.argv[2],'r')

class res(object):
    _registry = []
    def __init__(self,i,name,atom):
        self._registry.append(self)
        self.i = i
        self.name = name
        self.atoms = [atom]
        self.clusters = []
        self.energy = 0

# get atom numbers for each CA atom from pdb
atom_number = 0
for line in pdb:
    if line[:4] == 'ATOM':
        atom_number += 1
        atom_name = line[12:16].replace(' ','')
        resn = line[17:20].replace(' ','')
        resi = int(line[22:26])
        if atom_name == 'CA':
            r = res(resi,resn,atom_number)

# import decomp_list         
data = []
all_data = []
energy = []
rigid_clusters = []
for line in decomp:
    if 'END' in line:
        rc = []
        data = data[:atom_number] # remove overflow
        all_data.append(data)
        cluster_count = Counter(data)
        for c in cluster_count:
            if int(cluster_count[c]) >= 15:
                rc.append(c)
        rigid_clusters.append(rc)
        data = []
    elif 'A:' in line:
        energy.append(float(line[10:21]))
    elif 'A' not in line and 'B' not in line:
        line = line.replace('\n','').split(':')
        for l in line:
            data.append(int(l))  

# append rigid clusters to residues
for d in all_data:
    for r in res._registry:
        r.clusters.append([d[a] for a in r.atoms])

# find rigid clusters
for x in range(len(all_data)):
    for r in res._registry:
        cluster = [y for y in r.clusters[x] if y != 0]
        for c in cluster:
            if c in rigid_clusters[x]:
                r.energy = energy[x]

#write output
out = open(os.path.basename(sys.argv[1]).split('.')[0]+'.decomp','w')
for r in res._registry:
	out.write(str(r.i) + ' ' + r.name +' ' + str(r.energy) + '\n')
out.close()
	
