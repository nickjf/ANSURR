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
        self.CA = atom
        self.clusters = []
        self.energy = 0

# get atom numbers for each CA atom from pdb
for line in pdb:
    if line[:4] == 'ATOM' or line[:6] == 'HETATM':
        atom_number = int(line[6:11])
        atom_name = line[12:16]
        resn = line[17:20].replace(' ','')
        resi = int(line[22:26])
        if atom_name == ' CA ': # this is occasionally a problem for ligands which have a CA atom, however, it's difficult to spot a ligand CA atom from a non-std residue CA atom, and doesn't effect scoring etc. Leave for now
            r = res(resi,resn,atom_number)

# import decomp_list         
data = []
all_data = []
energy = []
rigid_clusters = []
for line in decomp:
    if 'HEADER' not in line[:6]:
        if 'END' in line[:3]:
            rc = []
            data = data[:atom_number] # remove overflow
            all_data.append(data)
            cluster_count = Counter(data)
            for c in cluster_count:
                if int(cluster_count[c]) >= 15:
                    rc.append(c)
            rigid_clusters.append(rc)
            data = []
        elif 'A:' in line[:2]:
            energy.append(float(line[10:21]))
        elif 'A:' not in line[:2] and 'B:' not in line[:2]:
            line = line.replace('\n','').split(':')
            for l in line:
                data.append(int(l))  

# append rigid clusters to residues
for d in all_data:
     for r in res._registry:
        r.clusters.append(d[r.CA-1])

# find rigid clusters
for x in range(len(all_data)):
    for r in res._registry:
        if r.clusters[x] in rigid_clusters[x]:
            r.energy = energy[x]

#write output
out = open(os.path.basename(sys.argv[1]).split('.')[0]+'.decomp','w')
for r in res._registry:
	out.write(str(r.i) + ' ' + r.name +' ' + str(r.energy) + '\n')
out.close()
	
