import sys

pdb_in = open(sys.argv[1],'r')
pdb = sys.argv[1].split('.')[0]
free_ligands = int(sys.argv[2])
nonstandard_res = int(sys.argv[3])

standard_res = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']

pdb_lines = []
model = '1' # assume first model is model 1
ter = 0
nonstandardres = []
freeligands = []
for line in pdb_in:
    if line[:6] == 'EXPDTA':
        if 'NMR' not in line:
            print('this version of ANSURR is only for validated structures solved using SOLUTION NMR')
            quit()
    elif line[:5] == 'MODEL':
        model = line.split()[1]
    elif line[:4] == 'ATOM' or line[:6] == 'HETATM':
        resn = line[17:20]
        chain = line[21]
        if chain == ' ':
            chain = ''
        if resn not in standard_res and nonstandard_res == 0 and ter == 0:
            if resn not in nonstandardres:
                print("\n -> found a non-standard residue ("+resn+") which are ignored by default. To include non-standard residues re-run ANSURR with the -n flag.",end='')
                nonstandardres.append(resn)
        if resn not in standard_res and free_ligands == 0 and ter == 1:
            if resn not in freeligands:
                print("\n -> found a free ligand ("+resn+") which are ignored by default. To include free ligands re-run ANSURR with the -l flag.",end='')
                freeligands.append(resn)
        if ter == 0 and nonstandard_res == 1:
            pdb_lines.append('ATOM  '+line[6:]) # this labels HETATMS which aren't free ligands as ATOMS so that FIRST recognises them as part of the protein   
        elif ter == 1 and free_ligands == 1:
            pdb_lines.append(line)
        elif resn in standard_res:
            pdb_lines.append(line)
    elif 'TER' in line[:3]:
        ter = 1
    elif 'END' in line[:3]:# or 'TER' in line[0:3]:
        if len(pdb_lines) > 0:     
            out = open(pdb+chain+'_'+model+'.pdb','w')
            for l in pdb_lines:
                out.write(l)
            out.close()
            pdb_lines = []
            ter = 0
if len(freeligands) + len(nonstandardres) > 0:
    print('')
        
