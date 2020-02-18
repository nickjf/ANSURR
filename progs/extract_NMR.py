import sys

pdb_in = open(sys.argv[1],'r')
pdb = sys.argv[1].split('.')[0]

allowed_res = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']

pdb_lines = []
model = '1' # assume first model is model 1, in case MODEL record missing, such as for 3trx
for line in pdb_in:
    if line[:6] == 'EXPDTA':
        if 'SOLUTION NMR' not in line:
            print('this version of ANSURR will only validate structures solved using SOLUTION NMR')
            quit()
    elif line[:5] == 'MODEL':
        model = line.split()[1]
    elif line[:4] == 'ATOM':
        resn = line[17:20]
        if resn not in allowed_res:
            print('this version of ANSURR does not suppport non-standard residues')
            quit()
        chain = line[21]
        if chain == ' ':
            chain = ''
        pdb_lines.append(line)
    elif 'TER' in line[0:3] or 'END' in line[0:3]:
        if len(pdb_lines) > 0:     
            out = open(pdb+chain+'_'+model+'.pdb','w')
            for l in pdb_lines:
                out.write(l)
            out.close()
            pdb_lines = []
        
