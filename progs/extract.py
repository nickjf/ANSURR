import sys,pathlib,json

def make_monomers(pdbid,chains,model):
	pdb = 'combined/'+pdbid+chains+'_'+model+'.pdb'
	prev_chain = 'XXX'
	for line in open(pdb,'r'):
		if 'ATOM' in line[:4] or 'HETATM' in line[:6]:
			chain = line[21]
			resi = int(line[22:26])
			orig_resi = resi - resi_ref[model][chain]['new_first'] + resi_ref[model][chain]['orig_first']
			if prev_chain == 'XXX':
				prev_chain = chain
				out = open(pdbid+chain+'_'+model+'.pdb','w')
			elif chain != prev_chain:
				out.close()
				out = open(pdbid+chain+'_'+model+'.pdb','w')
			out.write(line[0:22]+''.join([' ']*(4-len(str(orig_resi))))+str(orig_resi)+line[26:])
			prev_chain = chain
	out.close()

def parse_pdb_lines(pdb_lines):
	selected_pdb_lines = {}
	for chain in pdb_lines:
		std_res = 0
		for line in pdb_lines[chain]:
			resn = line[17:20]
			if resn in standard_res and line[:4] == 'ATOM': # helps diagnose non-standard residues from free ligands, non-standard residues are HETATMS that come with standard residues too
				std_res = 1
				break
		if std_res == 1:
			for line in pdb_lines[chain]:
				resn = line[17:20]
				if resn in standard_res:
					if chain not in selected_pdb_lines:
						selected_pdb_lines[chain] = [line]
					else:
						selected_pdb_lines[chain].append(line)
				elif nonstandardres == 1:
					if chain not in selected_pdb_lines:
						selected_pdb_lines[chain] = ['ATOM  '+line[6:]]
					else:
						selected_pdb_lines[chain].append('ATOM  '+line[6:])
					if resn not in nonstandard:
						print(" -> non-standard residue "+resn.replace(' ','')+" will be included in rigidity calculations")
						nonstandard.append(resn)
				else:
					if resn not in nonstandard:
						print(" -> found a non-standard residue ("+resn.replace(' ','')+") which are ignored by default. To include non-standard residues, re-run ANSURR with the -n flag")
						nonstandard.append(resn)
		else:
			for line in pdb_lines[chain]:
				resn = line[17:20]
				if freeligands == 1:
					if chain not in selected_pdb_lines:
						selected_pdb_lines[chain] = [line]
					else:
						selected_pdb_lines[chain].append(line)
					if resn not in nonstandard:
						print(" -> free ligand "+resn.replace(' ','')+" will be included in rigidity calculations")
						nonstandard.append(resn)
				else:
					if resn not in nonstandard:
						print(" -> found a free ligand ("+resn.replace(' ','')+") which are ignored by default. To include free ligands, re-run ANSURR with the -l flag")
						nonstandard.append(resn)
	for chain in selected_pdb_lines:
		if len(selected_pdb_lines[chain]) > 0:
			if model not in models:  
				models[model] = {chain:selected_pdb_lines[chain]}
			elif chain not in models[model]:
				models[model][chain] = selected_pdb_lines[chain]
			else:
				models[model][chain].extend(selected_pdb_lines[chain])
	return models

	
pdb_in = open(sys.argv[1],'r')
pdb = sys.argv[1].split('.')[0]
freeligands = int(sys.argv[2])
nonstandardres = int(sys.argv[3])
if sys.argv[4] == '1':
	oligomers = 1
else:
	oligomers = 0

standard_res = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']

pdb_lines = {}
models = {}
model = '1' # assume first model is model 1
chain = ''
nonstandard = []
conect = []

for line in pdb_in:
	if line[:6] == 'EXPDTA':
		if 'NMR' not in line and 'SOLID' not in line: # some older solution NMR structures just reffered to as NMR
			print(' -> CRITICAL ERROR this version of ANSURR is only for validating structures solved using SOLUTION NMR, exiting')
			quit()
	elif line[:5] == 'MODEL':
		model = line.split()[1]
	elif line[:4] == 'ATOM' or line[:6] == 'HETATM':
		chain = '' if line[21] == ' ' else line[21]
		if chain not in pdb_lines:
			pdb_lines[chain] = [line]
		else:
			pdb_lines[chain].append(line)
	elif 'TER' in line[:3] or 'END' in line[:3]: # parse after TER to help diagnose free ligands from non-sandard residues
		if len(pdb_lines) > 0:   
			models = parse_pdb_lines(pdb_lines)
			pdb_lines = {}
	elif 'CONECT'in line[:6]:
		conect.append(line)

if len(pdb_lines) > 0:   # catch any extra stuff in case structure doesn't end with TER or END, e.g. for CNS output
	models = parse_pdb_lines(pdb_lines)
	pdb_lines = {}

for model in models: # this orders by resi number as sometimes pdbs are not ordered correctly (e.g. pdb 2lrl)
	for chain in models[model]:
		pdb_lines = models[model][chain]
		resi_pdb_lines = {}

		for line in pdb_lines:
		    resi = int(line[22:26])
		    if resi not in resi_pdb_lines:
		        resi_pdb_lines[resi] = [line]
		    else:
		        resi_pdb_lines[resi].append(line)
		sorted_resi_pdb_lines = {}
		for key in sorted(resi_pdb_lines.keys()) :
		    sorted_resi_pdb_lines[key] = resi_pdb_lines[key]    

		pdb_lines = [sorted_resi_pdb_lines[i] for i in sorted_resi_pdb_lines]
		models[model][chain] =  [item for sublist in pdb_lines for item in sublist]

chains_done = []
resi_ref = {}
old_new_atom_num = {} # 
for model in models:
	chains = ''
	if oligomers == 1:
		if len(models[model]) > 1:
			pathlib.Path("combined").mkdir(parents=True, exist_ok=True)
			for chain in models[model]:
				chains += chain
			out = open('combined/'+pdb+chains+'_'+model+'.pdb','w')
			count = 1
			num = 0
			resi_ref[model] = {}
			for chain in models[model]:
				prev_resi = -99999
				resi_ref[model][chain] = {'orig_first':int(models[model][chain][0][22:26]),'orig_last':int(models[model][chain][-1][22:26]),'new_first':'','new_last':''}
				for l in models[model][chain]:
					old_new_atom_num[l[6:11].strip()] = str(count)
					resi = int(l[22:26])
					if prev_resi != -99999:
						num += resi - prev_resi
					else:
						num += 1
					prev_resi = resi
					if resi_ref[model][chain]['new_first'] == '':
						resi_ref[model][chain]['new_first'] = num
					out.write(l[0:6]+format(str(count)," >5s")+l[11:22]+''.join([' ']*(4-len(str(num))))+str(num)+l[26:].replace('\n','')+'\n') #this makes sure there is a single new line char at end (needed for FIRST to run correctly)
					count +=1 
				resi_ref[model][chain]['new_last'] = num
			out.write('END\n')
			for c in conect:
				print(c)
				corrected_c = ''
				for i in c.split()[1:]:
					if i in old_new_atom_num:
						corrected_c += ' '+old_new_atom_num[i]
				if len(corrected_c.split()) > 1:
					out.write('CONECT '+corrected_c+'\n')
			out.close()
			make_monomers(pdb,chains,model)
			if chains not in chains_done:
				print(" -> chains "+chains+" combined into a single structure to calculate flexibility ")	
				chains_done.append(chains)
		else:
			oligomers = 0
			print(' -> found only a single chain ('+str([i for i in models[model]][0])+'), no need to combine')
	if oligomers == 0:
		for chain in models[model]:
			chains += chain
			out = open(pdb+chain+'_'+model+'.pdb','w')
			count = 1
			for l in models[model][chain]:
				old_new_atom_num[l[6:11].strip()] = str(count)
				out.write(l[0:6]+format(str(count)," >5s")+l[11:])
				count +=1 
			out.write('END\n')
			for c in conect:
				corrected_c = ''
				for i in c.split()[1:]:
					if i in old_new_atom_num:
						corrected_c += ' '+old_new_atom_num[i]
				if len(corrected_c.split()) > 1:
					out.write('CONECT '+corrected_c+'\n')
			out.close()
		if len(chains) > 1 and chains not in chains_done:
			print(" -> found "+str(len(chains))+" chains, which are extracted seperately by default. To combine chains when calculating flexibility, re-run ANSURR with the -o flag")
			chains_done.append(chains)
if len(conect) > 0:
	print(" -> found CONECT records")
if oligomers == 1:
	json.dump(resi_ref, open("resi_ref.tmp",'w'))



