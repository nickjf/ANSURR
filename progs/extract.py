import sys,pathlib,json

def make_monomers(pdb,chains):
	model = pdb.split('/')[1].split(chains)[0]
	num = pdb.split('/')[1].split(chains)[1].replace('_','').replace('.pdb','')
	prev_chain = 'XXX'
	for line in open(pdb,'r'):
		chain = line[21]
		if prev_chain == 'XXX':
			prev_chain = chain
			out = open(model+chain+'_'+num+'.pdb','w')
			out.write(line)
		elif chain != prev_chain:
			out.close()
			out = open(model+chain+'_'+num+'.pdb','w')
			out.write(line)
		else:
			out.write(line)
		prev_chain = chain
	out.close()
	
pdb_in = open(sys.argv[1],'r')
pdb = sys.argv[1].split('.')[0]
freeligands = int(sys.argv[2])
nonstandardres = int(sys.argv[3])
if sys.argv[4] == '1':
	oligomers = 1
else:
	oligomers = 0

standard_res = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']

pdb_lines = []
models = {}
model = '1' # assume first model is model 1
chain = ''
prev_chain = ''
nonstandard = []

for line in pdb_in:
	if line[:6] == 'EXPDTA':
		if 'NMR' not in line and 'SOLID' not in line: # some older solution NMR structures just reffered to as NMR
			print(' -> CRITICAL ERROR this version of ANSURR is only for validating structures solved using SOLUTION NMR, exiting')
			quit()
	elif line[:5] == 'MODEL':
		model = line.split()[1]
	elif line[:4] == 'ATOM' or line[:6] == 'HETATM':
		pdb_lines.append(line)
		chain = '' if chain == ' ' else line[21]
	elif 'TER' in line[:3] or 'END' in line[:3]:
		if len(pdb_lines) > 0:   
			selected_pdb_lines = []
			std_res = 0
			for line in pdb_lines:
				resn = line[17:20]
				if resn in standard_res:
					std_res = 1
					break
			if std_res == 1:
				for line in pdb_lines:
					resn = line[17:20]
					if resn in standard_res:
						selected_pdb_lines.append(line)
					elif nonstandardres == 1:
						selected_pdb_lines.append('ATOM  '+line[6:])
						if resn not in nonstandard:
							print(" -> non-standard residue "+resn.replace(' ','')+" will be included in rigidity calculations")
							nonstandard.append(resn)
					else:
						if resn not in nonstandard:
							print(" -> found a non-standard residue ("+resn.replace(' ','')+") which are ignored by default. To include non-standard residues, re-run ANSURR with the -n flag")
							nonstandard.append(resn)
			else:
				for line in pdb_lines:
					resn = line[17:20]
					if freeligands == 1:
						selected_pdb_lines.append(line)
						if resn not in nonstandard:
							print(" -> free ligand "+resn.replace(' ','')+" will be included in rigidity calculations")
							nonstandard.append(resn)
					else:
						if resn not in nonstandard:
							print(" -> found a free ligand ("+resn.replace(' ','')+") which are ignored by default. To include free ligands, re-run ANSURR with the -l flag")
							nonstandard.append(resn)
			
			if len(selected_pdb_lines) > 0:
				if model not in models:  
					models[model] = {chain:selected_pdb_lines}
				elif chain not in models[model]:
					models[model][chain] = selected_pdb_lines
				else:
					models[model][chain].extend(selected_pdb_lines)
			pdb_lines = []
	prev_chain = chain

chains_done = []
resi_ref = {}
for model in models:
	chains = ''
	if oligomers == 0:
		for chain in models[model]:
			chains += chain
			out = open(pdb+chain+'_'+model+'.pdb','w')
			for l in models[model][chain]:
				out.write(l)
			out.close()
		if len(chains) > 1 and chains not in chains_done:
			print(" -> found "+str(len(chains))+" chains, which are extracted seperately by default. To combine chains when calculating flexibility, re-run ANSURR with the -o flag")
			chains_done.append(chains)
	elif oligomers == 1:
		pathlib.Path("combined").mkdir(parents=True, exist_ok=True)
		for chain in models[model]:
			chains += chain
		out = open('combined/'+pdb+chains+'_'+model+'.pdb','w')
		count = 1
		num = 0
		prev_resi = -99999
		resi_ref[model] = {}
		for chain in models[model]:
			resi_ref[model][chain] = {'orig_first':int(models[model][chain][0][23:26]),'orig_last':int(models[model][chain][-1][23:26]),'new_first':'','new_last':''}
			for l in models[model][chain]:
				resi = int(l[23:26])
				if resi != prev_resi:
					num += 1
				prev_resi = resi
				if resi_ref[model][chain]['new_first'] == '':
					resi_ref[model][chain]['new_first'] = num
				out.write(l[0:6]+format(str(count)," >5s")+l[11:23]+''.join([' ']*(3-len(str(num))))+str(num)+l[26:])
				count +=1 
			resi_ref[model][chain]['new_last'] = num
		out.close()
		make_monomers('combined/'+pdb+chains+'_'+model+'.pdb',chains)
		if chains not in chains_done:
			print(" -> chains "+chains+" combined into a single structure to calculate flexibility ")	
			chains_done.append(chains)

if oligomers == 1:
	json.dump(resi_ref, open("resi_ref.txt",'w'))



