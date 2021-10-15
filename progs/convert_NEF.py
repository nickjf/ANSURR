from RCI_dicts import amino_acids, amino_acids_three_letter
import numpy as np
import sys, os

# this one uses loops, rather than save frames. it works, but prob should use save frames
file_basename = os.path.basename(sys.argv[1]).split('.')[0]

with open(sys.argv[1]) as f_in:
    lines = list(filter(None, (line.rstrip() for line in f_in)))   # ignore blank lines
    
fields = ['_nef_sequence','_nef_chemical_shift']
backbone_atoms = ['C','CA','CB','H','HA','N','QA','HA1','HA2','HA3','HAx','HAy']

chain_sequence = {}
chain_shifts = {}

loop = 0
loop_ID = 0   # this is useful in case a file has multiple sets of chemical shifts 
for line in lines:
    if not line.strip().startswith("#"):   # ignore commented out lines
        line = line.strip()
        if loop == 0:            # extract fields from each loop
            if line == 'loop_':
                f = []
                loop = 1
        elif loop == 1:
            if line.split('.')[0] in fields:
                field = line.split('.')[0]
                f.append(line)
                
            elif line == 'stop_':
                loop = 0
                loop_ID +=1
            elif len(f) > 0:
                line = line.split()
                if field == '_nef_sequence':
                    chain_ID = line[f.index('_nef_sequence.chain_code')]
                    sequence_ID = line[f.index('_nef_sequence.sequence_code')]
                    residue_name = line[f.index('_nef_sequence.residue_name')]
                    
                    if sequence_ID[-1].isnumeric(): # ignore multiples for sequence, just take the first
                        if chain_ID not in chain_sequence:
                            chain_sequence[chain_ID] = {sequence_ID:residue_name}
                        else:
                            chain_sequence[chain_ID][sequence_ID] = residue_name
                        
                elif field == '_nef_chemical_shift':
                    
                    chain_ID = line[f.index('_nef_chemical_shift.chain_code')]
                    sequence_ID = line[f.index('_nef_chemical_shift.sequence_code')]
                    atom_name = line[f.index('_nef_chemical_shift.atom_name')]
                    residue_name = line[f.index('_nef_chemical_shift.residue_name')]
                    shift_value = line[f.index('_nef_chemical_shift.value')]
                    
                    if atom_name in backbone_atoms and residue_name in amino_acids_three_letter: # just take backbone atoms from standard amino acids
                        if sequence_ID[-1].isnumeric(): # ignore multiples for sequence, just take the first
                            if loop_ID not in chain_shifts:
                                chain_shifts[loop_ID] = {}
                                chain_shifts[loop_ID][chain_ID] = {sequence_ID:{atom_name:shift_value}}
                            elif chain_ID not in chain_shifts[loop_ID]:
                                chain_shifts[loop_ID][chain_ID] = {sequence_ID:{atom_name:shift_value}}
                            elif sequence_ID not in chain_shifts[loop_ID][chain_ID]:
                                chain_shifts[loop_ID][chain_ID][sequence_ID] = {atom_name:shift_value}
                            elif atom_name in ['QA','HA1','HA2','HA3','HAx','HAy']: # average HAs i.e. for GLY
                                if 'HAs' not in chain_shifts[loop_ID][chain_ID][sequence_ID]:
                                    chain_shifts[loop_ID][chain_ID][sequence_ID]['HAs'] = [float(shift_value)]
                                else:
                                    chain_shifts[loop_ID][chain_ID][sequence_ID]['HAs'].append(float(shift_value))    
                            else:
                                chain_shifts[loop_ID][chain_ID][sequence_ID][atom_name] = shift_value


for loop_ID in chain_shifts:
    for chain_ID in chain_shifts[loop_ID]:
        for sequence_ID in chain_shifts[loop_ID][chain_ID]:
            if 'HAs' in chain_shifts[loop_ID][chain_ID][sequence_ID]:
                chain_shifts[loop_ID][chain_ID][sequence_ID]['HA'] = str(round(np.mean(chain_shifts[loop_ID][chain_ID][sequence_ID]['HAs']),2))
                chain_shifts[loop_ID][chain_ID][sequence_ID].pop('HAs', None)


shift_sets = {}  # for counting the number of sets of shifts and chains
i=1
for shift_id in chain_shifts:
    for chain in chain_shifts[shift_id]:
        lines_seq = []
        lines_shifts = []
        num_seq=1
        num_shifts=1
        seq_single_letter = ''
        old_sequence_ID = ''
        offset = ''
        if chain in chain_sequence:
            for sequence_ID in chain_shifts[shift_id][chain]:
                if offset == '':
                    offset = int(sequence_ID) - 1
                for shift in chain_shifts[shift_id][chain][sequence_ID]:
                    lines_shifts.append(str(num_shifts)+' '+str(int(sequence_ID)-offset)+' '+chain_sequence[chain][sequence_ID]+' '+shift+' '+chain_shifts[shift_id][chain][sequence_ID][shift] + ' 0 ' + str(int(sequence_ID)-offset))
                    num_shifts+=1
                x=1
                
                if old_sequence_ID == '':
                    old_sequence_ID = int(sequence_ID) - 1

                while True:
                    if int(sequence_ID) == old_sequence_ID + x:
                        seq_single_letter += amino_acids_three_letter[chain_sequence[chain][sequence_ID]]
                        lines_seq.append('  '+str(num_seq)+' . '+chain_sequence[chain][sequence_ID]+' . 99999 1')
                        num_seq+=1
                        break
                    else:
                        seq_single_letter += 'X'
                        lines_seq.append('  '+str(num_seq)+' . XXX . 99999 1')
                        x+=1
                        num_seq+=1

                        
                old_sequence_ID = old_sequence_ID + x
        
        if len(lines_shifts) > 0: 
            if shift_id not in shift_sets:
                shift_sets[shift_id] = [chain]
            else:
                shift_sets[shift_id].append(chain)
            out = open(file_basename+'_'+chain+'_'+str(i)+'.str','w')
            out.write('save_PROTEIN\n   _Entity.Polymer_seq_one_letter_code '+seq_single_letter+'\nsave_\n')
            out.write('loop_\n      _Entity_poly_seq.Hetero\n      _Entity_poly_seq.Mon_ID\n      _Entity_poly_seq.Num\n      _Entity_poly_seq.Comp_index_ID\n      _Entity_poly_seq.Entry_ID\n      _Entity_poly_seq.Entity_ID\n')
            for line in lines_seq:
                out.write(line+'\n')
            out.write('   stop_\n')
            out.write('loop_\n      _Atom_chem_shift.ID\n      _Atom_chem_shift.Seq_ID\n      _Atom_chem_shift.Comp_ID\n      _Atom_chem_shift.Atom_ID\n      _Atom_chem_shift.Val\n      _Atom_chem_shift.Val_err\n      _Atom_chem_shift.Assigned_chem_shift_list_ID\n')        
            for line in lines_shifts:
                out.write(line+'\n')
            out.write('   stop_\n')
            out.close()
    i+=1

c = 1
summary=''
for shift_id in shift_sets:
    summary += '\n -> set '+str(c)+', chain(s): '
    for chain in shift_sets[shift_id]:
        summary += chain+', '
    summary = summary[:-2]
    c+=1
    
print("found "+str(len(shift_sets))+' set(s) of shifts: '+summary)
