
from RCI_dicts import RC_values, RC_correction, RCI_weights, amino_acids
import sys
import os
import numpy as np

########################################### Objects ###########################################

class res(object):
    _registry = []
    def __init__(self,i,name):
        self._registry.append(self)
        self.i = i
        self.name = name
        self.shifts = [] 
        self.secondary_shifts = []
        self.secondary_shifts_smoothed = []
        self.neighbours = ['','','','']# positions refer to neighbour A,B,C,D
        self.assumed_RC = 0
        self.shift_types = ''
        self.RCI = 0
        self.RCI_smoothed = 0
        
class atom(object):
    _registry = []
    def __init__(self,atom_type,shift):
        self.atom_type = atom_type
        self.shift = shift
        
########################################### Functions ###########################################

def append_shift(res, atom_type, shift):
    try:
        shift = float(shift)
        res.shifts.append(atom(atom_type,shift))
    except:
        pass

def append_secshift(res, atom_type, shift):
    try:
        shift = float(shift)
        res.secondary_shifts.append(atom(atom_type,shift))
    except:
        pass
    
def check_version(file_in): # Check PANAV or NMRStarV3
    if '.panav' in file_in.name: # read in PANAV output
        extract_shifts_PANAV(file_in)
    else:
        #for line in file_in:
        #    if "_Entry.NMR_STAR_version" in line:
        #        extract_shifts_nmrstar3(file_in) # read in NMR-star3
        #        break
        extract_shifts_nmrstar3(file_in)

def extract_shifts_PANAV(file_in):
    for line in file_in:
        line = line.split('\t')
        resi = line[0].replace(' ','')[1:]
        try:
            resn = amino_acids[line[0].replace(' ','')[0]]
            resi = int(resi)
            r = res(resi,resn)
            if line[1] != '      ' or line[2] != '      ' or line[3] != '      ' or line[4] != '      ' or line[5] != '      ' or line[6] != '      ': 
                append_shift(r,'C',line[1])
                append_shift(r,'CA',line[2])
                append_shift(r,'CB',line[3])
                append_shift(r,'N',line[4])
                append_shift(r,'H',line[5])
                append_shift(r,'HA',line[6])
            else:       
                append_secshift(r,'C',0)
                append_secshift(r,'CA',0)
                append_secshift(r,'CB',0)
                append_secshift(r,'N',0)
                append_secshift(r,'H',0)
                append_secshift(r,'HA',0)
                r.assumed_RC = 1 
        except:
            pass

def extract_shifts_nmrstar3(file_in):
    seq = {}
    seq_labels = []
    while True:
        for line in file_in:
            if "_Entity.Polymer_seq_one_letter_code" in line:
                line = line.split()
                print(line)
                if len(line) = 2:
                    print(line[1])
    
    atoms = ['HA','HA2','HA3','H','C','CA','CB','N']
    resis = []
    for line in file_in:
        if 'stop_' in line:
            break
        elif line != '\n':
            line = line.split()
            atom_type = line[labels.index('Atom_ID')] 
            if atom_type in atoms:
                resi = int(line[labels.index('Seq_ID')])
                resn = line[labels.index('Comp_ID')]
                shift = float(line[labels.index('Val')])
                if resn in amino_acids.values():
                    for s in seq:
                        if s == resi:
                            break
                        if s not in [r.i for r in res._registry]:
                            r = res(s,seq[s])
                            append_secshift(r,'C',0)
                            append_secshift(r,'CA',0)
                            append_secshift(r,'CB',0)
                            append_secshift(r,'N',0)
                            append_secshift(r,'H',0)
                            append_secshift(r,'HA',0)
                            r.assumed_RC = 1
                    if resi not in resis:
                        r = res(resi,resn)
                    append_shift(r,atom_type,shift)
                    resis.append(resi)
    for r in res._registry:
        if r.name == 'GLY':
            average_gly(r)
                     
def average_gly(res):
    HA = []
    for s in res.shifts:
        if s.atom_type == 'HA2' or s.atom_type == 'HA3':
            HA.append(s.shift)
    if HA != []:
        append_shift(res,'HA',np.mean(HA))
        res.shifts = [s for s in res.shifts if s.atom_type not in ['HA2','HA3']]

def correct_shifts(res,neighbour_res,pos): # correct shifts based on position and neighbouring residue
    for s in res.shifts:
        #if s.atom_type != 'CB':
        s.shift -= RC_correction[s.atom_type][pos][neighbour_res] # I corrected exp shifts originally which is wrong, hence the minus sign here to fix
        
def calc_secondary_coil_shifts(res):
    for s in res.shifts: 
        res.secondary_shifts.append(atom(s.atom_type,np.absolute(s.shift - RC_values[s.atom_type][res.name])))
        
def floor_shifts(res):
    C_min = 0.04
    H_min = 0.01
    N_min = 0.1
    for s in res.secondary_shifts:
        if 'C' in s.atom_type:
            max(s.shift,C_min)
        elif 'H' in s.atom_type:
            max(s.shift,H_min)
        elif 'N' in s.atom_type:
            max(s.shift,N_min)
                    
def smooth_shifts(res):
    for s in res.secondary_shifts:
        shifts_temp = []
        shifts_temp.append(s.shift)
        for n in res.neighbours[1:3]:
            if n != '':
                for ns in n.secondary_shifts:
                    if s.atom_type == ns.atom_type:
                        shifts_temp.append(ns.shift)
        res.secondary_shifts_smoothed.append(atom(s.atom_type,np.mean(shifts_temp)))
        
def calc_RCI(res,weights):
    RCI = 0
    sum_weights = 0.0
    for s in res.secondary_shifts_smoothed:
        RCI += (10.0 * (weights[s.atom_type]) * s.shift)
        sum_weights += weights[s.atom_type]
    RCI = min(sum_weights / RCI,0.2)
        
    return RCI

def end_correction(RCI):
    end_length = 4
    RCI_corrected = []
    max_N = max(RCI[:end_length])
    max_N_pos = max([i for i, j in enumerate(RCI[:end_length]) if j == max_N])
    max_C = max(RCI[-end_length:])
    max_C_pos = max([i for i, j in enumerate(RCI[-end_length:]) if j == max_C]) + len(RCI) - end_length

    for r in enumerate(RCI): # Correct ends
        if r[0] in list(range(max_N_pos)):
            if (2*np.abs(max_N - r[1]) + r[1]) > 0.6:
                RCI_corrected.append(0.6)
            else:
                RCI_corrected.append((2*np.abs(max_N - r[1])) + r[1])
        elif r[0] in list(range(max_C_pos,len(RCI)+1)):
            if (2*np.abs(max_C - r[1]) + r[1]) > 0.6:
                RCI_corrected.append(0.6)
            else:
                RCI_corrected.append((2*np.abs(max_C - r[1])) + r[1])
        else:
            RCI_corrected.append(r[1])

    return RCI_corrected[0],RCI_corrected[1],RCI_corrected[2],RCI_corrected[3],RCI_corrected[-4],RCI_corrected[-3],RCI_corrected[-2],RCI_corrected[-1]

def smooth_RCI(res):
    RCI_smoothed = []
    RCI_smoothed.append(res.RCI)
    for n in res.neighbours[1:3]:
        if n != '':
            RCI_smoothed.append(n.RCI)
    res.RCI_smoothed = np.mean(RCI_smoothed)
    
########################################### Program starts here ###########################################

shifts_in = open(sys.argv[1],"r")                   # read in shift file
check_version(shifts_in)                            # check version NMR-Star or PANAV calibrated and extract shifts, make res object for each residue and append atom objects

for r in res._registry:
    for r2 in res._registry:                        # append neighbours
        if (r2.i == r.i-2):
            r.neighbours[0] = r2
        if (r2.i == r.i-1):
            r.neighbours[1] = r2
        if (r2.i == r.i+1):
            r.neighbours[2] = r2
        if (r2.i == r.i+2):
            r.neighbours[3] = r2
    for n in enumerate(r.neighbours):
        if n[1] != '':                              # skip if no neighbours i.e. at termini
            correct_shifts(r,n[1].name,n[0])        # apply RC correction according to neighbour residue types
    calc_secondary_coil_shifts(r)                   # calculate secondary shifts NB: RC correction already applied the step before
    floor_shifts(r)                                 # apply floor values to secondary shifts

for r in res._registry:
    if r.secondary_shifts != []:
        smooth_shifts(r)                            # 3-res smoothing of secondary shifts based on neighbours
        for s in r.secondary_shifts_smoothed:       # normalise secondary shifts
            if 'C' in s.atom_type:
                s.shift *= 2.5
            elif 'H' in s.atom_type:
                s.shift *= 10.0
            s.shift = max(s.shift,0.5)
        atom_list = []
        weight_hash = ''
        for s in r.secondary_shifts_smoothed:       # get atom_types for weight selection 
            atom_list.append(s.atom_type)
        atom_list = sorted(atom_list)
        for a in atom_list:                         # build weight hash
            weight_hash += a
        r.RCI = calc_RCI(r,RCI_weights[weight_hash])# calculate RCI
        r.num_shifts = len(atom_list)               # count number of shifts for output file
        r.shift_types = weight_hash                 # shift types for output file

                                                    # apply end_correction
res._registry[0].RCI,res._registry[1].RCI,res._registry[2].RCI,res._registry[3].RCI,res._registry[-4].RCI,res._registry[-3].RCI,res._registry[-2].RCI,res._registry[-1].RCI = end_correction([r.RCI for r in res._registry])

for r in res._registry:                             # smooth RCI values
    smooth_RCI(r)

out = open(os.path.basename(sys.argv[1]).split('.')[0]+'.rci','w')                          # write output
for r in res._registry:
    if r.assumed_RC == 1:
        out.write('{:>5}'.format(str(r.i)) + '{:>4}'.format(r.name) + ' ' + '{:<20.18f}'.format(r.RCI_smoothed) + '{:>2}'.format('0') + ' RC' + '\n')
    else:
        out.write('{:>5}'.format(str(r.i)) + '{:>4}'.format(r.name) + ' ' + '{:<20.18f}'.format(r.RCI_smoothed) + '{:>2}'.format(r.num_shifts) +' '+ r.shift_types + '\n')
out.close()
	
