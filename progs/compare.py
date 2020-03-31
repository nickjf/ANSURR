from __future__ import print_function
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import percentileofscore

def compare_seq(seq1,seq2):
    score = 0
    count = 0
    for s1 in seq1:
        if s1 == seq2[count]:
            score += 1
        count += 1
    return score

def fix_missing_res(d,report_break):
    count = d['resi'][0]
    pos = 0
    for i in d['resi']:
        if i != count:
            if report_break == 1:
                print('*WARNING break in structure at residue '+str(count-1)+'* ', end='')
                report_break = 2
            for k in d:
                if k == 'resi':
                    d[k].insert(pos,count)
                elif k == 'resn':
                    d[k].insert(pos,'XXX')
                else:
                    d[k].insert(pos,np.nan)
        elif report_break == 2:
            report_break = 1
        count += 1
        pos += 1
    return d

def calc_score(seq1,seq2):
    score = 0
    count = 0
    for s1 in seq1:
        if s1 == seq2[count]:
            score += 1
        count += 1
    return score

def find_breaks(nums):
    ranges = sum((list(t) for t in zip(nums, nums[1:]) if t[0]+1 != t[1]), [])
    iranges = iter(nums[0:1] + ranges + nums[-1:])
    breaks = [[n,int(next(iranges))+1] for n in iranges]
    return breaks

def align(RCI,FIRST,cut_off=0.5):
    if len(RCI['resn']) <= len(FIRST['resn']): # "a" is always the shorter seq
        a = RCI['resn']
        b = FIRST['resn']
        a_score = RCI['score']
        b_score = FIRST['score'] 
    else:
        a = FIRST['resn']
        b = RCI['resn']
        a_score = FIRST['score']
        b_score = RCI['score']

    max_score_N = 0
    for i in range(len(b)): # start from N
        a_test = a[-i-1:len(b)+1]
        a_test.extend([np.nan]*(len(b)-i-1))
        a_test = a_test[::-1]
        a_test.extend([np.nan]*(len(b)-len(a_test)))
        a_test = a_test[::-1]  
        score = calc_score(a_test,b)
        if score > max_score_N:
            max_score_N = score
            a_score_test = a_score[-i-1:len(b_score)+1]
            a_score_test.extend([np.nan]*(len(b_score)-i-1))
            a_score_test = a_score_test[::-1]
            a_score_test.extend([np.nan]*(len(b_score)-len(a_score_test)))
            a_score_test = a_score_test[::-1]
            chosen_i = i

    max_score_C = max_score_N # then start from C and look for better fit
    for i in range(1,len(a)):
        a_test = a[:-i]
        a_test = a_test[::-1]
        a_test.extend([np.nan]*(len(b)-len(a_test)))
        a_test = a_test[::-1]   
        score = calc_score(a_test,b)
        if score > max_score_C:
            max_score_C = score
            a_score_test = a_score[:-i]
            a_score_test = a_score_test[::-1]
            a_score_test.extend([np.nan]*(len(b_score)-len(a_score_test)))
            a_score_test = a_score_test[::-1]
            chosen_i = i

    if max_score_N >= max_score_C:
        if len(RCI['resi']) <= len(FIRST['resi']):
            for r in RCI:  
                RCI[r] = RCI[r][-chosen_i-1:len(b_score)+1] 
                RCI[r].extend([np.nan]*(len(b_score)-chosen_i-1))
                RCI[r] = RCI[r][::-1]   
                RCI[r].extend([np.nan]*(len(b_score)-len(RCI[r])))
                RCI[r] = RCI[r][::-1]
        else:
            for f in FIRST:
                FIRST[f] = FIRST[f][-chosen_i-1:len(b_score)+1] 
                FIRST[f].extend([np.nan]*(len(b_score)-chosen_i-1)) 
                FIRST[f] = FIRST[f][::-1]   
                FIRST[f].extend([np.nan]*(len(b_score)-len(FIRST[f])))
                FIRST[f] = FIRST[f][::-1]
        num_resn = len([i for i in FIRST['resn'] if i != 'XXX'])
        if max_score_N < cut_off * num_resn:
            print('sequence identity ('+str(round(100*max_score_N/num_resn,1))+'%) is below cut-off ('+str(100*cut_off)+'%), skipping')
            quit()
    else:
        if len(RCI['resi']) <= len(FIRST['resi']):
            for r in RCI:
                RCI[r] = RCI[r][:-chosen_i]
                RCI[r] = RCI[r][::-1]
                RCI[r].extend([np.nan]*(len(b_score)-len(RCI[r])))
                RCI[r] = RCI[r][::-1]
        else:
            for f in FIRST:
                FIRST[f] = FIRST[f][:-chosen_i]
                FIRST[f] = FIRST[f][::-1]
                FIRST[f].extend([np.nan]*(len(b_score)-len(FIRST[f])))
                FIRST[f] = FIRST[f][::-1]
        num_resn = len([i for i in FIRST['resn'] if i != 'XXX'])
        if max_score_C < cut_off * num_resn:
            print('sequence identity ('+str(round(100*max_score_C/num_resn,1))+'%) is below cut-off ('+str(100*cut_off)+'%), skipping')
            quit()
    RCI['resi'] = FIRST['resi']
    return RCI, FIRST

def trim(res1,res2): # lazt variable names fix this!
    s = 0
    e = 0
    sc = 0
    ec = 0
    for i in enumerate(res1):
        if str(i[1]) != 'nan': 
            if sc == 3:
                if i[1] != res2[i[0]]:
                    ec += 1
                else:
                    ec = 0
                if ec == 3:
                    e = i[0]
            else:
                if i[1] == res2[i[0]]:
                    sc += 1
                else:
                    sc = 0
                if sc == 3:
                    s = i[0]
    s = s-2
    if ec == 3:
        e = e-2
    else:
        e = len(res1) - ec
    return s,e
   
def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def get_segments(resi):
    prev_res = min(resi) if resi else None
    segments = list()
    for number in enumerate(resi):
        if number[1] != prev_res+1:
            segments.append([number[0]])
        elif len(segments[-1]) > 1:
            segments[-1][-1] = number[0]
        else:
            segments[-1].append(number[0])
        prev_res = number[1]
    return segments

def smooth(resi,data):
    data_smoothed = []
    segs = get_segments(resi)
    for s in segs:
        if len(s) > 2:
            data_smoothed_temp = movingaverage(data[s[0]:s[1]+1],3)
            N = (data[s[0]] + data[s[0]+1]) / 2.0
            data_smoothed_temp[0] = N
            C = (data[s[1]] + data[s[1]-1]) / 2.0
            data_smoothed_temp[-1] = C
            data_smoothed.extend(data_smoothed_temp)
        elif len(s) == 2:
            data_smoothed_temp = movingaverage(data[s[0]:s[1]+1],2)
            N = (data[s[0]] + data[s[0]+1]) / 2.0
            data_smoothed_temp[0] = N
            C = (data[s[1]] + data[s[1]-1]) / 2.0
            data_smoothed_temp[-1] = C
            data_smoothed.extend(data_smoothed_temp)
        else:
            data_smoothed.append(data[s[0]])
    return data_smoothed

def calc_RMSD(RCI_score_for_corr,FIRST_score_for_corr):
    diff = []
    for i in enumerate(RCI_score_for_corr):
        r = i[1]
        f = FIRST_score_for_corr[i[0]]
        diff.append(np.square(r-f))
    return np.sqrt(np.mean(diff))

def rescale_FIRST(FIRST):
    K = 315777.09
    T = 298.15
    Hartree = 0.00038
    FIRST['score'] = [np.exp((4.2* x * Hartree)*K/T) if not np.isnan(x) else np.nan for x in FIRST['score']]
    FIRST_nan = [i for i,x in enumerate(FIRST['score']) if np.isnan(x)]
    FIRST['score'] = smooth([FIRST['resi'][i] for i,x in enumerate(FIRST['score']) if not np.isnan(x)],[i for i in FIRST['score'] if not np.isnan(i)])
    for n in FIRST_nan:
        FIRST['score'].insert(n,np.nan)
    return FIRST['score']

def rescale_RCI(RCI_score): 
    offset = 0.024
    scale = 0.2-0.024
    RCI_score = [min(max((i-offset)/scale,0.0),1.0) for i in RCI_score]
    return RCI_score

def plot(RCI, FIRST):
    helices_x = [FIRST['resi'][i[0]] for i in enumerate(FIRST['ss']) if i[1] == 'H'] # FIRST resi can be nan!
    helices_y = [1.025]*len(helices_x)
    strands_x = [FIRST['resi'][i[0]] for i in enumerate(FIRST['ss']) if i[1] == 'E']
    strands_y = [1.025]*len(strands_x)

    plt.plot(RCI['resi'],RCI['score'])
    plt.plot(FIRST['resi'],FIRST['score'])
    plt.scatter(helices_x,helices_y,color='red',s=2)
    plt.scatter(strands_x,strands_y,color='blue',s=2)
    plt.xlabel('residue number',size=12)
    plt.ylabel('flexibility',size=12)
    plt.ylim(0,1.05)
    plt.title('structure: ' + PDB_ID + ' shifts: ' + SHIFT_ID + '\n correlation score: ' + str(round(spearman_noRC,1)) +' RMSD score: ' + str(round(RMSD_noRC,1)))
    plt.savefig(PDB_ID+'_'+SHIFT_ID+'.png',dpi=300,bbox_inches='tight')

def all_same(items):
    return all(x == items[0] for x in items)
    
# import data
FIRST_in = open(sys.argv[1],'r')
RCI_in = open(sys.argv[2],'r')
PDB_ID = os.path.basename(sys.argv[1]).split('.')[0]
SHIFT_ID= os.path.basename(sys.argv[2]).split('.')[0]
ANSURR_PATH = sys.argv[4]
rmsd_benchmark_in = open(ANSURR_PATH+'/lib/benchmark_rmsd','r')
corr_benchmark_in = open(ANSURR_PATH+'/lib/benchmark_corr','r')
print(" -> "+PDB_ID + '|'+ SHIFT_ID+' ',end='')

# secondary structure
ss_dict = {}
try:
    secondary_structure_in = open(sys.argv[3],'r')

    for line in secondary_structure_in:
        line = line.split()
        try:
            ss_dict[int(line[0])] = line[2]
        except:
            ss_dict[int(line[0])] = np.nan
except:
    pass

# read in data
RCI = {'resi':[],'resn':[],'score':[],'shifts':[],'shift_types':[]}
FIRST = {'resi':[],'resn':[],'score':[],'ss':[]}

for line in RCI_in:
    line = line.split()
    RCI['resi'].append(int(line[0]))
    RCI['resn'].append(line[1])
    RCI['score'].append(float(line[2]))
    if line[1] == 'GLY' or line[1] == 'PRO':
        RCI['shifts'].append(float(line[3])/5.0)
        RCI['shift_types'].append(line[4])
    else:
        RCI['shifts'].append(float(line[3])/6.0)
        RCI['shift_types'].append(line[4])
        
for line in FIRST_in:
    line = line.split()
    FIRST['resi'].append(int(line[0]))
    FIRST['resn'].append(line[1])
    FIRST['score'].append(float(line[2]))
    try:
        FIRST['ss'].append(ss_dict[int(line[0])])
    except:
        FIRST['ss'].append(np.nan)

corr_benchmark = []
for line in corr_benchmark_in:
    corr_benchmark.append(float(line))
    
rmsd_benchmark = []
for line in rmsd_benchmark_in:
    rmsd_benchmark.append(float(line))
    
# fix missing residues
FIRST = fix_missing_res(FIRST,1)
RCI = fix_missing_res(RCI,0)

# align
RCI,FIRST = align(RCI,FIRST)

# trim - look for where first 3 residues and last 3 residues in aligned RCI/FIRST output are and make these the beginning/end
start,end = trim(RCI['resn'],FIRST['resn'])
for r in RCI:
    RCI[r] = RCI[r][start:end] 
for f in FIRST:
    FIRST[f] = FIRST[f][start:end] 

# rescale
FIRST['score'] = rescale_FIRST(FIRST)
RCI['score'] = rescale_RCI(RCI['score'])

# various stats
RCI_nan = [i for i,x in enumerate(RCI['score']) if np.isnan(x)]
FIRST_nan = [i for i,x in enumerate(FIRST['score']) if np.isnan(x)]
RCI_noRC = [x for i,x in enumerate(RCI['score']) if i not in FIRST_nan and not np.isnan(x) and RCI['shift_types'][i] != 'RC']
FIRST_noRC = [x for i,x in enumerate(FIRST['score']) if i not in RCI_nan and not np.isnan(x) and RCI['shift_types'][i] != 'RC']

RMSD_noRC =  100.0 - percentileofscore(rmsd_benchmark,calc_RMSD(RCI_noRC,FIRST_noRC))

if all_same(FIRST_noRC) or all_same(RCI_noRC):
    spearman_noRC = np.nan
    print('*WARNING Spearman correlation coefficient cannot be determined, setting correlation score to NaN * ',end='')
else:
    spearman_noRC = percentileofscore(corr_benchmark,spearmanr(RCI_noRC,FIRST_noRC)[0])

av_perc_shifts = int(round(np.nanmean([RCI['shifts'][i[0]] for i in enumerate(RCI['resi']) if not np.isnan(i[1])])*100.0))

if av_perc_shifts < 75:
    print('*WARNING chemical shift completeness (' + str(av_perc_shifts) +'%)' +' is below recommended minimum (75%), RCI values may be unreliable* DONE')
else:
    print('DONE')

# write output file
RCI_FIRST_out = open(PDB_ID+'_'+SHIFT_ID+'.out','w')
for i in enumerate(FIRST['resi']):
    RCI_FIRST_out.write('{:>5}'.format(FIRST['resi'][i[0]])+' '+'{:<4}'.format(FIRST['resn'][i[0]])+'{:<23.20f}'.format(float(RCI['score'][i[0]]))+'{:<23.20f}'.format(float(FIRST['score'][i[0]]))+'{:<5.2f}'.format(RCI['shifts'][i[0]])+'{:<10}'.format(RCI['shift_types'][i[0]])+'\n')
RCI_FIRST_out.close()

# append to scores.out
scores = open('scores.out','a+')
scores.write('PDB: '+ '{:<11}'.format(PDB_ID) + ' SHIFTS: '+ '{:<12}'.format(SHIFT_ID) + ' SHIFT%: '+ '{:<4}'.format(av_perc_shifts) + ' CorrelationScore: '+ '{:<5.1f}'.format(spearman_noRC) + ' RMSDScore: '+ '{:<5.1f}'.format(RMSD_noRC)+ '\n')
scores.close()

# plot figs
plot(RCI,FIRST)
