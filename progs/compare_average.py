import glob
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import percentileofscore

def plot(resi,RCI, FIRST):
    #helices_x = [FIRST['resi'][i[0]] for i in enumerate(FIRST['ss']) if i[1] == 'H'] # FIRST resi can be nan!
    #helices_y = [1.025]*len(helices_x)
    #strands_x = [FIRST['resi'][i[0]] for i in enumerate(FIRST['ss']) if i[1] == 'E']
    #strands_y = [1.025]*len(strands_x)

    plt.plot(resi,RCI)
    plt.plot(resi,FIRST)
    #plt.scatter(helices_x,helices_y,color='red',s=2)
    #plt.scatter(strands_x,strands_y,color='blue',s=2)
    plt.xlabel('residue number',size=12)
    plt.ylabel('flexibility',size=12)
    plt.ylim(0,1.05)
    plt.title('structure: ' + PDB_ID + ' shifts: ' + SHIFT_ID + ' shift%: '+av_perc_shifts_out+'\n correlation score: ' + str(round(corr_score,1)) +' RMSD score: ' + str(round(RMSD_score,1)))
    #plt.title('structure: '  + ' shifts: ' + ' shift%: '+av_perc_shifts_out+'\n correlation score: ' + str(round(corr_score,1)) +' RMSD score: ' + str(round(RMSD_score,1)))
    plt.savefig(PDB_ID+'_'+SHIFT_ID+'average.png',dpi=300,bbox_inches='tight')
    
def calc_RMSD(RCI_score_for_corr,FIRST_score_for_corr):
    diff = []
    for i in enumerate(RCI_score_for_corr):
        r = i[1]
        f = FIRST_score_for_corr[i[0]]
        diff.append(np.square(r-f))
    return np.sqrt(np.mean(diff))

def all_same(items):
    return all(x == items[0] for x in items)

PDB_ID = os.path.basename(sys.argv[1]).split('.')[0]
SHIFT_ID= os.path.basename(sys.argv[2]).split('.')[0]
ANSURR_PATH = sys.argv[4]
rmsd_benchmark_in = open(ANSURR_PATH+'/lib/benchmark_rmsd','r')
corr_benchmark_in = open(ANSURR_PATH+'/lib/benchmark_corr','r')

corr_benchmark = []
for line in corr_benchmark_in:
    corr_benchmark.append(float(line))
    
rmsd_benchmark = []
for line in rmsd_benchmark_in:
    rmsd_benchmark.append(float(line))

o = glob.glob('ANSURR_output/out/*.out')

av = {}
for i in o:
    for line in open(i,'r'):
        line = line.split()
        resi = line[0]
        if resi not in av:
            av[resi] = {'resn':line[1],'rci':[float(line[2])],'first':[float(line[3])],'shift_perc':[float(line[4])],'shifts':[line[5]]}
        else:
            av[resi]['rci'].append(float(line[2]))
            av[resi]['first'].append(float(line[3]))
            av[resi]['shift_perc'].append(float(line[4]))

for r in av:
    av[r]['rci'] = np.nanmean(av[r]['rci'])
    av[r]['first'] = np.nanmean(av[r]['first'])
    av[r]['shift_perc'] = np.nanmean(av[r]['shift_perc'])
    
av_perc_shifts = int(round(100.0*np.nanmean([av[i]['shift_perc'] for i in av]),0))
    
resi = [int(i) for i in av]
RCI = [av[i]['rci'] for i in av]
FIRST = [av[i]['first'] for i in av]
shifts = [av[i]['shifts'] for i in av]

if av_perc_shifts < 75:
    av_perc_shifts_out = str(av_perc_shifts)+' (RCI may be unreliable!)'
    print('*WARNING chemical shift completeness (' + str(av_perc_shifts) +'%)' +' is below recommended minimum (75%), RCI may be unreliable* DONE')
else:
    av_perc_shifts_out = str(av_perc_shifts)
    print('DONE')

RCI_nan = [i for i,x in enumerate(RCI) if np.isnan(x)]
FIRST_nan = [i for i,x in enumerate(FIRST) if np.isnan(x)]
RCI_noRC = [x for i,x in enumerate(RCI) if i not in FIRST_nan and not np.isnan(x) and shifts[i] != 'RC']
FIRST_noRC = [x for i,x in enumerate(FIRST) if i not in RCI_nan and not np.isnan(x) and shifts[i] != 'RC']

RMSD_noRC =  calc_RMSD(RCI_noRC,FIRST_noRC)
RMSD_score = 100.0 - percentileofscore(rmsd_benchmark,RMSD_noRC)

if all_same(FIRST_noRC) or all_same(RCI_noRC):
    spearman_noRC = np.nan
    corr_score = np.nan
    print('*WARNING Spearman correlation coefficient cannot be determined, setting correlation score to NaN * ',end='')
else:
    spearman_noRC = spearmanr(RCI_noRC,FIRST_noRC)[0]
    corr_score = percentileofscore(corr_benchmark,spearman_noRC)

plot(resi,RCI,FIRST)
