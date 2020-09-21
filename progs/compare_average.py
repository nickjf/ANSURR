import glob
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import percentileofscore

def plot(resi,RCI, FIRST,n):
    #helices_x = [FIRST['resi'][i[0]] for i in enumerate(FIRST['ss']) if i[1] == 'H'] # FIRST resi can be nan!
    #helices_y = [1.025]*len(helices_x)
    #strands_x = [FIRST['resi'][i[0]] for i in enumerate(FIRST['ss']) if i[1] == 'E']
    #strands_y = [1.025]*len(strands_x)
    plt.figure(n)
    plt.plot(resi,RCI)
    plt.plot(resi,FIRST)
    #plt.scatter(helices_x,helices_y,color='red',s=2)
    #plt.scatter(strands_x,strands_y,color='blue',s=2)
    plt.xlabel('residue number',size=12)
    plt.ylabel('flexibility',size=12)
    plt.ylim(0,1.05)
    plt.title('structure: ' + PDB_ID+chain_shiftID.split('_')[0]+'_average' + ' shifts: ' + SHIFT_ID+'_'+chain_shiftID.split('_')[1] + ' shift%: '+av_perc_shifts_out+'\n correlation score: ' + str(round(corr_score,1)) +' RMSD score: ' + str(round(RMSD_score,1)))
    plt.tight_layout()
    plt.savefig(PDB_ID+chain_shiftID.split('_')[0] +'_'+SHIFT_ID+chain_shiftID.split('_')[1]+'_average.png',dpi=300,bbox_inches='tight')
    
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
    chain_shiftID = os.path.basename(i).split(PDB_ID)[1][0].replace('_','')+'_'+os.path.basename(i).replace('.out','').split(SHIFT_ID)[1].replace('_','')
    if chain_shiftID not in av:
        av[chain_shiftID] = {}
    for line in open(i,'r'):
        line = line.split()
        resi = line[0]
        if resi not in av[chain_shiftID]:
            av[chain_shiftID][resi] = {'resn':line[1],'rci':[float(line[2])],'first':[float(line[3])],'shift_perc':[float(line[4])],'shifts':[line[5]]}
        else:
            av[chain_shiftID][resi]['rci'].append(float(line[2]))
            av[chain_shiftID][resi]['first'].append(float(line[3]))
            av[chain_shiftID][resi]['shift_perc'].append(float(line[4]))

n = 0
for chain_shiftID in av:
    for r in av[chain_shiftID]:
        if all_same([str(i) for i in av[chain_shiftID][r]['rci']]): #nanmean complains about lists with all nans
            av[chain_shiftID][r]['rci'] = av[chain_shiftID][r]['rci'][0]
        else:
            av[chain_shiftID][r]['rci'] = np.nanmean(av[chain_shiftID][r]['rci'])
        if all_same([str(i) for i in av[chain_shiftID][r]['first']]): #nanmean complains about lists with all nans
            av[chain_shiftID][r]['first'] = av[chain_shiftID][r]['first'][0]
        else:
            av[chain_shiftID][r]['first'] = np.nanmean(av[chain_shiftID][r]['first'])
        av[chain_shiftID][r]['shift_perc'] = np.nanmean(av[chain_shiftID][r]['shift_perc'])
    
    av_perc_shifts = int(round(100.0*np.nanmean([av[chain_shiftID][i]['shift_perc'] for i in av[chain_shiftID]]),0))
    
    resi = [int(i) for i in av[chain_shiftID]]
    RCI = [av[chain_shiftID][i]['rci'] for i in av[chain_shiftID]]
    FIRST = [av[chain_shiftID][i]['first'] for i in av[chain_shiftID]]
    shifts = [av[chain_shiftID][i]['shift_perc'] for i in av[chain_shiftID]]

    if av_perc_shifts < 75:
        av_perc_shifts_out = str(av_perc_shifts)+' (RCI may be unreliable!)'
        #print('*WARNING chemical shift completeness (' + str(av_perc_shifts) +'%)' +' is below recommended minimum (75%), RCI may be unreliable* DONE')
    else:
        av_perc_shifts_out = str(av_perc_shifts)
        #print('DONE')

    RCI_nan = [i for i,x in enumerate(RCI) if np.isnan(x)]
    FIRST_nan = [i for i,x in enumerate(FIRST) if np.isnan(x)]
    RCI_noRC = [x for i,x in enumerate(RCI) if i not in FIRST_nan and not np.isnan(x) and shifts[i] >=0.17 ]
    FIRST_noRC = [x for i,x in enumerate(FIRST) if i not in RCI_nan and not np.isnan(x) and shifts[i] >=0.17]

    RMSD_noRC =  calc_RMSD(RCI_noRC,FIRST_noRC)
    RMSD_score = 100.0 - percentileofscore(rmsd_benchmark,RMSD_noRC)

    if all_same(FIRST_noRC) or all_same(RCI_noRC):
        spearman_noRC = np.nan
        corr_score = np.nan
        #print('*WARNING Spearman correlation coefficient cannot be determined, setting correlation score to NaN * ',end='')
    else:
        spearman_noRC = spearmanr(RCI_noRC,FIRST_noRC)[0]
        corr_score = percentileofscore(corr_benchmark,spearman_noRC)

    # write output file
    RCI_FIRST_out = open(PDB_ID+chain_shiftID.split('_')[0] +'_'+SHIFT_ID+chain_shiftID.split('_')[1]+'_average.out','w')
    for i in enumerate(resi):
        RCI_FIRST_out.write('{:>5}'.format(i[1])+' '+'{:<23.20f}'.format(float(RCI[i[0]]))+'{:<23.20f}'.format(float(FIRST[i[0]]))+'\n')
    RCI_FIRST_out.close()

    # append to scores.out
    scores = open('scores_average.out','a+')
    scores.write('PDB: '+ '{:<11}'.format(PDB_ID+chain_shiftID.split('_')[0]) + ' SHIFTS: '+ '{:<11}'.format(SHIFT_ID+'_'+chain_shiftID.split('_')[1]) + ' SHIFT%: '+ '{:3}'.format(av_perc_shifts) + ' Spearman: '+ '{: 5.3f}'.format(spearman_noRC) + ' CorrelationScore: '+ '{:4.1f}'.format(round(corr_score,1)) + ' RMSD: '+ '{:4.3f}'.format(RMSD_noRC) + ' RMSDScore: '+ '{:4.1f}'.format(round(RMSD_score,1))+ '\n')
    scores.close()
    plot(resi,RCI,FIRST,n)
    n+=1                                                                  
