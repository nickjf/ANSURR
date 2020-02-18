import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import adjust_text as at
import sys

data = {}
for line in open('scores.out','r'):
    if line[:4] == 'PDB:':
        line = line.split()
        pdb = line[1][:[pos for pos, char in enumerate(line[1]) if char == '_'][-1]] # accounts for "_" appearing in structure name
        model = int(line[1].split('_')[-1])
        shifts = line[3]
        if pdb not in data:
            data[pdb] = {'model':[model],'rmsd':[float(line[9])],'corr':[float(line[7])],'shifts':shifts}
        else:
            data[pdb]['model'].append(model)
            data[pdb]['rmsd'].append(float(line[9]))
            data[pdb]['corr'].append(float(line[7]))
            
for d in enumerate(data):
    plt.figure(d[0])
    plt.scatter(data[d[1]]['rmsd'], data[d[1]]['corr'], c='black',s=8)
    labels = [plt.text(data[d[1]]['rmsd'][i], data[d[1]]['corr'][i], data[d[1]]['model'][i],
                       ha='center', va='center', color='blue',size=6) for i in range(len(data[d[1]]['rmsd']))]
    plt.ylabel('correlation score',size=12)
    plt.xlabel('RMSD score',size=12)
    plt.axis('scaled')  
    plt.xlim(-5,105)
    plt.ylim(-5,105)
    plt.title('Structure: '+d[1]+' Shifts: '+data[d[1]]['shifts'])
    at.adjust_text(labels,expand_text=(1.75, 1.75), expand_points=(1.75, 1.75),arrowprops=dict(arrowstyle='->', color='red', alpha=0.6, linewidth=0.5))
    plt.tight_layout()
    plt.savefig(d[1]+'_'+data[d[1]]['shifts']+'.png',dpi=300)
