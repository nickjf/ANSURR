
import json, sys

resi_ref = json.load(open("resi_ref.tmp"))
decomp = sys.argv[1]
decomp_num = sys.argv[1].split('/')[-1]
num = decomp.split('.decomp')[0].split('_')[-1]

chains = ''
for chain in resi_ref[num]:
    chains += chain
    

for chain in resi_ref[num]:
    out = open(sys.argv[2]+chain+'_'+num+'.decomp','w')
    orig_first = resi_ref[num][chain]['orig_first']
    orig_last = resi_ref[num][chain]['orig_last']
    new_first = resi_ref[num][chain]['new_first']
    new_last = resi_ref[num][chain]['new_last']
    for line in open(decomp,'r'):
        resi= int(line.split()[0]) 
        if resi >= new_first and resi <= new_last:
            out.write(str(resi - new_first + orig_first)+ ' '+line.split()[1]+' '+line.split()[2]+'\n')
        if resi == new_last:
            out.close()
            break
    
            