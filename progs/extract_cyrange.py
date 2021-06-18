import sys, json

cyrange_output = sys.argv[1]

with open(cyrange_output, "r") as cyrange_output:
    cyrange_output = cyrange_output.readlines()

cyrange = {'default chain':[]}

for cyrange_out in cyrange_output:
    if "Optimal range" in cyrange_out:
        cyrange_out = cyrange_out[18:]
        gaps = int(cyrange_out.split(':')[1].split(',')[1].split()[0])
        if gaps == 0:
            if cyrange_out[0].isalpha():
                chain = cyrange_out[0]
                first_resi = int(cyrange_out.split(':')[0].split('..')[0][1:])
                last_resi = int(cyrange_out.split(':')[0].split('..')[1][1:])
                if chain not in cyrange:
                    cyrange[chain] = []
                cyrange[chain].extend(list(range(first_resi,last_resi+1)))
            else:
                first_resi = int(cyrange_out.split(':')[0].split('..')[0])
                last_resi = int(cyrange_out.split(':')[0].split('..')[1])
                cyrange['default chain'].extend(list(range(first_resi,last_resi+1)))
        else:
            for r in cyrange_out.split(':')[0].split(','):
                r = r.strip()
                if '..' in r:
                    if r[0].isalpha():
                        chain = r[0]
                        first_resi = int(r.split(':')[0].split('..')[0][1:])
                        last_resi = int(r.split(':')[0].split('..')[1][1:])
                        if chain not in cyrange:
                            cyrange[chain] = []
                        cyrange[chain].extend(list(range(first_resi,last_resi+1)))
                    else:
                        first_resi = int(r.split(':')[0].split('..')[0])
                        last_resi = int(r.split(':')[0].split('..')[1])
                        cyrange['default chain'].extend(list(range(first_resi,last_resi+1)))
                else:
                    if r[0].isalpha():
                        chain = r[0]
                        cyrange[chain].extend([int(r)])
                    else:
                        cyrange['default chain'].extend([int(r)])

chain_resi = {}
for chain in cyrange:
    resi = len(cyrange[chain])
    if resi > 0:
        chain_resi[chain] = resi
  
if len(chain_resi) > 0: 
    if len(chain_resi) == 1:
        print(' -> identified '+str(chain_resi['default chain'])+' well-defined residues')
    else:
        msg=' -> identified'
        for chain in chain_resi:
            msg+=' '+str(chain_resi[chain])+' ('+chain+'),'
        msg=msg[:-1]+" well-defined residues"
        print(msg)
    
    with open('cyrange.json', 'w') as outfile:
        json.dump(cyrange, outfile) 



        
