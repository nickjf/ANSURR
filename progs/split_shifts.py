import sys

shifts = sys.argv[1]
shiftID = shifts.split('.')[0]

details = 0
ents = {}
temp = []
for line in open(shifts,'r'):
    if "_Entity." in line and details == 0:
        details = 1
        temp.append(line)
    elif "save_" in line[:5] and details == 1:
        details = 0
        ents[entID] = [i for i in temp]
        temp = []
    elif details == 1:
        if "_Entity.ID" in line:
            entID = int(line.split()[1])
            ents[entID] = temp
        temp.append(line)            

temp = []
details = 0
ent_specific = {}
shifts_start = 0
count_to_entityID = -1 # variable to find where entityID appears in chem list (it can vary e.g for bmrb id 27632 compared to 25660)
entityID_found = 0
for line in open(shifts,'r'):
    if "_Assigned_chem_shift_list" in line and details == 0:
        details = 1
        temp.append(line)
    if "_Atom_chem_shift.ID" in line and shifts_start == 0 and entityID_found == 0:
        count_to_entityID += 1
    elif "_Atom_chem_shift.Entity_ID" in line and shifts_start == 0 and entityID_found == 0:
        entityID_found = 1   
    if count_to_entityID > -1 and entityID_found == 0:
        count_to_entityID += 1
    if "_Atom_chem_shift.Assigned_chem_shift_list_ID" in line and shifts_start == 0:
        shifts_start = 1
        temp.append(line)
    if "save_" in line[:5] and details == 1:
        for i in ent_specific:
            ents[i].extend(temp)
        for i in ent_specific:
            ents[i].extend(ent_specific[i])
        temp = []
        details = 0
        ent_specific = {}
        shifts_start = 0
        count_to_entityID = -1
        entityID_found = 0
    if details == 1 and shifts_start == 1:
        try:
            entID = int(line.split()[count_to_entityID])
            if entID not in ent_specific:
                ent_specific[entID] = [line]
            else:
                ent_specific[entID].append(line)
        except:
            pass
    elif details == 1:
        temp.append(line)

out_string = ' -> '+shiftID+'.str split into '
for e in ents:
    out = open(shiftID+'_'+str(e)+'.str','w')
    for l in ents[e]:
        out.write(l)
    out.close()
    out_string += shiftID+'_'+str(e)+'.str '
print(out_string)



        
        
        
