########################################### Dictionaries ###########################################

# Neighbouring residue chemical shift corrections(Tamiola 2010)
# A = i-2 
H_correction_A = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  
HA_correction_A = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  
C_correction_A = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  
CA_correction_A = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  
CB_correction_A = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  
N_correction_A = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  
# B = i+1 
H_correction_B = {'CYS': 0.293, 'ASP': 0.081, 'SER': 0.23, 'GLN': 0.225, 'LYS': 0.174, 'ILE': 0.25, 'PRO': 0.292, 'THR': 0.211, 'PHE': 0.079, 'ASN': 0.13, 'GLY': 0.0, 'HIS': 0.141, 'LEU': 0.105, 'ARG': 0.251, 'TRP': -0.555, 'ALA': 0.088, 'VAL': 0.259, 'GLU': 0.174, 'TYR': -0.052, 'MET': 0.205}
HA_correction_B = {'CYS': 0.025, 'ASP': -0.024, 'SER': 0.006, 'GLN': -0.009, 'LYS': -0.031, 'ILE': 0.01, 'PRO': -0.032, 'THR': -0.005, 'PHE': -0.04, 'ASN': -0.004, 'GLY': 0.0, 'HIS': -0.019, 'LEU': -0.042, 'ARG': -0.004, 'TRP': -0.168, 'ALA': -0.046, 'VAL': -0.007, 'GLU': -0.048, 'TYR': -0.035, 'MET': -0.018}
C_correction_B = {'CYS': -0.417, 'ASP': 0.059, 'SER': -0.128, 'GLN': -0.153, 'LYS': -0.117, 'ILE': -0.045, 'PRO': -0.04, 'THR': -0.13, 'PHE': -0.672, 'ASN': -0.052, 'GLY': 0.0, 'HIS': -0.243, 'LEU': -0.115, 'ARG': -0.246, 'TRP': -0.465, 'ALA': 0.068, 'VAL': -0.198, 'GLU': -0.032, 'TYR': -0.633, 'MET': -0.118}
CA_correction_B = {'CYS': -0.006, 'ASP': 0.29, 'SER': 0.155, 'GLN': 0.125, 'LYS': 0.195, 'ILE': -0.087, 'PRO': 0.021, 'THR': 0.143, 'PHE': -0.163, 'ASN': 0.293, 'GLY': 0.0, 'HIS': 0.069, 'LEU': 0.083, 'ARG': 0.029, 'TRP': 0.014, 'ALA': 0.161, 'VAL': 0.075, 'GLU': 0.171, 'TYR': -0.245, 'MET': 0.028}
CB_correction_B = {'CYS': 0.008, 'ASP': -0.182, 'SER': -0.153, 'GLN': -0.064, 'LYS': -0.102, 'ILE': -0.111, 'PRO': -0.151, 'THR': -0.072, 'PHE': 0.059, 'ASN': -0.205, 'GLY': 0.0, 'HIS': -0.093, 'LEU': -0.143, 'ARG': -0.038, 'TRP': -0.046, 'ALA': -0.135, 'VAL': -0.147, 'GLU': -0.119, 'TYR': 0.017, 'MET': -0.139}
N_correction_B = {'CYS': 3.1, 'ASP': 0.679, 'SER': 2.297, 'GLN': 1.7, 'LYS': 1.311, 'ILE': 4.205, 'PRO': 0.569, 'THR': 2.68, 'PHE': 2.238, 'ASN': 0.699, 'GLY': 0.0, 'HIS': 1.761, 'LEU': 0.835, 'ARG': 1.909, 'TRP': 0.795, 'ALA': -0.639, 'VAL': 4.507, 'GLU': 1.226, 'TYR': 2.729, 'MET': 1.407}
# C = i+1
H_correction_C = {'CYS': 0.375, 'ASP': -0.003, 'SER': -0.002, 'GLN': -0.016, 'LYS': -0.053, 'ILE': -0.071, 'PRO': -0.028, 'THR': 0.042, 'PHE': -0.133, 'ASN': -0.023, 'GLY': 0.0, 'HIS': -0.136, 'LEU': -0.086, 'ARG': -0.039, 'TRP': 0.032, 'ALA': -0.056, 'VAL': -0.052, 'GLU': -0.014, 'TYR': -0.093, 'MET': -0.066}
HA_correction_C = {'CYS': -0.008, 'ASP': -0.009, 'SER': 0.03, 'GLN': -0.039, 'LYS': -0.011, 'ILE': 0.0, 'PRO': 0.276, 'THR': 0.065, 'PHE': -0.052, 'ASN': -0.033, 'GLY': 0.0, 'HIS': -0.059, 'LEU': -0.021, 'ARG': -0.017, 'TRP': -0.082, 'ALA': -0.042, 'VAL': 0.021, 'GLU': -0.005, 'TYR': -0.044, 'MET': -0.018}
C_correction_C = {'CYS': -0.61, 'ASP': -0.853, 'SER': -0.483, 'GLN': -0.551, 'LYS': -0.52, 'ILE': -0.559, 'PRO': -2.874, 'THR': -0.327, 'PHE': -0.891, 'ASN': -0.849, 'GLY': 0.0, 'HIS': -0.677, 'LEU': -0.612, 'ARG': -0.625, 'TRP': -0.694, 'ALA': -0.816, 'VAL': -0.662, 'GLU': -0.503, 'TYR': -0.848, 'MET': -0.271}
CA_correction_C = {'CYS': 0.047, 'ASP': -0.095, 'SER': -0.151, 'GLN': -0.098, 'LYS': -0.142, 'ILE': -0.129, 'PRO': -2.423, 'THR': -0.122, 'PHE': -0.247, 'ASN': -0.067, 'GLY': 0.0, 'HIS': -0.248, 'LEU': -0.184, 'ARG': -0.258, 'TRP': -0.201, 'ALA': -0.187, 'VAL': -0.214, 'GLU': -0.071, 'TYR': -0.148, 'MET': -0.061}
CB_correction_C = {'CYS': -0.17, 'ASP': 0.079, 'SER': 0.087, 'GLN': 0.004, 'LYS': 0.029, 'ILE': 0.093, 'PRO': -0.426, 'THR': 0.086, 'PHE': 0.023, 'ASN': -0.013, 'GLY': 0.0, 'HIS': 0.041, 'LEU': -0.053, 'ARG': 0.046, 'TRP': 0.069, 'ALA': 0.018, 'VAL': 0.044, 'GLU': 0.065, 'TYR': 0.019, 'MET': -0.205}
N_correction_C = {'CYS': 0.785, 'ASP': -0.373, 'SER': -0.053, 'GLN': -0.087, 'LYS': 0.068, 'ILE': -0.167, 'PRO': 1.072, 'THR': 0.076, 'PHE': -0.484, 'ASN': -0.219, 'GLY': 0.0, 'HIS': -0.225, 'LEU': -0.397, 'ARG': -0.065, 'TRP': -0.273, 'ALA': -0.004, 'VAL': -0.047, 'GLU': -0.1, 'TYR': -0.387, 'MET': -0.145}
# D = i+2 
H_correction_D = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  
HA_correction_D = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  
C_correction_D = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  
CA_correction_D = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  
CB_correction_D = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  
N_correction_D = {"ALA":0.0,"CYS":0.0,"ASP":0.0,"GLU":0.0,"PHE":0.0,"GLY":0.0,"HIS":0.0,"ILE":0.0,"LYS":0.0,"LEU":0.0,"MET":0.0,"ASN":0.0,"PRO":0.0,"GLN":0.0,"ARG":0.0,"SER":0.0,"THR":0.0,"VAL":0.0,"TRP":0.0,"TYR":0.0}  

RC_correction = {'HA':[HA_correction_A, HA_correction_B, HA_correction_C, HA_correction_D],'H':[H_correction_A, H_correction_B, H_correction_C, H_correction_D],'CA':[CA_correction_A, CA_correction_B, CA_correction_C, CA_correction_D],'CB':[CB_correction_A, CB_correction_B, CB_correction_C, CB_correction_D],'C':[C_correction_A, C_correction_B, C_correction_C, C_correction_D],'N':[N_correction_A, N_correction_B, N_correction_C, N_correction_D]}

# random coil reference values (Tamiola 2010)
H_RC = {'CYS': 8.41, 'ASP': 8.217, 'SER': 8.215, 'GLN': 8.258, 'LYS': 8.221, 'ILE': 7.963, 'THR': 8.047, 'PHE': 8.107, 'ASN': 8.366, 'GLY': 8.307, 'HIS': 8.31, 'LEU': 8.088, 'ARG': 8.232, 'TRP': 7.725, 'ALA': 8.158, 'VAL': 8.037, 'GLU': 8.304, 'TYR': 8.026, 'MET': 8.209}
HA_RC = {'CYS': 4.447, 'ASP': 4.537, 'SER': 4.392, 'GLN': 4.254, 'LYS': 4.237, 'ILE': 4.076, 'PRO': 4.339, 'THR': 4.252, 'PHE': 4.573, 'ASN': 4.632, 'GLY': 3.98, 'HIS': 4.585, 'LEU': 4.26, 'ARG': 4.239, 'TRP': 4.567, 'ALA': 4.224, 'VAL': 4.009, 'GLU': 4.222, 'TYR': 4.504, 'MET': 4.425}
C_RC = {'CYS': 174.927, 'ASP': 176.987, 'SER': 175.236, 'GLN': 176.51, 'LYS': 177.224, 'ILE': 176.897, 'PRO': 177.542, 'THR': 175.122, 'PHE': 176.368, 'ASN': 175.825, 'GLY': 174.63, 'HIS': 175.349, 'LEU': 178.037, 'ARG': 176.821, 'TRP': 174.549, 'ALA': 178.418, 'VAL': 176.772, 'GLU': 177.125, 'TYR': 176.284, 'MET': 176.953}
CA_RC = {'CYS': 58.327, 'ASP': 54.331, 'SER': 58.352, 'GLN': 55.84, 'LYS': 56.412, 'ILE': 61.247, 'PRO': 63.18, 'THR': 61.926, 'PHE': 57.934, 'ASN': 53.231, 'GLY': 45.236, 'HIS': 55.964, 'LEU': 55.26, 'ARG': 56.088, 'TRP': 57.5, 'ALA': 52.599, 'VAL': 62.347, 'GLU': 56.65, 'TYR': 57.761, 'MET': 55.591}
CB_RC = {'CYS': 28.085, 'ASP': 41.089, 'SER': 63.766, 'GLN': 29.509, 'LYS': 32.921, 'ILE': 38.563, 'PRO': 32.072, 'THR': 69.794, 'PHE': 39.66, 'ASN': 38.79, 'MET': 32.69, 'HIS': 29.719, 'LEU': 42.212, 'ARG': 30.691, 'TRP': 29.38, 'ALA': 19.102, 'VAL': 32.674, 'GLU': 30.225, 'TYR': 38.75}
N_RC = {'CYS': 119.068, 'ASP': 120.207, 'SER': 115.935, 'GLN': 120.224, 'LYS': 121.353, 'ILE': 120.512, 'PRO': 136.612, 'THR': 114.024, 'PHE': 120.138, 'ASN': 118.668, 'GLY': 108.783, 'HIS': 118.93, 'LEU': 121.877, 'ARG': 121.288, 'TRP': 120.733, 'ALA': 123.906, 'VAL': 120.403, 'GLU': 120.769, 'TYR': 120.228, 'MET': 120.002}

RC_values = {'N':N_RC,'C':C_RC,'CA':CA_RC,'CB':CB_RC,'H':H_RC,'HA':HA_RC}

# RCI weighting coefficients (Whishart Server as of May 2018)
RCI_weights = {'CA':{'CA':0.6,'CB':0,'C':0,'N':0,'HA':0,'H':0}, 'C':{'CA':0,'CB':0,'C':0.3,'N':0,'HA':0,'H':0}, 'N':{'CA':0,'CB':0,'C':0,'N':0.6,'HA':0,'H':0}, 'CB':{'CA':0,'CB':0.6,'C':0,'N':0,'HA':0,'H':0}, 'H':{'CA':0,'CB':0,'C':0,'N':0,'HA':0,'H':0.6}, 'HA':{'CA':0,'CB':0,'C':0,'N':0,'HA':0.6,'H':0}, 'HN':{'CA':0,'CB':0,'C':0,'N':0.71,'HA':0,'H':0.51}, 'HAN':{'CA':0,'CB':0,'C':0,'N':0.27,'HA':0.77,'H':0}, 'CN':{'CA':0,'CB':0,'C':0.77,'N':0.27,'HA':0,'H':0}, 'HHA':{'CA':0,'CB':0,'C':0,'N':0,'HA':0.67,'H':0.03}, 'CHA':{'CA':0,'CB':0,'C':0.55,'N':0,'HA':0.69,'H':0}, 'CH':{'CA':0,'CB':0,'C':0.85,'N':0,'HA':0,'H':0.25}, 'CBN':{'CA':0,'CB':0.71,'C':0,'N':0.65,'HA':0,'H':0}, 'CBHA':{'CA':0,'CB':0.06,'C':0,'N':0,'HA':0.69,'H':0}, 'CBH':{'CA':0,'CB':0.8,'C':0,'N':0,'HA':0,'H':0.4}, 'CCB':{'CA':0,'CB':0.27,'C':0.77,'N':0,'HA':0,'H':0}, 'CAN':{'CA':0.77,'CB':0,'C':0,'N':0.27,'HA':0,'H':0}, 'CAHA':{'CA':0.6,'CB':0,'C':0,'N':0,'HA':0.69,'H':0}, 'CAH':{'CA':0.85,'CB':0,'C':0,'N':0,'HA':0,'H':0.25}, 'CCA':{'CA':0.71,'CB':0,'C':0.65,'N':0,'HA':0,'H':0}, 'CACB':{'CA':0.75,'CB':0.38,'C':0,'N':0,'HA':0,'H':0}, 'CHAN':{'CA':0,'CB':0,'C':0.72,'N':0.34,'HA':0.83,'H':0}, 'CHN':{'CA':0,'CB':0,'C':0.93,'N':0.27,'HA':0,'H':0.2}, 'HHAN':{'CA':0,'CB':0,'C':0,'N':0.27,'HA':0.77,'H':0}, 'CBHAN':{'CA':0,'CB':0.07,'C':0,'N':0.27,'HA':0.82,'H':0}, 'CBHN':{'CA':0,'CB':0.81,'C':0,'N':0.7,'HA':0,'H':0.26}, 'CHHA':{'CA':0,'CB':0,'C':0.76,'N':0,'HA':0.82,'H':0.18}, 'CCBN':{'CA':0,'CB':0.36,'C':0.82,'N':0.31,'HA':0,'H':0}, 'CBHHA':{'CA':0,'CB':0.07,'C':0,'N':0,'HA':0.76,'H':0.04}, 'CCBHA':{'CA':0,'CB':0.06,'C':0.61,'N':0,'HA':0.76,'H':0}, 'CCBH':{'CA':0,'CB':0.3,'C':0.87,'N':0,'HA':0,'H':0.23}, 'CAHAN':{'CA':0.75,'CB':0,'C':0,'N':0.37,'HA':0.82,'H':0}, 'CAHN':{'CA':0.83,'CB':0,'C':0,'N':0.2,'HA':0,'H':0.23}, 'CCAN':{'CA':0.79,'CB':0,'C':0.74,'N':0.32,'HA':0,'H':0}, 'CAHHA':{'CA':0.77,'CB':0,'C':0,'N':0,'HA':0.82,'H':0.2}, 'CCAHA':{'CA':0.64,'CB':0,'C':0.59,'N':0,'HA':0.78,'H':0}, 'CCAH':{'CA':0.76,'CB':0,'C':0.76,'N':0,'HA':0,'H':0.26}, 'CACBN':{'CA':0.8,'CB':0.43,'C':0,'N':0.32,'HA':0,'H':0}, 'CACBHA':{'CA':0.7,'CB':0.17,'C':0,'N':0,'HA':0.74,'H':0}, 'CACBH':{'CA':0.84,'CB':0.36,'C':0,'N':0,'HA':0,'H':0.2}, 'CCACB':{'CA':0.81,'CB':0.31,'C':0.73,'N':0,'HA':0,'H':0}, 'CCACBHA':{'CA':0.7,'CB':0.14,'C':0.63,'N':0,'HA':0.82,'H':0}, 'CCACBH':{'CA':0.83,'CB':0.34,'C':0.8,'N':0,'HA':0,'H':0.28}, 'CCACBN':{'CA':0.83,'CB':0.36,'C':0.75,'N':0.34,'HA':0,'H':0}, 'CACBHHA':{'CA':0.81,'CB':0.17,'C':0,'N':0,'HA':0.83,'H':0.18}, 'CACBHAN':{'CA':0.78,'CB':0.18,'C':0,'N':0.36,'HA':0.82,'H':0}, 'CACBHN':{'CA':0.92,'CB':0.4,'C':0,'N':0.22,'HA':0,'H':0.2}, 'CCAHHA':{'CA':0.73,'CB':0,'C':0.72,'N':0,'HA':0.85,'H':0.25}, 'CCAHAN':{'CA':0.68,'CB':0,'C':0.67,'N':0.42,'HA':0.85,'H':0}, 'CCAHN':{'CA':0.85,'CB':0,'C':0.81,'N':0.23,'HA':0,'H':0.25}, 'CAHHAN':{'CA':0.77,'CB':0,'C':0,'N':0.34,'HA':0.84,'H':0.1}, 'CCBHHA':{'CA':0,'CB':0.05,'C':0.78,'N':0,'HA':0.85,'H':0.18}, 'CCBHAN':{'CA':0,'CB':0.05,'C':0.74,'N':0.37,'HA':0.86,'H':0}, 'CCBHN':{'CA':0,'CB':0.29,'C':0.89,'N':0.26,'HA':0,'H':0.2}, 'CBHHAN':{'CA':0,'CB':0.07,'C':0,'N':0.27,'HA':0.82,'H':0}, 'CHHAN':{'CA':0,'CB':0,'C':0.77,'N':0.34,'HA':0.86,'H':0.07}, 'CCACBHHA':{'CA':0.79,'CB':0.13,'C':0.73,'N':0,'HA':0.86,'H':0.25}, 'CACBHHAN':{'CA':0.83,'CB':0.19,'C':0,'N':0.37,'HA':0.84,'H':0.1}, 'CCAHHAN':{'CA':0.76,'CB':0,'C':0.76,'N':0.41,'HA':0.89,'H':0.18}, 'CCBHHAN':{'CA':0,'CB':0.06,'C':0.8,'N':0.36,'HA':0.89,'H':0.08}, 'CCACBHAN':{'CA':0.74,'CB':0.15,'C':0.7,'N':0.47,'HA':0.87,'H':0}, 'CCACBHN':{'CA':0.87,'CB':0.39,'C':0.8,'N':0.24,'HA':0,'H':0.26}, 'CCACBHHAN':{'CA':0.81,'CB':0.15,'C':0.78,'N':0.43,'HA':0.9,'H':0.18}}

# dictionary to convert residue abbreviations
amino_acids = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}   
amino_acids_three_letter = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}   

