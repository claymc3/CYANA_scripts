## 03/12/2020 Mary Clay PhD
## Script for preparing cyana input files 
## name.seq
## name_dihed.aco
## hbond.upl/hbond.lol
## name.upl

import pandas as pd 
import numpy as np
import os
import sys
import glob
import re


pdb_columns = ['name', 'resn', 'resid', 'X', 'Y', 'Z','nuc']
# Read in PDB file one line at a time, if the first four letter ar ATOM or HETA then it will parse the data into the 
# data frame, using the atome index int he PDB as the row index in the data frame. 
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "PTR":'Y', "TPO":"T", "SEP":'S' }
A_dict = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR' }
Atoms_dict = {'I':['CD1'], 'L':['CD1','CD2'], 'V':['CG1','CG2'], 'M':['CE'], 'A':['CB'], 'T':['CG2'], 'W':['NE1'], 'F':['CE1','CE2'], 'Y':['CE1','CE2']}

cwd = os.getcwd()
#####
if len(sys.argv)==1:
	print('''

Usage: 
	PrepCyan [Sequence] [index(s)] [pdb] [upl_extras] [residues] [TALOS] 

Required Input:
	Sequence 		Fasta style or three letter code formats accepted no need to
					remove header. Fusions indicated by + before start of fused 
					sequence. If this is not located in current directory
					provide the path. 

	Index(s)		Starting index for provided sequence. For fusions provide 
					starting indexes separated by commas (eg 1,230). The + will 
					trigger the use of the correct starting index. 

	PDB 			Reference PDB used to generate initial upl file. If this is 
					not located in current directory provide path 

	UPL Extras 		Can specify use of side chain heavy atoms of I, L, V, M, A,
					T, Y, F, or W to use in upl generation. 
					I = I-CD1, L = L-CD1, L-CD2, V = V-CG1, V-CG2, M = M-CE, 
					A = A-CB, T = T-CG2 F = F-CE1, F-CE2, Y = Y-CE1, Y-CE2, 
					W = W-NE1

	residues 		Residues to use in upl generation and rmsd calumniation.
					20-200,300-490

	TALOS 			Path to results of TALOSN analysis. The pred.tab will be 
					used to generate the dihid.aco input file. The predss.tab
					will be used to filter out unstructured regions from the 
					upl, hbond.lol, and hbond.upl files 


Assuming you have saved your .prot and .peak files are in current location

name is taken from the prot file in this directory and any peak files in 
this directory will be listing int he CYANA.calc file. 

OutPut:
	CYAN.calc
	init.cya
	name.seq
	name.upl 
	name_hbond.upl
	name_hbond.lol
	name_dihed.aco
''')
	exit()



szFileName = sys.argv[1]
Indexes = [int(x) for x in sys.argv[2].split(',')]
in_pdb  = sys.argv[3]
atoms = sys.argv[4]
residues = sys.argv[5]
TALOSdir = sys.argv[6]

prot = glob.glob('*.prot')[0]
name = prot.split('.')[0]
talosSS = os.path.join(TALOSdir +'/predSS.tab')

indihed = os.path.join(TALOSdir +'/pred.tab')
print("Using " + talosSS)
print("Using " + indihed)
# Indexes = [int(x) for x in index.split(',')]

#------------------------------------------------------------------------------
# Read in sequence and generate name.seq file 
#
sflines = [line.rstrip() for line in open(szFileName).readlines() if line.rstrip() and ">" not in line and "#" not in line]
Sequence = ''
### For 3 letter code sequences numbered or not
if re.search('[A-Z][A-Z][A-Z] * ([0-9]*)', sflines[0]) or len(sflines[0]) == 3:
	for line in sflines:
		if re.search('[A-Z][A-Z][A-Z] * ([0-9]*)', line):
			if re.search('[A-Z][A-Z][A-Z] * ([0-9]*)', line).group().split()[0].upper() in AAA_dict.keys():
				Sequence = Sequence + AAA_dict[re.search('[A-Z][A-Z][A-Z] * ([0-9]*)', line).group().split()[0].upper()]
		if re.search('[A-Z][A-Z][A-Z]', line) and not re.search('[A-Z][A-Z][A-Z] * ([0-9]*)', line):
			if re.search('[A-Z][A-Z][A-Z]', line).group().upper() in AAA_dict.keys():
				Sequence = Sequence + AAA_dict[re.search('[A-Z][A-Z][A-Z]', line).group().upper()]
		if '+' in line:
			Sequence = Sequence + '+'
### For single letter code fasta style sequences 
if len(sflines[0]) != 3 and not re.search('[A-Z][A-Z][A-Z] * ([0-9]*)', sflines[0]):
	for line in sflines:
		Sequence = Sequence + line
Sequence = Sequence.replace('++','+')

seqfile = open(name + '.seq','w')
for x in range(len(Sequence.split('+'))):
	if x != 0:
		seqfile.write("## + \n")
	for i in range(len(Sequence.split('+')[x])):
		index = i + int(Indexes[x])
		seqfile.write('%s %5d\n'% (A_dict[Sequence.split('+')[x][i]], index))
seqfile.close()
#
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Read in predSS.tab from TALOS run and great dictionary indicating secondary
# structure of residue to use for preparation of hbond and manual upl/lol 
#
talos_lines = [line.strip() for line in open(talosSS).readlines() if line.strip() and re.search(' *([0-9]*) [A-Z]', line)]
SecStrDict = {}
for line in talos_lines:
	res = line.split()[1] + line.split()[0]
	SecStrDict[res] = line.split()[-1].upper().replace('C','L')

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# Read in TALOS pred.tab and create dihed.aco file 

# dihed = glob.glob("*.aco")[0]
dihed_lines = [line.strip() for line in open(indihed).readlines() if line.strip() and re.search(' *([0-9]*) [A-Z]', line)]
# (resid, slc, phi, psi, dphi, dpsi, s2, count, cs_cound, Class)
aco = open('dihed.aco','w')
aco.write("\n")
scale = {"Strong":2.0, "Generous":3.0}
for line in dihed_lines:
	dline = line.split()
	if dline[-1] in scale.keys():
		dphi = float(dline[4])
		if dphi<10: dphi = 15.0
		if dphi>35: dphi = 35
		dpsi = float(dline[5])
		if dpsi<10: dpsi = 15.0
		if dpsi>35: dpsi = 35.0
		aco.write("#  " + line + "\n")
		if dline[1] != 'P':
			#print "%4d  %4s  PHI  %8.1f%8.1f" % (int(dline[0]), A_dict[dline[1]], float(dline[2])-scale[dline[-1]]*dphi, float(dline[2])+ scale[dline[-1]]*dphi)
			aco.write("%4d  %4s  PHI  %8.1f%8.1f\n" % (int(dline[0]), A_dict[dline[1]], float(dline[2])-scale[dline[-1]]*dphi/2, float(dline[2])+scale[dline[-1]]*dphi/2))
		#print "%4d  %4s  PSI  %8.1f%8.1f\n" % (int(dline[0]), A_dict[dline[1]], float(dline[3])-scale[dline[-1]]*dpsi, float(dline[2])+scale[dline[-1]]*dpsi)
		aco.write("%4d  %4s  PSI  %8.1f%8.1f\n\n" % (int(dline[0]), A_dict[dline[1]], float(dline[3])-scale[dline[-1]]*dpsi/2, float(dline[3])+scale[dline[-1]]*dpsi/2))
aco.close()

#------------------------------------------------------------------------------
# Fetch the .peak files and create a list 
#
peaks = glob.glob(os.path.join(cwd + '/*.peaks'))
peaks_out = ''
for peak in peaks:
	peaks_out = peaks_out + peak.split('/')[-1] + ','
peaks_out=peaks_out[:-1]

#------------------------------------------------------------------------------
# Generate list of allowed residue numbers
#
resn = residues.split(',')
allowed_resi = []
for i in range(len(resn)):
	allowed_resi.extend(range(int(resn[i].split('-')[1])+1)[int(resn[i].split('-')[0]):])

#------------------------------------------------------------------------------
# Generate list of allowed atoms in addition to N and O used for hydrogen bond limit generation 
#
allowed_atoms =["ALA-O", "ALA-N", "ARG-O", "ARG-N", "ASN-O", "ASN-N", "ASP-O", "ASP-N", "CYS-O", "CYS-N", "GLU-O", "GLU-N", "GLN-O", "GLN-N", "GLY-O", "GLY-N", "HIS-O", "HIS-N", "ILE-O", "ILE-N", "LEU-O", "LEU-N", "LYS-O", "LYS-N", "MET-O", "MET-N", "PHE-O", "PHE-N", "PRO-O", "PRO-N", "SER-O", "SER-N", "THR-O", "THR-N", "TRP-O", "TRP-N", "TYR-O", "TYR-N", "VAL-O", "VAL-N"]

for i in range(len(atoms)):
	for x in range(len(Atoms_dict[atoms[i]])):
		allowed_atoms.append(A_dict[atoms[i]] + '-' + Atoms_dict[atoms[i]][x])

#------------------------------------------------------------------------------
# Generating pandas data fram from PDB file using only specified residues and atoms
#
PDB_df = pd.DataFrame(columns=pdb_columns)
with open(in_pdb) as In_pdb:
	for line in In_pdb:
		if line[0:4] == "ATOM" or line[0:4] == 'HETA':
			if int(line[22:26].strip()) in allowed_resi:
				group = line[17:20].strip() + '-' + line[12:16].strip()
				if group in allowed_atoms:
					index =  AAA_dict[line[17:20].strip()] + line[22:26].strip() + "-" + line[12:16].strip()
					PDB_df.loc[index, 'name'] = line[12:16].strip()
					PDB_df.loc[index, 'resn'] = line[16:20].strip()
					PDB_df.loc[index, 'resid'] = int(line[22:26].strip())
					PDB_df.loc[index, 'X'] = float(line[30:38])
					PDB_df.loc[index, 'Y'] = float(line[38:46])
					PDB_df.loc[index, 'Z'] = float(line[46:54])
					PDB_df.loc[index, 'nuc'] = line[12:16].strip()[0]
					PDB_df.loc[index, 'chain'] = line[21]
					if AAA_dict[line[17:20].strip()] + line[22:26].strip() in SecStrDict.keys():
						PDB_df.loc[index, 'SecStr'] = SecStrDict[AAA_dict[line[17:20].strip()] + line[22:26].strip()]

#------------------------------------------------------------------------------
# Find Hydrogen bonding connections 
#
constrained = []
hblol = open('hbond.lol','w')
hbupl = open('hbond.upl', 'w')
### For Helical residues from TALOS
hbond_N = PDB_df[(PDB_df['nuc'] == 'N') & (PDB_df['SecStr'] == 'H')].index.tolist()
hbond_O = PDB_df[(PDB_df['nuc'] == 'O') & (PDB_df['SecStr'] == 'H')].index.tolist()
hblol.write("## Helical residues from TALOS \n")
hbupl.write("## Helical residues from TALOS \n")
for O in hbond_O:
	for N in hbond_N:
		if abs(PDB_df.loc[O,'resid'] - PDB_df.loc[N,'resid']) == 4: 
			dist = round(np.sqrt(((PDB_df.loc[O,'X'] - PDB_df.loc[N,'X'])**2) + ((PDB_df.loc[O,'Y'] - PDB_df.loc[N,'Y'])**2) + ((PDB_df.loc[O,'Z'] - PDB_df.loc[N,'Z'])**2)),1)
			if dist >= 2.7 and dist <= 3.3: 
				constrained.append(str(PDB_df.loc[O,'resid']) + '-' + str(PDB_df.loc[N,'resid']))
				hbupl.write("%4s %s  %-3s   %4s %s  %-3s     3.10\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
				hbupl.write("%4s %s  %-3s   %4s %s  %-3s     2.10\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))
				hblol.write("%4s %s  %-3s   %4s %s  %-3s     2.70\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
				hblol.write("%4s %s  %-3s   %4s %s  %-3s     1.80\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))
### For Betta Sheet residues from TALOS
hblol.write("## Betta Sheet residues from TALOS \n")
hbupl.write("## Betta Sheet residues from TALOS \n")				
hbond_N = PDB_df[(PDB_df['nuc'] == 'N') & (PDB_df['SecStr'] == 'E')].index.tolist()
hbond_O = PDB_df[(PDB_df['nuc'] == 'O') & (PDB_df['SecStr'] == 'E')].index.tolist()
for O in hbond_O:
	for N in hbond_N:
		if abs(PDB_df.loc[O,'resid'] - PDB_df.loc[N,'resid']) != 0: 
			dist = round(np.sqrt(((PDB_df.loc[O,'X'] - PDB_df.loc[N,'X'])**2) + ((PDB_df.loc[O,'Y'] - PDB_df.loc[N,'Y'])**2) + ((PDB_df.loc[O,'Z'] - PDB_df.loc[N,'Z'])**2)),1)
			if dist >= 2.7 and dist <= 3.3: 
				constrained.append(str(PDB_df.loc[O,'resid']) + '-' + str(PDB_df.loc[N,'resid']))
				hbupl.write("%4s %s  %-3s   %4s %s  %-3s     3.10\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
				hbupl.write("%4s %s  %-3s   %4s %s  %-3s     2.10\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))
				hblol.write("%4s %s  %-3s   %4s %s  %-3s     2.70\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
				hblol.write("%4s %s  %-3s   %4s %s  %-3s     1.80\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))

# ### For Helical residues NOT in TALOS
# hblol.write("## Helical residues NOT in TALOS\n")
# hbupl.write("## Helical residues NOT in TALOS \n")
# hbond_N = PDB_df[(PDB_df['nuc'] == 'N')].index.tolist()
# hbond_O = PDB_df[(PDB_df['nuc'] == 'O')].index.tolist()
# for O in hbond_O:
# 	for N in hbond_N:
# 		if str(PDB_df.loc[O,'resid']) + '-' + str(PDB_df.loc[N,'resid']) not in constrained:
# 			if abs(PDB_df.loc[O,'resid'] - PDB_df.loc[N,'resid']) == 4: 
# 				dist = round(np.sqrt(((PDB_df.loc[O,'X'] - PDB_df.loc[N,'X'])**2) + ((PDB_df.loc[O,'Y'] - PDB_df.loc[N,'Y'])**2) + ((PDB_df.loc[O,'Z'] - PDB_df.loc[N,'Z'])**2)),1)
# 				if dist >= 2.7 and dist <= 3.3: 
# 					constrained.append(str(PDB_df.loc[O,'resid']) + '-' + str(PDB_df.loc[N,'resid']))
# 					hbupl.write("%4s %s  %-3s   %4s %s  %-3s     3.40\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
# 					hbupl.write("%4s %s  %-3s   %4s %s  %-3s     2.40\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))
# 					hblol.write("%4s %s  %-3s   %4s %s  %-3s     2.70\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
# 					hblol.write("%4s %s  %-3s   %4s %s  %-3s     1.80\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))
# ### For Betta Sheet residues 
# hblol.write("## Betta Sheet residues not in TALOS \n")
# hbupl.write("## Betta Sheet residues not in TALOS \n")				
# for O in hbond_O:
# 	for N in hbond_N:
# 		if str(PDB_df.loc[O,'resid']) + '-' + str(PDB_df.loc[N,'resid']) not in constrained:
# 			if abs(PDB_df.loc[O,'resid'] - PDB_df.loc[N,'resid']) != 0: 
# 				dist = round(np.sqrt(((PDB_df.loc[O,'X'] - PDB_df.loc[N,'X'])**2) + ((PDB_df.loc[O,'Y'] - PDB_df.loc[N,'Y'])**2) + ((PDB_df.loc[O,'Z'] - PDB_df.loc[N,'Z'])**2)),1)
# 				if dist >= 2.7 and dist <= 3.3: 
# 					constrained.append(str(PDB_df.loc[O,'resid']) + '-' + str(PDB_df.loc[N,'resid']))
# 					hbupl.write("%4s %s  %-3s   %4s %s  %-3s     3.40\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
# 					hbupl.write("%4s %s  %-3s   %4s %s  %-3s     2.40\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))
# 					hblol.write("%4s %s  %-3s   %4s %s  %-3s     2.70\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
# 					hblol.write("%4s %s  %-3s   %4s %s  %-3s     1.80\n" % (PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))
hblol.close()
hbupl.close()
hblol.close()
hbupl.close()
print("Generated hbond.upl and hbond.lol")
#------------------------------------------------------------------------------
# Prepare upl constraints for initial calculation only considering structured elements
# based on TALOS predSS.tab
#
NN_used = []
upl = open(name + '.upl','w')
upl.write('## N-N distances from structured regions\n')
N_list = PDB_df[(PDB_df['nuc'] == 'N') & (PDB_df['SecStr'] != 'L')].index.tolist()
NA_list = PDB_df[(PDB_df['nuc'] == 'N')].index.tolist()
C_list = PDB_df[(PDB_df['nuc'] == 'C') & (PDB_df['SecStr'] != 'L')].index.tolist()
CA_list = PDB_df[(PDB_df['nuc'] == 'C')].index.tolist()
for i in range(len(N_list)):
	for x in range(len(N_list)):
		diff = abs(int(PDB_df.loc[N_list[i],'resid']) - int(PDB_df.loc[N_list[x],'resid']))
		if diff >= 3:
			dist = np.sqrt(((PDB_df.loc[N_list[i],'X'] - PDB_df.loc[N_list[x],'X'])**2) + ((PDB_df.loc[N_list[i],'Y'] - PDB_df.loc[N_list[x],'Y'])**2) + ((PDB_df.loc[N_list[i],'Z'] - PDB_df.loc[N_list[x],'Z'])**2))
			if dist < 6.0:
				constraint = N_list[x] + '-' + N_list[i]
				if constraint not in NN_used:
					NN_used.append(N_list[i] + '-' + N_list[x])
					NN_out = "%4s %s  %-3s   %4s %s  %-3s   %6.2f\n" % (PDB_df.loc[N_list[i],'resid'],PDB_df.loc[N_list[i],'resn'],PDB_df.loc[N_list[i],'name'], PDB_df.loc[N_list[x],'resid'],PDB_df.loc[N_list[x],'resn'],PDB_df.loc[N_list[x],'name'],dist)
					upl.write(NN_out)
upl.write('## N-N distances from Unstructured regions\n')
for i in range(len(NA_list)):
	for x in range(len(NA_list)):
		if NA_list[i] + '-' + NA_list[x] not in NN_used:
			diff = abs(int(PDB_df.loc[NA_list[i],'resid']) - int(PDB_df.loc[NA_list[x],'resid']))
			if diff >= 3:
				dist = np.sqrt(((PDB_df.loc[NA_list[i],'X'] - PDB_df.loc[NA_list[x],'X'])**2) + ((PDB_df.loc[NA_list[i],'Y'] - PDB_df.loc[NA_list[x],'Y'])**2) + ((PDB_df.loc[NA_list[i],'Z'] - PDB_df.loc[NA_list[x],'Z'])**2))
				if dist < 6.0:
					constraint = NA_list[x] + '-' + NA_list[i]
					if constraint not in NN_used:
						constraint2 = NA_list[i] + '-' + NA_list[x]
						NN_used.append(constraint2)
						NN_out = "%4s %s  %-3s   %4s %s  %-3s   %6.2f\n" % (PDB_df.loc[NA_list[i],'resid'],PDB_df.loc[NA_list[i],'resn'],PDB_df.loc[NA_list[i],'name'], PDB_df.loc[NA_list[x],'resid'],PDB_df.loc[NA_list[x],'resn'],PDB_df.loc[NA_list[x],'name'],dist)
						upl.write(NN_out)


print("Made %3.0f NN upl constraints " % (len(NN_used)))
NC_used = []
upl.write('## N-C distances from structured regions\n')
for i in range(len(N_list)):
	for x in range(len(C_list)):
		diff = abs(int(PDB_df.loc[N_list[i],'resid']) - int(PDB_df.loc[C_list[x],'resid']))
		if diff >= 3:
			dist = np.sqrt(((PDB_df.loc[N_list[i],'X'] - PDB_df.loc[C_list[x],'X'])**2) + ((PDB_df.loc[N_list[i],'Y'] - PDB_df.loc[C_list[x],'Y'])**2) + ((PDB_df.loc[N_list[i],'Z'] - PDB_df.loc[C_list[x],'Z'])**2))
			if dist < 6.0:
				NC_used.append(N_list[i] + '-' + C_list[x])
				NC_out = "%4s %s  %-3s   %4s %s  %-3s   %6.2f\n" % (PDB_df.loc[N_list[i],'resid'],PDB_df.loc[N_list[i],'resn'],PDB_df.loc[N_list[i],'name'], PDB_df.loc[C_list[x],'resid'],PDB_df.loc[C_list[x],'resn'],PDB_df.loc[C_list[x],'name'],dist)
				upl.write(NC_out)

upl.write('## N-C distances from Unstructured regions\n')
for i in range(len(NA_list)):
	for x in range(len(CA_list)):
		if NA_list[i] + '-' + CA_list[x] not in NC_used:
			diff = abs(int(PDB_df.loc[NA_list[i],'resid']) - int(PDB_df.loc[CA_list[x],'resid']))
			if diff >= 3:
				dist = np.sqrt(((PDB_df.loc[NA_list[i],'X'] - PDB_df.loc[CA_list[x],'X'])**2) + ((PDB_df.loc[NA_list[i],'Y'] - PDB_df.loc[CA_list[x],'Y'])**2) + ((PDB_df.loc[NA_list[i],'Z'] - PDB_df.loc[CA_list[x],'Z'])**2))
				if dist < 6.0:
					NC_used.append(NA_list[i] + '-' + CA_list[x])
					NC_out = "%4s %s  %-3s   %4s %s  %-3s   %6.2f\n" % (PDB_df.loc[NA_list[i],'resid'],PDB_df.loc[NA_list[i],'resn'],PDB_df.loc[NA_list[i],'name'], PDB_df.loc[CA_list[x],'resid'],PDB_df.loc[CA_list[x],'resn'],PDB_df.loc[CA_list[x],'name'],dist)
					upl.write(NC_out)

print("Made %3.0f NC upl constraints " % (len(NC_used)))
CC_used = []
upl.write('## C-C distances from structured regions\n')
for i in range(len(C_list)):
	for x in range(len(C_list)):
		diff = abs(int(PDB_df.loc[C_list[i],'resid']) - int(PDB_df.loc[C_list[x],'resid']))
		if diff >= 3:
			dist = np.sqrt(((PDB_df.loc[C_list[i],'X'] - PDB_df.loc[C_list[x],'X'])**2) + ((PDB_df.loc[C_list[i],'Y'] - PDB_df.loc[C_list[x],'Y'])**2) + ((PDB_df.loc[C_list[i],'Z'] - PDB_df.loc[C_list[x],'Z'])**2))
			if dist < 6.0:
				constraint = C_list[x] + '-' + C_list[i]
				if constraint not in CC_used:
					CC_used.append(C_list[i] + '-' + C_list[x])
					CC_out = "%4s %s  %-3s   %4s %s  %-3s   %6.2f\n" % (PDB_df.loc[C_list[i],'resid'],PDB_df.loc[C_list[i],'resn'],PDB_df.loc[C_list[i],'name'], PDB_df.loc[C_list[x],'resid'],PDB_df.loc[C_list[x],'resn'],PDB_df.loc[C_list[x],'name'],dist)
					upl.write(CC_out)
upl.write('## C-C distances from UNstructured regions\n')
for i in range(len(CA_list)):
	for x in range(len(CA_list)):
		if CA_list[i] + '-' + CA_list[x] not in CC_used:
			diff = abs(int(PDB_df.loc[CA_list[i],'resid']) - int(PDB_df.loc[CA_list[x],'resid']))
			if diff >= 3:
				dist = np.sqrt(((PDB_df.loc[CA_list[i],'X'] - PDB_df.loc[CA_list[x],'X'])**2) + ((PDB_df.loc[CA_list[i],'Y'] - PDB_df.loc[CA_list[x],'Y'])**2) + ((PDB_df.loc[CA_list[i],'Z'] - PDB_df.loc[CA_list[x],'Z'])**2))
				if dist < 6.0:
					constraint = CA_list[x] + '-' + CA_list[i]
					if constraint not in CC_used:
						CC_used.append(CA_list[i] + '-' + CA_list[x])
						CC_out = "%4s %s  %-3s   %4s %s  %-3s   %6.2f\n" % (PDB_df.loc[CA_list[i],'resid'],PDB_df.loc[CA_list[i],'resn'],PDB_df.loc[CA_list[i],'name'], PDB_df.loc[CA_list[x],'resid'],PDB_df.loc[CA_list[x],'resn'],PDB_df.loc[CA_list[x],'name'],dist)
						upl.write(CC_out)

print("Made %3.0f CC upl constraints " % (len(CC_used)))
upl.close()

ssa = open('ssa.cya'.'w')
ssa.close()

CYANA = open('CALC.cya','w')
clac_text = '''
peaks       := %s      # names of NOESY peak lists
prot        := %s.prot                   # names of chemical shift lists
constraints := %s.upl,%s,%s,%s              # additional (non-NOE) constraints
tolerance   := 0.01,0.01,0.01,0.01         # chemical shift tolerances
calibration :=                          # NOE calibration parameters
structures  := 100,20                   # number of initial, final structures
steps       := 20000                    # number of torsion angle dynamics steps
rmsdrange   := %s                       # residue range for RMSD calculation
randomseed  := 52541                   # random number generator seed

upl_values  := 2.4,7.0
#weight_rdc   = 0.02               # weight for RDC restraints
#cut_rdc      = 0.2                # cutoff for RDC violation output

#subroutine KEEP
#   peaks select "*, * number=10000..40000"
#end

ssa
noeassign peaks=$peaks prot=$prot autoaco # keep=KEEP 
''' % (peaks_out, name, name, 'hbond.lol', 'hbond.upl',name + '_dihed.aco',residues.replace("-",".."))
CYANA.write(clac_text)
CYANA.close()

Init_Cyana = open('init.cya','w')
init_text = '''name:=%s
rmsdrange:=%s
cyanalib
read lib special.lib append
#read lib cyana_Zn2.lib append
nproc:=20
read seq $name.seq
#molecules define 1..146 201..346
#molecule identity
#weight_ide=0.09
#molecule symdist "CA 1..146" "CA 201..346"
#weight_sym=0.025
#rdcdistances        # default bond lengths for dipole types''' % (name,residues.replace('-','...'))
Init_Cyana.write(init_text)
Init_Cyana.close()

special = open('special.lib','w')
special_text = '''ATOMTYPES    20      0.10
	 1 PSEUD    -10.00    0    0
	 2 H_ALI      1.00    0    1
	 3 H_AMI      0.95    1    1
	 4 H_ARO      1.00    0    1
	 5 H_SUL      1.00    1    1
	 6 H_OXY      1.00    1    1
	 7 C_ALI      1.60    0    6
	 8 C_BYL      1.50    0    6
	 9 C_ARO      1.60    0    6
	10 C_VIN      1.60    0    6
	11 N_AMI      1.45    0    7
	12 N_AMO      1.50   -1    7
	13 O_BYL      1.30   -1    8
	14 O_HYD      1.30   -1    8
	15 O_EST      1.30   -1    8
	16 S_OXY      1.80    0   16
	17 S_RED      1.80    0   16
	18 P_ALI      1.80    0   15
	19 METAL      1.80    0   50
	20 DUMMY    -10.00    0  999
 

RESIDUE   HIST     5   21    3   20
	 1 OMEGA    0    0    0.0000    2    1    3    4    0
	 2 PHI      0    0    0.0000    1    3    5   19    0
	 3 CHI1     0    0    0.0000    3    5    7   11   18
	 4 CHI2     0    0    0.0000    5    7   11   12   18
	 5 PSI      0    0    0.0000    3    5   19   21    0
	 1 C    C_BYL    0    0.0000    0.0000    0.0000    0.0000    2    3    0    0    0
	 2 O    O_BYL    0    0.0000   -0.6701    0.0000   -1.0332    1    0    0    0    0
	 3 N    N_AMI    0    0.0000    1.3293   -0.0000    0.0000    1    4    5    0    0
	 4 H    H_AMI    0    0.0000    1.8070    0.0001    0.8554    3    0    0    0    0
	 5 CA   C_ALI    0    0.0000    2.0938   -0.0009   -1.2422    3    6    7   19    0
	 6 HA   H_ALI    0    0.0000    2.6919    0.8975   -1.2641    5    0    0    0    0
	 7 CB   C_ALI    0    0.0000    3.0198   -1.2164   -1.2950    5    8    9   11    0
	 8 HB2  H_ALI    0    0.0000    2.7704   -1.8133   -2.1598    7    0    0    0   10
	 9 HB3  H_ALI    0    0.0000    4.0426   -0.8784   -1.3795    7    0    0    0   10
	10 QB   PSEUD    0    0.0000    3.4065   -1.3458   -1.7697    0    0    0    0    0
	11 CG   C_VIN    0    0.0000    2.9232   -2.0959   -0.0866    7   12   13    0    0
	12 ND1  N_AMO    0    0.0000    2.0787   -1.8317    0.9710   11   14    0    0    0
	13 CD2  C_ARO    0    0.0000    3.5702   -3.2425    0.2279   11   15   16    0    0
	14 CE1  C_ARO    0    0.0000    2.2117   -2.7768    1.8844   12   15   17    0    0
	15 NE2  N_AMO    0    0.0000    3.1105   -3.6454    1.4575   13   14   18    0    0
	16 HD2  H_ARO    0    0.0000    4.3113   -3.7475   -0.3758   13    0    0    0    0
	17 HE1  H_ARO    0    0.0000    1.6766   -2.8310    2.8206   14    0    0    0    0
	18 HE2  H_AMI    0    0.0000    3.3996   -4.4461    1.9422   15    0    0    0    0
	19 C    C_BYL    0    0.0000    1.1644   -0.0016   -2.4521    5   20   21    0    0
	20 O    O_BYL    0    0.0000    1.6186   -0.0031   -3.5967   19    0    0    0    0
	21 N    N_AMI    0    0.0000   -0.1390   -0.0014   -2.1911   19    0    0    0    0

RESIDUE   HIS+     5   22    3   21
	 1 OMEGA    0    0    0.0000    2    1    3    4    0
	 2 PHI      0    0    0.0000    1    3    5   20    0
	 3 CHI1     0    0    0.0000    3    5    7   11   19
	 4 CHI2     0    0    0.0000    5    7   11   12   19
	 5 PSI      0    0    0.0000    3    5   20   22    0
	 1 C    C_BYL    0    0.0000    0.0000    0.0000    0.0000    2    3    0    0    0
	 2 O    O_BYL    0    0.0000   -0.6700   -0.0000   -1.0333    1    0    0    0    0
	 3 N    N_AMI    0    0.0000    1.3294    0.0000    0.0000    1    4    5    0    0
	 4 H    H_AMI    0    0.0000    1.8069   -0.0000    0.8555    3    0    0    0    0
	 5 CA   C_ALI    0    0.0000    2.0939   -0.0005   -1.2421    3    6    7   20    0
	 6 HA   H_ALI    0    0.0000    2.6919    0.8979   -1.2638    5    0    0    0    0
	 7 CB   C_ALI    0    0.0000    3.0200   -1.2160   -1.2951    5    8    9   11    0
	 8 HB2  H_ALI    0    0.0000    2.7707   -1.8127   -2.1601    7    0    0    0   10
	 9 HB3  H_ALI    0    0.0000    4.0428   -0.8779   -1.3794    7    0    0    0   10
	10 QB   PSEUD    0    0.0000    3.4068   -1.3453   -1.7698    0    0    0    0    0
	11 CG   C_VIN    0    0.0000    2.9234   -2.0957   -0.0868    7   12   13    0    0
	12 ND1  N_AMO    0    0.0000    2.0788   -1.8318    0.9708   11   14   15    0    0
	13 CD2  C_ARO    0    0.0000    3.5705   -3.2423    0.2275   11   16   17    0    0
	14 HD1  H_AMI    0    0.0000    1.4735   -1.0649    1.0419   12    0    0    0    0
	15 CE1  C_ARO    0    0.0000    2.2118   -2.7769    1.8841   12   16   18    0    0
	16 NE2  N_AMO    0    0.0000    3.1107   -3.6453    1.4571   13   15   19    0    0
	17 HD2  H_ARO    0    0.0000    4.3118   -3.7471   -0.3762   13    0    0    0    0
	18 HE1  H_ARO    0    0.0000    1.6767   -2.8314    2.8203   15    0    0    0    0
	19 HE2  H_AMI    0    0.0000    3.3998   -4.4460    1.9417   16    0    0    0    0
	20 C    C_BYL    0    0.0000    1.1646   -0.0011   -2.4522    5   21   22    0    0
	21 O    O_BYL    0    0.0000    1.6189   -0.0024   -3.5968   20    0    0    0    0
	22 N    N_AMI    0    0.0000   -0.1389   -0.0011   -2.1913   20    0    0    0    0

RESIDUE   PROU     2   20    3   19
	 1 OMEGA    0    0    0.0000    2    1    3    4    0
	 2 PSI      0    0    0.0000    3    5   18   20    0
	 1 C    C_BYL    0    0.0000    0.0000    0.0000    0.0000    2    3    0    0    0
	 2 O    O_BYL    0    0.0000   -0.6561    0.0000   -1.0418    1    0    0    0    0
	 3 N    N_AMI    0    0.0000    1.3413    0.0000    0.0000    1    4    5    0    0
	 4 CD   C_ALI    0    0.0000    2.1884    0.0000    1.2045    3   11   15   16    0
	 5 CA   C_ALI    0    0.0000    2.1329   -0.0006   -1.2338    3    6    7   18    0
	 6 HA   H_ALI    0    0.0000    2.1323    0.9710   -1.7060    5    0    0    0    0
	 7 CB   C_ALI    0    0.0000    3.5459   -0.3336   -0.7485    5    8    9   11    0
	 8 HB2  H_ALI    0    0.0000    3.7115   -1.3995   -0.8189    7    0    0    0   10
	 9 HB3  H_ALI    0    0.0000    4.2709    0.1897   -1.3535    7    0    0    0   10
	10 QB   PSEUD    0    0.0000    3.9914   -0.6041   -1.0860    0    0    0    0    0
	11 CG   C_ALI    0    0.0000    3.5800    0.1317    0.6663    4    7   12   13    0
	12 HG2  H_ALI    0    0.0000    4.2544   -0.4858    1.2398   11    0    0    0   14
	13 HG3  H_ALI    0    0.0000    3.8907    1.1653    0.7050   11    0    0    0   14
	14 QG   PSEUD    0    0.0000    4.0726    0.3402    0.9724    0    0    0    0    0
	15 HD2  H_ALI    0    0.0000    2.0852   -0.9272    1.7486    4    0    0    0   17
	16 HD3  H_ALI    0    0.0000    1.9534    0.8409    1.8403    4    0    0    0   17
	17 QD   PSEUD    0    0.0000    2.0193   -0.0430    1.7951    0    0    0    0    0
	18 C    C_BYL    0    0.0000    1.6513   -1.0488   -2.2311    5   19   20    0    0
	19 O    O_BYL    0    0.0000    1.7796   -0.8720   -3.4421   18    0    0    0    0
	20 N    N_AMI    0    0.0000    1.0964   -2.1400   -1.7139   18    0    0    0    0

RESIDUE   PROO     6   20    3   19
	 1 OMEGA    0    0    0.0000    2    1    3    4    0
	 2 CHI3     0    0    0.0000    1    3    4    5    7
	 3 PHI      0    0    0.0000    1    3    8   18    0
	 4 CHI1     0    0    0.0000    3    8   10   14   17
	 5 CHI2     0    0    0.0000    8   10   14   15   17
	 6 PSI      0    0    0.0000    3    8   18   20    0
	 1 C    C_BYL    0    0.0000    0.0000    0.0000    0.0000    2    3    0    0    0
	 2 O    O_BYL    0    0.0000   -0.6561    0.0000   -1.0418    1    0    0    0    0
	 3 N    N_AMI    0    0.0000    1.3413    0.0000    0.0000    1    4    8    0    0
	 4 CD   C_ALI    0    0.0000    2.1884    0.0004    1.2045    3    5    6   14    0
	 5 HD2  H_ALI    0    0.0000    1.8873   -0.7796    1.8880    4    0    0    0    7
	 6 HD3  H_ALI    0    0.0000    2.1624    0.9634    1.6929    4    0    0    0    7
	 7 QD   PSEUD    0    0.0000    2.0249    0.0921    1.7911    0    0    0    0    0
	 8 CA   C_ALI    0    0.0000    2.1329   -0.0010   -1.2338    3    9   10   18    0
	 9 HA   H_ALI    0    0.0000    1.8441    0.8061   -1.8911    8    0    0    0    0
	10 CB   C_ALI    0    0.0000    3.5628    0.2285   -0.7380    8   11   12   14    0
	11 HB2  H_ALI    0    0.0000    4.2546   -0.3236   -1.3586   10    0    0    0   13
	12 HB3  H_ALI    0    0.0000    3.7967    1.2817   -0.7775   10    0    0    0   13
	13 QB   PSEUD    0    0.0000    4.0254    0.4796   -1.0674    0    0    0    0    0
	14 CG   C_ALI    0    0.0000    3.5636   -0.2772    0.6633    4   10   15   16    0
	15 HG2  H_ALI    0    0.0000    3.7625   -1.3380    0.6708   14    0    0    0   17
	16 HG3  H_ALI    0    0.0000    4.3073    0.2496    1.2428   14    0    0    0   17
	17 QG   PSEUD    0    0.0000    4.0352   -0.5439    0.9570    0    0    0    0    0
	18 C    C_BYL    0    0.0000    2.0353   -1.3249   -1.9841    8   19   20    0    0
	19 O    O_BYL    0    0.0000    1.0984   -1.5471   -2.7510   18    0    0    0    0
	20 N    N_AMI    0    0.0000    3.0080   -2.2009   -1.7577   18    0    0    0    0

RESIDUE   PCA      1   17    1   16
	 1 PSI      0    0    0.0000    1    3   15   17    0
	 1 N    N_AMI    0    0.0000    0.0000    0.0000    0.0000    2    3   13    0    0
	 2 H    H_AMI    0    0.0000   -0.4918    0.0000   -0.8376    1    0    0    0    0
	 3 CA   C_ALI    0    0.0000    1.4581    0.0000    0.0000    1    4    5   15    0
	 4 HA   H_ALI    0    0.0000    1.8644    0.9834   -0.0185    3    0    0    0    0
	 5 CB   C_ALI    0    0.0000    1.7448   -0.6589    1.3591    3    6    7    9    0
	 6 HB2  H_ALI    0    0.0000    2.7220   -0.3882    1.7345    5    0    0    0    8
	 7 HB3  H_ALI    0    0.0000    1.6387   -1.7332    1.3338    5    0    0    0    8
	 8 QB   PSEUD    0    0.0000    2.1804   -1.0607    1.5342    0    0    0    0    0
	 9 CG   C_ALI    0    0.0000    0.6171   -0.0187    2.2053    5   10   11   13    0
	10 HG2  H_ALI    0    0.0000    0.8536    0.9863    2.5047    9    0    0    0   12
	11 HG3  H_ALI    0    0.0000    0.3482   -0.6224    3.0589    9    0    0    0   12
	12 QG   PSEUD    0    0.0000    0.6009    0.1820    2.7818    0    0    0    0    0
	13 CD   C_BYL    0    0.0000   -0.5160   -0.0066    1.1838    1    9   14    0    0
	14 OE   O_BYL    0    0.0000   -1.6966   -0.0140    1.4547   13    0    0    0    0
	15 C    C_BYL    0    0.0000    1.8951   -0.7085   -1.2687    3   16   17    0    0
	16 O    O_BYL    0    0.0000    1.6547   -0.2537   -2.3656   15    0    0    0    0
	17 N    N_AMI    0    0.0000    2.5253   -1.8067   -1.0445   15    0    0    0    0

RESIDUE   PTR      7   31    3   30
	 1 OMEGA    0    0    0.0000    2    1    3    4    0
	 2 PHI      0    0    0.0000    1    3    5   29    0
	 3 CHI1     0    0    0.0000    3    5    7   14   28
	 4 CHI2     0    0    0.0000    5    7   14   15   28
	 5 CHI6     0    0    0.0000   17   19   24   25   28
	 6 CHI7     0    0    0.0000   19   24   25   28   28
	 7 PSI      0    0    0.0000    3    5   29   31    0
	 1 C    C_BYL    0    0.0000    0.0000    0.0000    0.0000    2    3    0    0    0
	 2 O    O_BYL    0    0.0000   -0.6703    0.0000   -1.0328    1    0    0    0    0
	 3 N    N_AMI    0    0.0000    1.3283    0.0000    0.0000    1    4    5    0    0
	 4 H    H_AMI    0    0.0000    1.8067    0.0012    0.8552    3    0    0    0    0
	 5 CA   C_ALI    0    0.0000    2.0926    0.0007   -1.2417    3    6    7   29    0
	 6 HA   H_ALI    0    0.0000    2.6339    0.9338   -1.2975    5    0    0    0    0
	 7 CB   C_ALI    0    0.0000    3.0946   -1.1550   -1.2490    5    8    9   14    0
	 8 HB2  H_ALI    0    0.0000    3.7086   -1.0846   -2.1335    7    0    0    0   10
	 9 HB3  H_ALI    0    0.0000    3.7232   -1.0829   -0.3735    7    0    0    0   10
	10 QB   PSEUD    0    0.0000    3.7159   -1.0838   -1.2535    0    0    0    0    0
	11 QD   PSEUD    0    0.0000    2.3703   -2.6769   -1.2423    0    0    0    0   13
	12 QE   PSEUD    0    0.0000    1.3124   -4.8971   -1.2330    0    0    0    0   13
	13 QR   PSEUD    0    0.0000    1.8414   -3.7873   -1.2373    0    0    0    0    0
	14 CG   C_VIN    0    0.0000    2.4445   -2.5201   -1.2433    7   15   22    0    0
	15 CD1  C_ARO    0    0.0000    2.5614   -3.3621   -0.1441   14   16   17    0    0
	16 HD1  H_ARO    0    0.0000    3.1267   -3.0292    0.7144   15    0    0    0   11
	17 CE1  C_ARO    0    0.0000    1.9698   -4.6105   -0.1345   15   18   19    0    0
	18 HE1  H_ARO    0    0.0000    2.0716   -5.2509    0.7293   17    0    0    0   12
	19 CZ   C_VIN    0    0.0000    1.2489   -5.0319   -1.2323   17   20   24    0    0
	20 CE2  C_ARO    0    0.0000    1.1182   -4.2145   -2.3356   19   21   22    0    0
	21 HE2  H_ARO    0    0.0000    0.5541   -4.5449   -3.1953   20    0    0    0   12
	22 CD2  C_ARO    0    0.0000    1.7144   -2.9684   -2.3370   14   20   23    0    0
	23 HD2  H_ARO    0    0.0000    1.6144   -2.3257   -3.1998   22    0    0    0   11
	24 OH   O_HYD    0    0.0000    0.6577   -6.2745   -1.2276   19   25    0    0    0
	25 P    P_ALI    0    0.0000    0.7301   -7.3207   -0.0192   24   26   27   28    0
	26 O1P  O_BYL    0    0.0000    1.5265   -6.7553    1.0928   25    0    0    0    0
	27 O2P  O_BYL    0    0.0000   -0.6332   -7.6225    0.4714   25    0    0    0    0
	28 O3P  O_BYL    0    0.0000    1.3657   -8.5772   -0.4748   25    0    0    0    0
	29 C    C_BYL    0    0.0000    1.1680   -0.1054   -2.4506    5   30   31    0    0
	30 O    O_BYL    0    0.0000    1.6224   -0.3027   -3.5762   29    0    0    0    0
	31 N    N_AMI    0    0.0000   -0.1313    0.0279   -2.2068   29    0    0    0    0


CSTABLE      65
	 1 HIST N     6491  119.62    4.02  105.00  136.48
	 2 HIST H     7180    8.24    0.69    5.24   12.39
	 3 HIST CA    6056   56.49    2.33   46.01   66.98
	 4 HIST HA    5556    4.61    0.44    1.93    8.90
	 5 HIST CB    5614   30.22    2.10   18.75   43.30
	 6 HIST QB   10044    3.08    0.38    0.18    8.70
	 7 HIST CG      78  131.48    3.31  122.67  137.19
	 8 HIST ND1    127  195.73   18.04  164.64  229.14
	 9 HIST CD2   1972  120.46    3.45  112.07  159.95
	10 HIST CE1   1530  137.57    2.39  104.67  144.54
	11 HIST HD2   3579    7.03    0.44    3.46    9.01
	12 HIST HE1   2960    7.97    0.51    3.21   10.26
	13 HIST HE2    126    9.76    2.61    6.58   16.53
	14 HIST C     4297  175.28    1.98  166.90  182.80
	15 HIS+ N     6491  119.62    4.02  105.00  136.48
	16 HIS+ H     7180    8.24    0.69    5.24   12.39
	17 HIS+ CA    6056   56.49    2.33   46.01   66.98
	18 HIS+ HA    5556    4.61    0.44    1.93    8.90
	19 HIS+ CB    5614   30.22    2.10   18.75   43.30
	20 HIS+ QB   10044    3.08    0.38    0.18    8.70
	21 HIS+ CG      78  131.48    3.31  122.67  137.19
	22 HIS+ ND1    127  195.73   18.04  164.64  229.14
	23 HIS+ CD2   1972  120.46    3.45  112.07  159.95
	24 HIS+ HD1    265    8.78    2.70    3.79   17.20
	25 HIS+ CE1   1530  137.57    2.39  104.67  144.54
	26 HIS+ NE2    146  182.19   14.39  161.70  226.29
	27 HIS+ HD2   3579    7.03    0.44    3.46    9.01
	28 HIS+ HE1   2960    7.97    0.51    3.21   10.26
	29 HIS+ HE2    126    9.76    2.61    6.58   16.53
	30 HIS+ C     4297  175.28    1.98  166.90  182.80
	31 PROU N       26  134.49    7.01  106.00  142.10
	32 PROU CD     103   50.06    1.02   46.50   52.60
	33 PROU CA     117   62.85    1.37   59.50   66.30
	34 PROU HA     116    4.47    0.44    2.27    5.37
	35 PROU CB      95   31.74    1.54   25.90   36.50
	36 PROU QB     214    2.01    0.44   -0.15    2.91
	37 PROU CG      85   26.77    0.91   24.10   28.60
	38 PROU QG     197    1.89    0.39   -0.29    2.55
	39 PROU QD     206    3.69    0.38    1.67    4.61
	40 PROU C       46  175.91    1.79  172.10  180.60
	41 PROO N       26  134.49    7.01  106.00  142.10
	42 PROO CD     103   50.06    1.02   46.50   52.60
	43 PROO QD     206    3.69    0.38    1.67    4.61
	44 PROO CA     117   62.85    1.37   59.50   66.30
	45 PROO HA     116    4.47    0.44    2.27    5.37
	46 PROO CB      95   31.74    1.54   25.90   36.50
	47 PROO QB     214    2.01    0.44   -0.15    2.91
	48 PROO CG      85   26.77    0.91   24.10   28.60
	49 PROO QG     197    1.89    0.39   -0.29    2.55
	50 PROO C       46  175.91    1.79  172.10  180.60
	51 PTR  N       72  121.40    4.29  113.30  130.60
	52 PTR  H       86    8.44    0.92    6.36   10.55
	53 PTR  CA      93   57.07    2.39   50.80   62.17
	54 PTR  HA      91    4.80    0.57    3.79    6.70
	55 PTR  CB      80   38.75    2.30   33.50   45.00
	56 PTR  QB     176    2.96    0.36    1.62    4.10
	57 PTR  QD     170    7.01    0.26    6.28    7.77
	58 PTR  QE     162    6.72    0.22    6.04    7.31
	59 PTR  CG      27  129.40    1.56  125.70  132.50
	60 PTR  CD1    125  132.32    1.25  129.70  137.70
	61 PTR  CE1    122  117.33    1.11  114.30  119.80
	62 PTR  CZ      28  156.22    1.69  150.20  158.70
	63 PTR  CE2    122  117.33    1.11  114.30  119.80
	64 PTR  CD2    125  132.32    1.25  129.70  137.70
	65 PTR  C       41  175.16    1.82  170.20  178.50
'''
special.write(special_text)
special.close()
