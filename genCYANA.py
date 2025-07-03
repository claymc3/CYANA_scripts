## 03/12/2020 Mary Clay PhD
## Updated to python 3 06/27/2022 Mary Clay PhD
## Script for preparing cyana input files 
## name.seq
## dihed.aco
## hbond.upl/hbond.lol
## name.upl


import pandas as pd 
import numpy as np
import os
import sys
import glob
import re
from datetime import datetime

# datetime object containing current date and time
now = datetime.now()
MasterDict = {}
# dd/mm/YY H:M:S
dt_string = now.strftime("%Y-%m-%d %H:%M")
pdb_columns = ['name', 'resn', 'resid', 'X', 'Y', 'Z','nuc']
# Read in PDB file one line at a time, if the first four letter ar ATOM or HETA then it will parse the data into the 
# data frame, using the atome index int he PDB as the row index in the data frame. 
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "CYSS":"C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H","HIST": "H","HISE": "H","HIS+": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "PROO":"P","PROU":"P","CPRO":"P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S',"ADE":"A","RADE":"A","CYT":"C","RCYT":"C","GUA":"G","RGUA":"G","THY":"T","URA":"U"}
A_dict = {'A': 'ALA', 'C': 'CYS', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE', 'H': 'HIS', 'K': 'LYS', 'M': 'MET', 'L': 'LEU', 'N': 'ASN', 'Q': 'GLN', 'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR' }
Atoms_dict = {'I':['CD1'], 'L':['CD1','CD2'], 'V':['CG1','CG2'], 'M':['CE'], 'A':['CB'], 'T':['CG2'], 'W':['NE1','HE1'], 'F':['CE1','CE2','HE1','HE2'], 'Y':['CE1','CE2','HE1','HE2'],'K':['CA','CB','CG','CD','CE','HA','HB3','HG3','HD3','HE2'],'R':['CA','CB','CG','CD','HA','HB3','HG3','HD3'],'H':['CE1','HE1','CD2','HD2']}


cwd = os.getcwd()
#####
if len(sys.argv) == 1:
	print('''

Usage: 
	gencya [pdb] [chain] [Labeling] [residues] [TALOS]

Required Input:

	PDB 			Reference PDB used to generate initial upl file. If this is 
					not located in current directory provide path. 

	Chain:			Which chain in PDB should be used

	Labeling		What side chains are labeled and can be used to generate 
					heavy atom – heavy atom model based upls 
					Supported labeling of: I, L, V, M, A,T, Y, F, or W
					I = I-CD1, L = L-CD1, L-CD2, V = V-CG1, V-CG2, M = M-CE,
					A = A-CB, T = T-CG2 F = F-CE1, F-CE2, Y = Y-CE1, Y-CE2,
					W = W-NE1, K, R, H

	residues 		Residues to use in upl generation and rmsd calumniation.
					20-200,300-490
					Residues of dynamic regions, and regions that do or are 
					expected to show variations between states should be excluded. 

	TALOS 			Path to results of TALOSN analysis. The pred.tab will be 
					used to generate the dihid.aco input file. The predss.tab
					will be used to filter out unstructured regions from the 
					upl, hbond.lol, and hbond.upl files 

Assuming you have saved your .proj.seq in the current location
Assuming you have saved your .prot and .peak files are in current location

name is taken from the first prot file found in this directory.
All peaks and prot finles found in this directory will be listed in the CALC.cya 
All upl, lol, and aco files created or already in this directory will be listed 
in CALC.cya

OutPut:
	CALC.cya
	init.cya
	name.seq
	name_modle.upl
	hbond.upl
	hbond.lol
	dihed.aco
	inital.aco
	gencya_log
''')
	exit()
Vnum = '2.1'
szFileName = glob.glob('*.seq')[0]
print(szFileName)
in_pdb  = sys.argv[1]
chain = sys.argv[2]
atoms = sys.argv[3]
residues = sys.argv[4]
TALOSdir = sys.argv[5]
prot = glob.glob('*.prot')[0]
name = prot.split('.')[0]
if szFileName == name + '.seq': pass 
else: os.rename(szFileName, name + '.seq')

talosSS = os.path.join(TALOSdir +'/predSS.tab')
talosS2 = os.path.join(TALOSdir +'/predS2.tab')
indihed = os.path.join(TALOSdir +'/pred.tab')
inchi1 = os.path.join(TALOSdir +'/predChi1.tab')
print("Generated using V2.1")
print("Using " + talosSS)
print("Using " + indihed)
# print("Using " + talosS2)
log = open('gencya_log', 'a')
log.write('## Generated using genCYANA_{:} on {:} \n'.format(Vnum,dt_string))
log.write('gencya {:} {:} {:} {:} {:}\n'.format(in_pdb,chain,atoms,residues,TALOSdir))

#------------------------------------------------------------------------------
# Read the prot files and figure what assignments are available
# 
Assignments = []
prots = glob.glob('*.prot')
for prot in prots:
	for line in open(cwd + '/' + prot.replace('-final.prot','.prot')).readlines():
		if '#' != line[0] and line.strip():
			resi = line.strip().split()[4]
			atom = line.strip().split()[3]
			assign = resi + '-' + atom
			if assign not in Assignments:
				Assignments.append(assign)

#------------------------------------------------------------------------------
# Read in sequence and generate name.seq file and translate the number to correct
# amino acid abrivation

num2AAA = {}
num_seq = []
for line in open(name + '.seq').readlines():
		num_seq.append((AAA_dict[line.strip().split()[0]], line.strip().split()[1]))
		num2AAA[line.strip().split()[1]] = line.strip().split()[0]
#
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Read in predSS.tab from TALOS run and great dictionary indicating secondary
# structure of residue to use for preparation of hbond and manual upl/lol 
#
talos_lines = [line.strip() for line in open(talosSS).readlines() if line.strip() and not re.search('[A-Z]', line[0])]
SecStrDict = {}
for line in talos_lines:
	res = line.split()[1] + line.split()[0]
	SecStrDict[res] = line.split()[-1].upper().replace('C','L')

#------------------------------------------------------------------------------
# Read in predSS.tab from TALOS run and great dictionary indicating secondary
# structure of residue to use for preparation of hbond and manual upl/lol 
#
talosS2_lines = [line.strip() for line in open(talosS2).readlines() if line.strip() and not re.search('[A-Z]', line[0])]
S2Dict = {}
for line in talosS2_lines:
	res = line.split()[1] + line.split()[0]
	S2Dict[res] = float(line.split()[-1].replace('9999.0000','0.600'))

#------------------------------------------------------------------------------
# Read in TALOS pred.tab and create dihed.aco file 

# dihed = glob.glob("*.aco")[0]
dihed_lines = [line.strip() for line in open(indihed).readlines() if line.strip() and not re.search('[A-Z]', line[0])]
chi1_lines = [line.strip() for line in open(inchi1).readlines() if line.strip() and not re.search('[A-Z]', line[0])]
# (resid, slc, phi, psi, dphi, dpsi, s2, count, cs_cound, Class)
aco = open('dihed.aco','w')
# aco.write('## Chi1 values from TALOS\n')
# for line in chi1_lines:
# 	cline = line.strip().split()
# 	if cline[6] != 'na':
# 		dchi = float(cline[-1])
# 		if dchi<10: dchi = 10.0
# 		if dchi>35: dchi = 35
# 		chi1 = float(cline[7])
# 		if chi1 < 0: chi1 = chi1 + 360.0
# 		aco.write("#  " + line + "\n")
# 		aco.write("{:>4s}  {:<4s} CHI1  {:8.1f}{:8.1f}\n\n".format(cline[0], A_dict[cline[1]], chi1-dchi, chi1+dchi))
aco.write('## Generated using PrepCyana_{:}\n## Dihedrals extracted from {:}\n'.format(Vnum,indihed))

aco.write('\n## Preformatted entries for chi1/chi2 restraints that should be included as strucutre is refined\n')
for (res, resn) in num_seq:
	if res == 'I':
		aco.write('## {:>5}  ILE  CHI1    360.0   360.0\n## {:>5}  ILE  CHI21   360.0   360.0\n'.format(resn,resn))
aco.write('\n\n')
for (res, resn) in num_seq:
	if res == 'L':
		aco.write('## {:>5}  LEU  CHI1    360.0   360.0\n## {:>5}  LEU  CHI2    360.0   360.0\n'.format(resn,resn))
aco.write('\n\n')
for (res, resn) in num_seq:
	if res == 'M':
		aco.write('## {:>5}  MET  CHI1    360.0   360.0\n## {:>5}  MET  CHI2    360.0   360.0\n'.format(resn,resn))
aco.write('\n\n')
for (res, resn) in num_seq:
	if res == 'F':
		aco.write('## {:>5}  PHE  CHI1    360.0   360.0\n## {:>5}  PHE  CHI2    360.0   360.0\n'.format(resn,resn))
aco.write('\n\n')
for (res, resn) in num_seq:
	if res == 'Y':
		aco.write('## {:>5}  TYR  CHI1    360.0   360.0\n## {:>5}  TYR  CHI2    360.0   360.0\n'.format(resn,resn))
aco.write('\n\n')
for (res, resn) in num_seq:
	if res == 'W':
		aco.write('## {:>5}  TRP  CHI1    360.0   360.0\n## {:>5}  TRP  CHI2    360.0   360.0\n'.format(resn,resn))
aco.write('\n\n')
if 'K' in sys.argv[3]:
	for (res, resn) in num_seq:
		if res == 'K':
			aco.write('## {:>5}  LYS  CHI1    360.0   360.0\n## {:>5}  LYS  CHI2    360.0   360.0\n'.format(resn,resn))
	aco.write('\n\n')
if 'R' in sys.argv[3]:
	for (res, resn) in num_seq:
		if res == 'R':
			aco.write('## {:>5}  ARG  CHI1    360.0   360.0\n## {:>5}  ARG  CHI2    360.0   360.0\n'.format(resn,resn))
	aco.write('\n\n')

aco.write('### Other rotamers to check ###')
for (res, resn) in num_seq:
	if res == 'H':
		aco.write('## {:>5}  {:<4} CHI1    360.0   360.0\n## {:>5}  {:<4} CHI2    360.0   360.0\n'.format(resn,num2AAA[resn],resn,num2AAA[resn]))
aco.write('\n\n')
for (res, resn) in num_seq:
	if res == 'D':
		aco.write('## {:>5}  ASP  CHI1    360.0   360.0\n## {:>5}  ASP  CHI2    360.0   360.0\n'.format(resn,resn))
aco.write('\n\n')
for (res, resn) in num_seq:
	if res == 'N':
		aco.write('## {:>5}  ASN  CHI1    360.0   360.0\n## {:>5}  ASN  CHI2    360.0   360.0\n'.format(resn,resn))
aco.write('\n\n')
for (res, resn) in num_seq:
	if res == 'E':
		aco.write('## {:>5}  GLU  CHI1    360.0   360.0\n## {:>5}  GLU  CHI2    360.0   360.0\n'.format(resn,resn))
aco.write('\n\n')
for (res, resn) in num_seq:
	if res == 'Q':
		aco.write('## {:>5}  GLN  CHI1    360.0   360.0\n## {:>5}  GLN  CHI2    360.0   360.0\n'.format(resn,resn))
aco.write('\n\n')


aco.write('## Phi/Psi values from TALOS\n')
phicount, psicount = 0,0 
for line in dihed_lines:
	dline = line.split()
	if dline[-1] in ['Strong','Generous']:
		dphi = float(dline[4])
		dphi = 20.0
		if dphi<10: dphi = 20.0
		if dphi>35: dphi = 30.0
		dpsi = float(dline[5])
		dpsi = 20.0
		if dpsi<10: dpsi = 20.0
		if dpsi>35: dpsi = 30.0
		aco.write("#  " + line + "\n")
		if dline[1] != 'P':
			phicount+=1
			#print "%4d  %4s  PHI  %8.1f%8.1f" % (int(dline[0]), A_dict[dline[1]], float(dline[2])-scale[dline[-1]]*dphi, float(dline[2])+ scale[dline[-1]]*dphi)
			aco.write("{:>5}  {:<4}  PHI  {:8.1f}{:8.1f}\n".format(int(dline[0]), num2AAA[dline[0]], float(dline[2])-dphi, float(dline[2])+dphi))
		#print "%4d  %4s  PSI  %8.1f%8.1f\n" % (int(dline[0]), A_dict[dline[1]], float(dline[3])-scale[dline[-1]]*dpsi, float(dline[2])+scale[dline[-1]]*dpsi)
		aco.write("{:>5}  {:<4}  PSI  {:8.1f}{:8.1f}\n\n".format(int(dline[0]), num2AAA[dline[0]], float(dline[3])-dpsi, float(dline[3])+dpsi))
		psicount+= 1
for line in dihed_lines:
	dline = line.split()
	if dline[-1] == 'Dyn':
		dphi = 30.0
		dpsi = 30.0
		aco.write("#  " + line + "\n")
		if dline[1] != 'P':
			aco.write("#{:>5}  {:<4}  PHI  {:8.1f}{:8.1f}\n".format(int(dline[0]), num2AAA[dline[0]], float(dline[2])-dphi, float(dline[2])+dphi))
		aco.write("#{:>5}  {:<4}  PSI  {:8.1f}{:8.1f}\n\n".format(int(dline[0]), num2AAA[dline[0]], float(dline[3])-dpsi, float(dline[3])+dpsi))

aco.close()
log.write('Extracted {:} PHI angles and {:} PSI angles form TALOS\n'.format(phicount,psicount))



#------------------------------------------------------------------------------
# Prep and inital aco file to dias LEU, ILE, and MET residues to favorable 
# chi1/chi2 regions 

aco2 = open('inital.aco','w')
for (res, resn) in num_seq:
	if res == 'L':
		aco2.write(' {:>5}  LEU   CHI1    155.0   212.0 9.00E-01 type=2\n {:>5}  LEU   CHI2     40.0    88.0 9.00E-01 type=2\n {:>5}  LEU   CHI1    264.0   324.0 9.00E-01 type=2 OR\n {:>5}  LEU   CHI2    143.0   205.0 9.00E-01 type=2\n\n'.format(resn,resn,resn,resn))
	if res == 'I':
		aco2.write(' {:>5}  ILE   CHI1     40.0    80.0 9.00E-01 type=2\n {:>5}  ILE   CHI21   145.0   195.0 9.00E-01 type=2\n {:>5}  ILE   CHI1    168.0   212.0 9.00E-01 type=2 OR\n {:>5}  ILE   CHI21   145.0   185.0 9.00E-01 type=2\n {:>5}  ILE   CHI1    275.0   325.0 9.00E-01 type=2 OR\n {:>5}  ILE   CHI21   138.0   198.0 9.00E-01 type=2\n {:>5}  ILE   CHI1    278.0   325.0 9.00E-01 type=2 OR\n {:>5}  ILE   CHI21   274.0   322.0 9.00E-01 type=2\n {:>5}  ILE   CHI1    172.0   212.0 9.00E-01 type=2 OR\n {:>5}  ILE   CHI21    48.0    82.0 9.00E-01 type=2\n\n'.format(resn,resn,resn,resn,resn,resn,resn,resn,resn,resn))
	if res == 'M':
		aco2.write(' {:>5}  MET   CHI1     45.0    84.0 9.00E-01 type=2\n {:>5}  MET   CHI2    157.0   210.0 9.00E-01 type=2\n {:>5}  MET   CHI1    157.0   208.0 9.00E-01 type=2 OR\n {:>5}  MET   CHI2    150.0   208.0 9.00E-01 type=2\n {:>5}  MET   CHI1    265.0   322.0 9.00E-01 type=2 OR\n {:>5}  MET   CHI2    146.0   212.0 9.00E-01 type=2\n {:>5}  MET   CHI1    266.0   320.0 9.00E-01 type=2 OR\n {:>5}  MET   CHI2    270.0   325.0 9.00E-01 type=2\n {:>5}  MET   CHI1    161.0   208.0 9.00E-01 type=2 OR\n {:>5}  MET   CHI2     42.0    86.0 9.00E-01 type=2\n\n'.format(resn,resn,resn,resn,resn,resn,resn,resn,resn,resn))
aco2.close()
ssa = open('ssa.cya','w')
ssa.close()
#------------------------------------------------------------------------------
# Fetch the .peak files and create a list 
#
peaks = glob.glob(os.path.join(cwd + '/*.peaks'))
peaks_out = ''
for peak in peaks:
	peaks_out = peaks_out + peak.split('/')[-1] + ','
peaks_out=peaks_out[:-1]

prots = glob.glob(os.path.join(cwd + '/*.prot'))
prots_out = ''
for prot in prots:
	prots_out = prots_out + prot.split('/')[-1] + ','
prots_out=prots_out[:-1]

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
allowed_atoms =[["ALA-O", "ALA-N", "ALA-C", "ARG-O", "ARG-N", "ARG-C", "ASN-O", "ASN-N", "ASN-C", "ASP-O", "ASP-N", "ASP-C", "CYS-O", "CYS-N", "CYS-C", "GLU-O", "GLU-N", "GLU-C", "GLN-O", "GLN-N", "GLN-C", "GLY-O", "GLY-N", "GLY-C", "HIS-O", "HIS-N", "HIS-C", "ILE-O", "ILE-N", "ILE-C", "LEU-O", "LEU-N", "LEU-C", "LYS-O", "LYS-N", "LYS-C", "MET-O", "MET-N", "MET-C", "PHE-O", "PHE-N", "PHE-C", "PRO-O", "PRO-N", "SER-O", "SER-N", "SER-C", "SEP-O", "SEP-N", "SEP-C", "TPO-O", "TPO-N", "TPO-C", "THR-O", "THR-N", "THR-C", "TRP-O", "TRP-N", "TRP-C", "TYR-O", "TYR-N", "TYR-C", "VAL-O", "VAL-N", "VAL-C"]
for i in range(len(atoms)):
	for x in range(len(Atoms_dict[atoms[i]])):
		allowed_atoms.append(A_dict[atoms[i]] + '-' + Atoms_dict[atoms[i]][x])

#------------------------------------------------------------------------------
# Generating pandas data frame from PDB file using only specified residues and atoms
#
PDB_df = pd.DataFrame(columns=pdb_columns)
with open(in_pdb) as In_pdb:
	for line in In_pdb:
		if line[0:4] == "ATOM" or line[0:4] == 'HETA':
			chainid = line[21]
			if chainid  == ' ':chainid = line[72]
			if int(line[22:26].strip()) in allowed_resi and chainid == chain:
				group = line[17:20].strip() + '-' + line[12:16].strip()
				if group in allowed_atoms:
					index =  AAA_dict[line[17:20].strip()] + line[22:26].strip() + "-" + line[12:16].strip()
					PDB_df.loc[index, 'name'] = line[12:16].strip()
					PDB_df.loc[index, 'resn2'] = line[16:20].strip()
					PDB_df.loc[index, 'resn'] = num2AAA[line[22:26].strip()]
					PDB_df.loc[index, 'resid'] = line[22:26].strip()
					PDB_df.loc[index, 'X'] = float(line[30:38])
					PDB_df.loc[index, 'Y'] = float(line[38:46])
					PDB_df.loc[index, 'Z'] = float(line[46:54])
					PDB_df.loc[index, 'nuc'] = line[12:16].strip()[0]
					PDB_df.loc[index, 'chain'] = line[21]
					if AAA_dict[line[17:20].strip()] + line[22:26].strip() in SecStrDict.keys():
						PDB_df.loc[index, 'SecStr'] = SecStrDict[AAA_dict[line[17:20].strip()] + line[22:26].strip()]
					if AAA_dict[line[17:20].strip()] + line[22:26].strip() in S2Dict.keys():
						PDB_df.loc[index, 'S2'] = S2Dict[AAA_dict[line[17:20].strip()] + line[22:26].strip()]

#------------------------------------------------------------------------------
# Find Hydrogen bonding connections 
#
constrained = []
hblol = open('hbond.lol','w')
hbupl = open('hbond.upl', 'w')
### For Helical residues from TALOS
helical_hb, sheet_hb = 0,0
hbond_N = PDB_df[(PDB_df['nuc'] == 'N') & (PDB_df['SecStr'] == 'H')].index.tolist()
hbond_O = PDB_df[(PDB_df['nuc'] == 'O') & (PDB_df['SecStr'] == 'H')].index.tolist()
hblol.write('##Generated using PrepCyana_{:}\n## Distances extracted from {:}\n'.format(Vnum,in_pdb))
hbupl.write('##Generated using PrepCyana_{:}\n## Distances extracted from {:}\n'.format(Vnum,in_pdb))
hblol.write("## Helical residues from TALOS \n")
hbupl.write("## Helical residues from TALOS \n")
for O in hbond_O:
	for N in hbond_N:
		if abs(int(PDB_df.loc[O,'resid']) - int(PDB_df.loc[N,'resid'])) == 4: 
			dist = round(np.sqrt(((PDB_df.loc[O,'X'] - PDB_df.loc[N,'X'])**2) + ((PDB_df.loc[O,'Y'] - PDB_df.loc[N,'Y'])**2) + ((PDB_df.loc[O,'Z'] - PDB_df.loc[N,'Z'])**2)),1)
			if dist >= 2.7 and dist <= 3.3: 
				helical_hb+= 1
				constrained.append(str(PDB_df.loc[O,'resid']) + '-' + str(PDB_df.loc[N,'resid']))
				hbupl.write("{:>5}  {:<4}  {:<4}  {:>5}  {:<4}  {:<4}     3.10\n".format(PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
				hbupl.write("{:>5}  {:<4}  {:<4}  {:>5}  {:<4}  {:<4}     2.10\n".format(PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))
				hblol.write("{:>5}  {:<4}  {:<4}  {:>5}  {:<4}  {:<4}     2.70\n".format(PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
				hblol.write("{:>5}  {:<4}  {:<4}  {:>5}  {:<4}  {:<4}     1.80\n".format(PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))
### For Betta Sheet residues from TALOS
hblol.write("## Betta Sheet residues from TALOS \n")
hbupl.write("## Betta Sheet residues from TALOS \n")				
hbond_N = PDB_df[(PDB_df['nuc'] == 'N') & (PDB_df['SecStr'] == 'E')].index.tolist()
hbond_O = PDB_df[(PDB_df['nuc'] == 'O') & (PDB_df['SecStr'] == 'E')].index.tolist()
for O in hbond_O:
	for N in hbond_N:
		if abs(int(PDB_df.loc[O,'resid']) - int(PDB_df.loc[N,'resid'])) > 3: 
			dist = round(np.sqrt(((PDB_df.loc[O,'X'] - PDB_df.loc[N,'X'])**2) + ((PDB_df.loc[O,'Y'] - PDB_df.loc[N,'Y'])**2) + ((PDB_df.loc[O,'Z'] - PDB_df.loc[N,'Z'])**2)),1)
			if dist >= 2.7 and dist <= 3.3: 
				sheet_hb+= 1
				constrained.append(str(PDB_df.loc[O,'resid']) + '-' + str(PDB_df.loc[N,'resid']))
				hbupl.write("{:>5}  {:<4}  {:<4}  {:>5}  {:<4}  {:<4}     3.10\n".format(PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
				hbupl.write("{:>5}  {:<4}  {:<4}  {:>5}  {:<4}  {:<4}     2.10\n".format(PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))
				hblol.write("{:>5}  {:<4}  {:<4}  {:>5}  {:<4}  {:<4}     2.70\n".format(PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'], PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name']))
				hblol.write("{:>5}  {:<4}  {:<4}  {:>5}  {:<4}  {:<4}     1.80\n".format(PDB_df.loc[O,'resid'],PDB_df.loc[O,'resn'],PDB_df.loc[O,'name'].replace('N','H'), PDB_df.loc[N,'resid'],PDB_df.loc[N,'resn'],PDB_df.loc[N,'name'].replace('N','H')))
hblol.close()
hbupl.close()
hblol.close()
hbupl.close()
log.write('Generated {:} unique Helical hbonds\nGenerated {:} unique Beta sheet hbonds\n'.format(helical_hb,sheet_hb))
print("Generated hbond.upl and hbond.lol")
#------------------------------------------------------------------------------
# Prepare upl constraints for initial calculation only considering structured elements
# based on TALOS predSS.tab
#
NN_used, NN_lines = [],[]
upl = open('{:}_{:}.upl'.format(in_pdb.replace('.pdb',''),chain),'w')
upl.write('##Generated using PrepCyana_{:}\n## Distances extracted from {:}\n'.format(Vnum,in_pdb))
N_list = PDB_df[(PDB_df['nuc'] == 'N')  & (PDB_df['resn'] != 'PRO')].index.tolist()
C_list = PDB_df[(PDB_df['nuc'] == 'C')].index.tolist()
for i in range(len(N_list)):
	for x in range(len(N_list)):
		diff = abs(int(PDB_df.loc[N_list[i],'resid']) - int(PDB_df.loc[N_list[x],'resid']))
		if PDB_df.loc[N_list[i],'SecStr'] == 'H' or PDB_df.loc[N_list[x],'SecStr'] == 'H': ridiff = 3
		else: ridiff = 3
		ridiff = 4
		if diff >= ridiff:
			dist = np.round(np.sqrt(((PDB_df.loc[N_list[i],'X'] - PDB_df.loc[N_list[x],'X'])**2) + ((PDB_df.loc[N_list[i],'Y'] - PDB_df.loc[N_list[x],'Y'])**2) + ((PDB_df.loc[N_list[i],'Z'] - PDB_df.loc[N_list[x],'Z'])**2)),1)
			if dist < 3.1: 
				print('N-N distance shorter than 3.1, likely van der Waals clash refine PDB model and try again')
				#exit()
			if dist < 5.0:
				constraint = N_list[x] + '-' + N_list[i]
				atom1 = str(PDB_df.loc[N_list[i],'resid']) + '-' + PDB_df.loc[N_list[i],'name']
				atom2 = str(PDB_df.loc[N_list[x],'resid']) + '-' + PDB_df.loc[N_list[x],'name']
				if constraint not in NN_used:
					NN_used.append(N_list[i] + '-' + N_list[x])
					NN_out = "{:>5}  {:<4}  {:<4}   {:>5}  {:<4}  {:<4}   {:6.2f}\n".format(PDB_df.loc[N_list[i],'resid'],PDB_df.loc[N_list[i],'resn'],PDB_df.loc[N_list[i],'name'], PDB_df.loc[N_list[x],'resid'],PDB_df.loc[N_list[x],'resn'],PDB_df.loc[N_list[x],'name'],dist)
					if atom1 in Assignments and atom2 in Assignments: pass
					if atom1 not in Assignments:
						NN_out = NN_out.replace('\n',' # missing {:}\n'.format(N_list[i]))
					if atom2 not in Assignments:
						NN_out = NN_out.replace('\n',' # missing {:}\n'.format(N_list[x]))
					if 'missing' in NN_out:
						NN_out = '#' + NN_out
					NN_lines.append(NN_out)

print("Made {:3.0f} NN upl constraints ".format(len(NN_used)))
log.write("Generated  {:} NN upl constraints\n".format(len(NN_used)))
NC_used, NC_lines, NC_ids = [], [], []
for i in range(len(N_list)):
	for x in range(len(C_list)):
		diff = abs(int(PDB_df.loc[N_list[i],'resid']) - int(PDB_df.loc[C_list[x],'resid']))
		if diff >= 3:
			dist = np.round(np.sqrt(((PDB_df.loc[N_list[i],'X'] - PDB_df.loc[C_list[x],'X'])**2) + ((PDB_df.loc[N_list[i],'Y'] - PDB_df.loc[C_list[x],'Y'])**2) + ((PDB_df.loc[N_list[i],'Z'] - PDB_df.loc[C_list[x],'Z'])**2)),1)
			if dist < 3.1:
				NC_out = "{:>5}  {:<4}  {:<4}   {:>5}  {:<4}  {:<4}   {:6.2f}\n".format(PDB_df.loc[N_list[i],'resid'],PDB_df.loc[N_list[i],'resn'],PDB_df.loc[N_list[i],'name'], PDB_df.loc[C_list[x],'resid'],PDB_df.loc[C_list[x],'resn'],PDB_df.loc[C_list[x],'name'],dist)
				print('{:} distance shorter than 3.1, likely van der Waals clash refine PDB model and try again'.format(NC_out))
				#exit()
			if dist < 6.0:
				atom1 = str(PDB_df.loc[N_list[i],'resid']) + '-' + PDB_df.loc[N_list[i],'name']
				atom2 = str(PDB_df.loc[C_list[x],'resid']) + '-' + PDB_df.loc[C_list[x],'name']
				NC_used.append(N_list[i] + '-' + C_list[x])
				NC_out = "{:>5}  {:<4}  {:<4}   {:>5}  {:<4}  {:<4}   {:6.2f}\n".format(PDB_df.loc[N_list[i],'resid'],PDB_df.loc[N_list[i],'resn'],PDB_df.loc[N_list[i],'name'], PDB_df.loc[C_list[x],'resid'],PDB_df.loc[C_list[x],'resn'],PDB_df.loc[C_list[x],'name'],dist)
				if atom1 in Assignments and atom2 in Assignments: pass
				if atom1 not in Assignments:
					NC_out = NC_out.replace('\n',' # missing {:}\n'.format(N_list[i]))
				if atom2 not in Assignments:
					NC_out = NC_out.replace('\n',' # missing {:}\n'.format(C_list[x]))
				NC_lines.append(NC_out)
				if (str(PDB_df.loc[N_list[i],'resid']),str(PDB_df.loc[C_list[x],'resid'])) not in NC_ids:
					NC_ids.append((str(PDB_df.loc[N_list[i],'resid']),str(PDB_df.loc[C_list[x],'resid'])))

CC_used,CC_lines,CC_ids = [],[],[]
for i in range(len(C_list)):
	for x in range(len(C_list)):
		diff = abs(int(PDB_df.loc[C_list[i],'resid']) - int(PDB_df.loc[C_list[x],'resid']))
		if diff >= 3:
			dist = np.round(np.sqrt(((PDB_df.loc[C_list[i],'X'] - PDB_df.loc[C_list[x],'X'])**2) + ((PDB_df.loc[C_list[i],'Y'] - PDB_df.loc[C_list[x],'Y'])**2) + ((PDB_df.loc[C_list[i],'Z'] - PDB_df.loc[C_list[x],'Z'])**2)),1)
			if dist < 3.1:
				CC_out = "{:>5}  {:<4}  {:<4}   {:>5}  {:<4}  {:<4}   {:6.2f}\n".format(PDB_df.loc[C_list[i],'resid'],PDB_df.loc[C_list[i],'resn'],PDB_df.loc[C_list[i],'name'], PDB_df.loc[C_list[x],'resid'],PDB_df.loc[C_list[x],'resn'],PDB_df.loc[C_list[x],'name'],dist)
				print('{:} distance shorter than 3.1, likely van der Waals clash refine PDB model and try again'.format(CC_out))
				#exit()
			if dist < 6.0:
				constraint = C_list[x] + '-' + C_list[i]
				atom1 = str(PDB_df.loc[C_list[i],'resid']) + '-' + PDB_df.loc[C_list[i],'name']
				atom2 = str(PDB_df.loc[C_list[x],'resid']) + '-' + PDB_df.loc[C_list[x],'name']
				if constraint not in CC_used:
					CC_used.append(C_list[i] + '-' + C_list[x])
					CC_out = "{:>5}  {:<4}  {:<4}   {:>5}  {:<4}  {:<4}   {:6.2f}\n".format(PDB_df.loc[C_list[i],'resid'],PDB_df.loc[C_list[i],'resn'],PDB_df.loc[C_list[i],'name'], PDB_df.loc[C_list[x],'resid'],PDB_df.loc[C_list[x],'resn'],PDB_df.loc[C_list[x],'name'],dist)
					if atom1 in Assignments and atom2 in Assignments: pass
					if atom1 not in Assignments:
						CC_out = CC_out.replace('\n',' # missing {:}\n'.format(C_list[i]))
					if atom2 not in Assignments:
						CC_out = CC_out.replace('\n',' # missing {:}\n'.format(C_list[x]))
					CC_lines.append(CC_out)
					if (str(PDB_df.loc[C_list[i],'resid']),str(PDB_df.loc[C_list[x],'resid'])) not in CC_ids:
						CC_ids.append((str(PDB_df.loc[C_list[i],'resid']),str(PDB_df.loc[C_list[x],'resid'])))



#------------------------------------------------------------------------------
# Filter the NC and CC distance restraints so that only the shortest connection
# to geminal methyls and Aromatic CE1/CE2 is kept. 
# This means that 
#		712 LEU   CD1   716  LEU   CD1     4.20
#		712 LEU   CD1   716  LEU   CD2     4.29
#		712 LEU   CD2   716  LEU   CD1     4.46
#		712 LEU   CD2   716  LEU   CD2     5.50
# becomes 
#		712 LEU   CD1   716  LEU   CD1     4.20
#------------------------------------------------------------------------------

upl.write('### N-N Distances\n')
upl.writelines(NN_lines)
NC_outlines = []
for (nid,cid) in NC_ids:
	temp,temp2 = [],[]
	for line in NC_lines:
		if line.split()[0] == str(nid) and line.split()[3] == str(cid):
			temp.append(line)
			temp2.append(float(line.split()[6]))
	if len(temp) >= 1:
		NC_outlines.append(temp[temp2.index(min(temp2))])
CC_outlines = []
for (cid1,cid2) in CC_ids:
	temp,temp2 = [],[]
	for line in CC_lines:
		if line.split()[0] == str(cid1) and line.split()[3] == str(cid2):
			temp.append(line)
			temp2.append(float(line.split()[6]))
	if len(temp) >= 1:
		CC_outlines.append(temp[temp2.index(min(temp2))])

NC_methyl, NC_Aro, CC_methyl, CC_Aro,CC_LYS,NC_LYS,CC_ARG,NC_ARG,CC_HIS,NC_HIS,CC_HIST,NC_HIST = [], [], [], [], [], [], [], [], [], [], [], []
NC_sorted,CC_sorted = [],[]
MasterDict['NC_methyl'] = NC_methyl;MasterDict['NC_Aro'] = NC_Aro;MasterDict['CC_methyl'] = CC_methyl;MasterDict['CC_Aro'] = CC_Aro;MasterDict['CC_LYS'] = CC_LYS;MasterDict['NC_LYS'] = NC_LYS;MasterDict['CC_ARG'] = CC_ARG;MasterDict['NC_ARG'] = NC_ARG;MasterDict['CC_HIS'] = CC_HIS;MasterDict['NC_HIS'] = NC_HIS;MasterDict['CC_HIST'] = CC_HIST;MasterDict['NC_HIST'] = NC_HIST

for line in NC_outlines:
	if line.split()[4] in ['PHE','TYR','PTR']:
		if 'missing' in line: line = '#' + line
		NC_Aro.append(line)
		NC_sorted.append(line)
for line in NC_outlines:
	if line.split()[4] in ['LYS','ARG','HIS','HIST']:
		nc_list = MasterDict['NC_{:}'.format(line.split()[4])]
		if 'missing' in line: line = '#' + line
		nc_list.append(line)
		NC_sorted.append(line)
for line in NC_outlines:
	if line not in NC_sorted:
		if 'missing' in line: line = '#' + line
		NC_methyl.append(line)
upl.write('### N-C Methyl Distances\n')
upl.writelines(NC_methyl)
upl.write('### N-C Aromatic Distances\n')
upl.writelines(NC_Aro)
if 'K' in sys.argv[3] and len(NC_LYS) >=1:
	upl.write('### N-C Lysine Distances\n')
	upl.writelines(NC_LYS)
if 'R' in sys.argv[3] and len(NC_ARG) >=1:
	upl.write('### N-C Arginine Distances\n')
	upl.writelines(NC_ARG)
if 'H' in sys.argv[3] and len(NC_HIS) >=1:
	upl.write('### N-C Histadine Distances\n')
	upl.writelines(NC_HIS)
	if len(NC_HIST) >=1:
		upl.writelines(NC_HIST)

for line in CC_outlines:
	if line.split()[1] in ['PHE','TYR','PTR'] and line not in CC_sorted:
		if 'missing' in line: line = '#' + line
		CC_Aro.append(line)
		CC_sorted.append(line)
	if line.split()[4] in ['PHE','TYR','PTR'] and line not in CC_sorted:
		if 'missing' in line: line = '#' + line
		CC_Aro.append(line)
		CC_sorted.append(line)
	if line.split()[1] in ['LYS','ARG','HIS','HIST'] and line not in CC_sorted:
		cc_list = MasterDict['CC_{:}'.format(line.split()[1])]
		if 'missing' in line: line = '#' + line
		cc_list.append(line)
		CC_sorted.append(line)
	if line.split()[4] in ['LYS','ARG','HIS','HIST'] and line not in CC_sorted:
		cc_list = MasterDict['CC_{:}'.format(line.split()[4])]
		if 'missing' in line: line = '#' + line
		cc_list.append(line)
		CC_sorted.append(line)
for line in CC_outlines:
	if line not in CC_sorted:
		if 'missing' in line: 
			line = '#' + line
		CC_methyl.append(line)

upl.write('### C-C Methyl-Methyl Distances\n')
upl.writelines(CC_methyl)
upl.write('### C-C Methyl-Aromatic Distances\n')
upl.writelines(CC_Aro)
if 'K' in sys.argv[3] and len(CC_LYS) >=1:
	upl.write('### C-C Lysine Distances\n')
	upl.writelines(CC_LYS)
if 'R' in sys.argv[3] and  len(CC_ARG) >=1:
	upl.write('### C-C Arginine Distances\n')
	upl.writelines(CC_ARG)
if 'H' in sys.argv[3] and len(CC_HIS) >=1:
	upl.write('### C-C Histadine Distances\n')
	upl.writelines(CC_HIS)
	if len(CC_HIST) >=1:
		upl.writelines(CC_HIST)

upl.write('\n\n')
print("Original %3.0f NC upl constraints " % (len(NC_lines)))
print("Made %3.0f NC upl constraints " % (len(NC_outlines)))
log.write("Generated {:} NC upl constraints\n".format(len(NC_outlines)))
print("Original %3.0f CC upl constraints " % (len(CC_lines)))
print("Made %3.0f CC upl constraints " % (len(CC_outlines)))
log.write("Generated  {:} CC upl constraints\n\n".format(len(CC_outlines)))
upl.close()

CYANA = open('CALC.cya','w')
clac_text = '''
peaks       := {:}      # names of NOESY peak lists
prot        := {:}      # names of chemical shift lists
constraints := {:},{:},{:},{:}      # additional (non-NOE) constraints
tolerance   := 0.01,0.01,0.01,0.01      # chemical shift tolerances
calibration :=      # NOE calibration parameters
dref        := 4.5
structures  := 100,20      # number of initial, final structures
steps       := 20000       # number of torsion angle dynamics steps
rmsdrange   := {:}      # residue range for RMSD calculation
randomseed  := 52541      # random number generator seed

upl_values  := 2.4,7.0

#subroutine KEEP
#   peaks select "*, * number=10000..40000"
#end

ssa
noeassign peaks=$peaks prot=$prot autoaco # keep=KEEP 
'''.format(peaks_out, prots_out, '{:}_{:}.upl'.format(in_pdb.replace('.pdb',''),chain), 'hbond.lol', 'hbond.upl','dihed.aco,inital.aco,manual.upl',residues.replace("-",".."))
CYANA.write(clac_text)
CYANA.close()

manupl = open('manual.upl','w')
manupl.write('### Manualy generated UPL entries based on NOESY Data ###\n')
manupl.write('#{:>5}  {:<4}  {:<4}   {:>5}  {:<4}  {:<4}   {:6.1f}\n'.format('123','ALA','CB','456','ARG','N',6.0))
manupl.write('#{:>5}  {:<4}  {:<4}   {:>5}  {:<4}  {:<4}   {:6.1f}\n'.format('123','ALA','QB','456','ARG','H',6.0))
manupl.write('#{:>5}  {:<4}  {:<4}   {:>5}  {:<4}  {:<4}   {:6.1f}\n'.format('123','PHE','CE1','456','ARG','N',6.0))
manupl.write('#{:>5}  {:<4}  {:<4}   {:>5}  {:<4}  {:<4}   {:6.1f}\n'.format('123','PHE','QE','456','ARG','H',6.0))
manupl.write('#{:>5}  {:<4}  {:<4}   {:>5}  {:<4}  {:<4}   {:6.1f}\n'.format('123','PHE','HE1','456','ARG','H',6.0))
manupl.close()
Init_Cyana = open('init.cya','w')
init_text = '''name:={:}
rmsdrange:={:}
cyanalib
read lib special.lib append
#read lib cyana_Zn2.lib append
nproc:=60
read seq $name.seq
#molecules define 1..146 201..346
#molecule identity
#weight_ide=0.09
#molecule symdist "CA 1..146" "CA 201..346"
#weight_sym=0.025'''.format(name,residues.replace('-','..'))
Init_Cyana.write(init_text)

Init_Cyana.write("###lock stereospecific assignments based on ProS sample assignments\n\n### Leucines\n\n")
for (res, resn) in num_seq:
	if res == 'L':
		Init_Cyana.write('# atom stereo "QD1  QD2  {:>5}"  #{:}{:}\n'.format(resn,res,resn))
Init_Cyana.write('\n\n### Valines\n\n')
for (res, resn) in num_seq:
	if res == 'V':
		Init_Cyana.write('# atom stereo "QG1  QG2  {:>5}"  #{:}{:}\n'.format(resn,res,resn))
if 'K' in sys.argv[3]:
	Init_Cyana.write('\n\n### Sail Lysine\n\n')
	for (res, resn) in num_seq:
		if res == 'K':
			Init_Cyana.write(' atom stereo "HB3  HG3  HD3  HE2  {:>5}"  #{:}{:}\n'.format(resn,res,resn))
if 'R' in sys.argv[3]:
	Init_Cyana.write('\n\n### Sail Arginine\n\n')
	for (res, resn) in num_seq:
		if res == 'R':
			Init_Cyana.write(' atom stereo "HB3  HG2  HD2  {:>5}"  #{:}{:}\n'.format(resn,res,resn))


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

RESIDUE   SEP      6   18    3   17
      1 OMEGA    0    0    0.0000    2    1    3    4    0
      2 PHI      0    0    0.0000    1    3    5   16    0
      3 CHI1     0    0    0.0000    3    5    6    7   14
      4 CHI2     0    0    0.0000    5    6    7    8   11
      5 CHI3     0    0    0.0000    6    7    8    9   11
      6 PSI      0    0    0.0000    3    5   16   18    0
      1 C    C_BYL    0    0.0000    2.6673    0.9171    0.8235    2    3    0    0    0
      2 O    O_BYL    0    0.0000    2.2658    1.5682   -0.1398    1    0    0    0    0
      3 N    N_AMI    0    0.0000    1.8550    0.4210    1.7510    1    4    5    0    0
      4 H    H_AMI    0    0.0000    2.2353   -0.0900    2.4957    3    0    0    0    0
      5 CA   C_ALI    0    0.0000    0.4010    0.6200    1.6870    3    6   15   16    0
      6 CB   C_ALI    0    0.0000   -0.1390    0.0150    0.3910    5    7   12   13    0
      7 OG   O_EST    0    0.0000    0.4770    0.6550   -0.7270    6    8    0    0    0
      8 P    P_ALI    0    0.0000   -0.1350   -0.0270   -2.0500    7    9   10   11    0
      9 O1P  O_BYL    0    0.0000   -1.6010    0.1720   -2.0740    8    0    0    0    0
     10 O2P  O_HYD    0    0.0000    0.5200    0.6490   -3.3560    8    0    0    0    0
     11 O3P  O_HYD    0    0.0000    0.1910   -1.6030   -2.0410    8    0    0    0    0
     12 HB2  H_ALI    0    0.0000    0.0820   -1.0510    0.3670    6    0    0    0   14
     13 HB3  H_ALI    0    0.0000   -1.2180    0.1630    0.3440    6    0    0    0   14
     14 QB   PSEUD    0    0.0000   -0.5680   -0.4440    0.3555    0    0    0    0    0
     15 HA   H_ALI    0    0.0000    0.1790    1.6870    1.7110    5    0    0    0    0
     16 C    C_BYL    0    0.0000   -0.2490   -0.0530    2.8670    5   17   18    0    0
     17 O    O_BYL    0    0.0000    0.4272   -0.6527    3.7012   16    0    0    0    0
     18 N    N_AMI    0    0.0000   -1.5723    0.0440    2.9422   16    0    0    0    0

RESIDUE   TPO      9   23    3   22
   1 OMEGA    0    0    0.0000    2    1    3    4    0
   2 PHI      0    0    0.0000    1    3    5   21    0
   3 CHI1     0    0    0.0000    3    5    6    7   19
   4 CHI2     0    0    0.0000    5    6    7    8   11
   5 CHI3     0    0    0.0000    5    6   12   13   18
   6 CHI4     0    0    0.0000    6   12   13   14   18
   7 CHI5     0    0    0.0000   12   13   15   16   16
   8 CHI6     0    0    0.0000   12   13   17   18   18
   9 PSI      0    0    0.0000    3    5   21   23    0
   1 C    C_BYL    0    0.0000    2.4657   -1.1866    2.5242    2    3    0    0    0
   2 O    O_BYL    0    0.0000    3.2750   -0.3086    2.2291    1    0    0    0    0
   3 N    N_AMI    0    0.0000    1.1530   -1.0400    2.3770    1    4    5    0    0
   4 H    H_AMI    0    0.0000    0.5583   -1.7766    2.6303    3    0    0    0    0
   5 CA   C_ALI    0    0.0000    0.5720    0.1990    1.8440    3    6   20   21    0
   6 CB   C_ALI    0    0.0000    1.1110    0.4490    0.4340    5    7   12   19    0
   7 CG2  C_ALI    0    0.0000    2.6340    0.5800    0.4850    6    8    9   10    0
   8 HG21 H_ALI    0    0.0000    3.0650   -0.3390    0.8810    7    0    0    0   11
   9 HG22 H_ALI    0    0.0000    3.0180    0.7580   -0.5180    7    0    0    0   11
  10 HG23 H_ALI    0    0.0000    2.9060    1.4150    1.1310    7    0    0    0   11
  11 QG2  PSEUD    0    0.0000    2.9963    0.6113    0.4980    0    0    0    0    0
  12 OG1  O_EST    0    0.0000    0.7550   -0.6450   -0.4120    6   13    0    0    0
  13 P    P_ALI    0    0.0000   -0.1420   -0.0390   -1.6030   12   14   15   17    0
  14 O1P  O_BYL    0    0.0000    0.6440    0.9680   -2.3500   13    0    0    0    0
  15 O2P  O_HYD    0    0.0000   -0.5800   -1.2240   -2.6010   13   16    0    0    0
  16 HOP2 H_OXY    0    0.0000   -1.1140   -0.8190   -3.2980   15    0    0    0    0
  17 O3P  O_HYD    0    0.0000   -1.4560    0.6560   -0.9850   13   18    0    0    0
  18 HOP3 H_OXY    0    0.0000   -1.9380   -0.0330   -0.5090   17    0    0    0    0
  19 HB   H_ALI    0    0.0000    0.6800    1.3690    0.0390    6    0    0    0    0
  20 HA   H_ALI    0    0.0000    0.8440    1.0340    2.4900    5    0    0    0    0
  21 C    C_BYL    0    0.0000   -0.9270    0.0700    1.7940    5   22   23    0    0
  22 O    O_BYL    0    0.0000   -1.4832   -0.9616    2.1674   21    0    0    0    0
  23 N    N_AMI    0    0.0000   -1.5914    1.1231    1.3294   21    0    0    0    0

CSTABLE      81
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
	66 SEP  H       98    8.73    0.36    7.18   10.05
	67 SEP  N       73  118.32    1.87  113.00  123.35
	68 SEP  CA      68   57.57    1.56   55.23   63.70
	69 SEP  HA      47    4.50    0.16    4.07    4.83
	70 SEP  CB      68   65.72    1.17   63.29   72.07
	71 SEP  QB      39    4.06    0.19    3.64    4.84
	72 SEP  C       42  173.82    1.15  169.38  176.94 
	73 TPO  H       34    8.87    0.53    7.86   10.23
	74 TPO  N       29  120.39    5.12  103.22  132.72
	75 TPO  CA      27   61.35    1.67   59.38   67.00
	76 TPO  HA      10    4.37    0.32    3.85    5.03
	77 TPO  CB      28   72.70    2.48   69.28   81.93
	78 TPO  HB      11    4.25    0.19    3.91    4.49
	79 TPO  CG2      5   22.28    2.38   20.50   26.46
	80 TPO  QG2      2    1.26    0.14    1.16    1.35
	81 TPO  C       12  173.47    1.06  172.06  175.74
'''
special.write(special_text)
special.close()
import viewCYANA as viewcya

# os.system('python3 precya {:} {:}'.format(in_pdb,TALOSdir))