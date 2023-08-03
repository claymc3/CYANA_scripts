import pandas as pd 
import os
import numpy as np
import sys
from itertools import combinations
from statistics import median,mean



if len(sys.argv)==1:
	print('''

Usage: 
	p2b2noe [pdb] [labeling] [type] [distance cutoff] [chain]

Required Input:

	PDB 			PDB file to generate list of possible NOEs. If this is 
					not located in current directory provide path. 
					Assumes pdb contains only residues of interest.

	labeling 		Specifies use of Methyl or Aromatic Labeling. 
					I, L, V, M, A, T, Y, or F
					I = I-CD1, L = L-CD1, L-CD2, V = V-CG1, V-CG2, M = M-CE, 
					A = A-CB, T = T-CG2 F = F-CE1, F-CE2, Y = Y-CE1, Y-CE2, 

	Type 			Mono methyl (m) or Dymethyl (d)

	Distance cutoff Maximum distance to be listed

	Chain		Specify which chain to use to generate list.
''')
	exit() 

in_pdb  = sys.argv[1]
labeling = sys.argv[2]
labeling_type = sys.argv[3]
max_dist = float(sys.argv[4])
chains = sys.argv[5]

########
pdb_columns = ['name', 'resn', 'resid', 'Chain', 'X', 'Y', 'Z', 'Element']
# Read in PDB file one line at a time, if the first four letter ar ATOM or HETA then it will parse the data into the 
# data frame, using the atome index int he PDB as the row index in the data frame. 
AAA2A = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V' }
A2AAA = { "I": "ILE", "A": "ALA", "V": "VAL", "L": "LEU", "M": "MET", "T" : "THR","F" : "PHE", "Y": "TYR" }
Methylpseudo = {"LEU" : ['QD1','QD2'], "ILE" : ['QD1'], "VAL" : ['QG1','QG2'], "THR" : ['QG2'], "MET" : ['QE'], "ALA" : ['QB']}
Aromaticpseudo = {"PHE" : ['HE1','HE2'], "TYR": ['HE1','HE2']}
Methylcarbons = {"LEU" : ['CD1','CD2'], "ILE" : ['CD1'], "VAL" : ['CG1','CG2'], "THR" : ['CG2'], "MET" : ['CE'], "ALA" : ['CB']}
Methylprotons = {"LEU" : ['HD11', 'HD12', 'HD13','HD21', 'HD22', 'HD23'], "ILE" : ['HD11', 'HD12', 'HD13'], "VAL" : ['HG11', 'HG12', 'HG13','HG21', 'HG22', 'HG23'], "THR" : ['HG21', 'HG22', 'HG23'], "MET" : ['HE1', 'HE2', 'HE3'], "ALA" : ['HB1', 'HB2', 'HB3']}
Aromaticcarbons = {"PHE" : ['CE1', 'CE2'], "TYR": ['CE1','CE2']}
Aromaticprotons = {"PHE" : ['HE1', 'HE2'], "TYR": ['HE1','HE2']}
Pseudo2Prot = {'ALAQB':['HB1', 'HB2', 'HB3'], 'ARGQB':['HB2', 'HB3'], 'ARGQG':['HG2', 'HG3'], 'ARGQD':['HD2', 'HD3'], 'ARGQH1':['HH11', 'HH12'], 'ARGQH2':['HH21', 'HH22'], 'ASNQB':['HB2', 'HB3'], 'ASNQD2':['HD21', 'HD22'], 'ASPQB':['HB2', 'HB3'], 'CYSQB':['HB2', 'HB3'], 'CYSSQB':['HB2', 'HB3'], 'GLNQB':['HB2', 'HB3'], 'GLNQG':['HG2', 'HG3'], 'GLNQE2':['HE21', 'HE22'], 'GLUQB':['HB2', 'HB3'], 'GLUQG':['HG2', 'HG3'], 'GLYQA':['HA2', 'HA3'], 'HISQB':['HB2', 'HB3'], 'ILEQG1':['HG11', 'HG12', 'HG13'], 'ILEQG2':['HG21', 'HG22', 'HG23'], 'ILEQD1':['HD11', 'HD12', 'HD13'], 'LEUQB':['HB2', 'HB3'], 'LEUQD1':['HD11', 'HD12', 'HD13'], 'LEUQD2':['HD21', 'HD22', 'HD23'], 'LEUQQD':['HD11', 'HD12', 'HD13','HD21', 'HD22', 'HD23'], 'LYSQG':['HG2', 'HG3'], 'LYSQD':['HD2', 'HD3'], 'LYSQE':['HE2', 'HE3'], 'LYSQZ':['HZ1','HZ2', 'HZ3'], 'METQB':['HB2', 'HB3'], 'METQG':['HG2', 'HG3'], 'METQE':['HE1', 'HE2', 'HE3'], 'PHEQB':['HB2', 'HB3'], 'PHEQD':['HD1', 'HD2'], 'PHEQE':['HE1', 'HE2'],'PROQB':['HB2', 'HB3'], 'PROQG':['HG2', 'HG3'], 'PROQD':['HD2', 'HD3'], 'SERQB':['HB2', 'HB3'], 'THRQG2':['HG21', 'HG22', 'HG23'], 'TRPQB':['HB2', 'HB3'], 'TYRQB':['HB2', 'HB3'], 'TYRQD':['HD1', 'HD2'], 'TYRQE':['HE1', 'HE2'], 'VALQB':['HB2', 'HB3'], 'VALQG1':['HG11', 'HG12', 'HG13'], 'VALQG2':['HG21', 'HG22', 'HG23'], 'VALQQG':['HG11', 'HG12', 'HG13','HG21', 'HG22', 'HG23']}
Aromaticpseudo = {"PHE" : ['QE'], "TYR": ['QE']}
# AllowedProtons = {
# 'ALA':['H','HA','QB'],
# 'ARG':['H','HA','HB2','HB3','HB2', 'HB3','HG2', 'HG3','HD2', 'HD3','HH11', 'HH12','HH21', 'HH22'],
# 'ASN':['H','HA','HB2', 'HB3','HD21', 'HD22'],
# 'ASP':['H','HA','HB2', 'HB3'],
# 'CYS':['H','HA','HB2', 'HB3'],
# 'GLN':['H','HA','HB2', 'HB3','HG2', 'HG3','HE21', 'HE22'],
# 'GLU':['H','HA','HB2', 'HB3','HG2', 'HG3'],
# 'GLY':['H','HA2','HA3'],
# 'GQA':['HA2', 'HA3'],
# 'HIS':['H','HA','HB2', 'HB3','HD1','HE1','HD2','HE2'],
# 'ILE':['H','HA','HB','QG2','HG12','HG13','QD1'],
# 'LEU':['H','HA','HB2', 'HB3','QD1','QD2','HG'],
# 'LYS':['H','HA','HG2', 'HG3','HD2', 'HD3','HE2', 'HE3','QZ'],
# 'MET':['H','HA','HB2', 'HB3','HG2', 'HG3','QE'],
# 'PHE':['H','HA','HB2', 'HB3','HD1', 'HD2','HE1', 'HE2'],
# 'PRO':['HA','HB2', 'HB3','HG2', 'HG3','HD2', 'HD3'],
# 'SER':['H','HA','HB2', 'HB3'],
# 'THR':['H','HA','HB','QG2'],
# 'VAL':['H','HA','HB','QG1','QG2'],
# 'TRP':['H','HA','HB2', 'HB3','HD1','HE1','HZ2','HH2','HZ3','HE3'],
# 'TYR':['H','HA','HB2', 'HB3','HD1', 'HD2','HE1', 'HE2']}
AllowedProtons = {
'ALA':['H','HA','QB'],
'ARG':['H','HA','HB2','HB3','HB2', 'HB3','HG2', 'HG3','HD2', 'HD3'],
'ASN':['H','HA','HB2', 'HB3'],
'ASP':['H','HA','HB2', 'HB3'],
'CYS':['H','HA','HB2', 'HB3'],
'GLN':['H','HA','HB2', 'HB3','HG2', 'HG3'],
'GLU':['H','HA','HB2', 'HB3','HG2', 'HG3'],
'GLY':['H','HA2','HA3'],
'GQA':['HA2', 'HA3'],
'HIS':['H','HA','HB2', 'HB3','HD1','HE1','HD2','HE2'],
'ILE':['H','HA','HB','QG2','HG12','HG13','QD1'],
'LEU':['H','HA','HB2', 'HB3','QD1','QD2','HG'],
'LYS':['H','HA','HG2', 'HG3','HD2', 'HD3','HE2', 'HE3'],
'MET':['H','HA','HB2', 'HB3','HG2', 'HG3','QE'],
'PHE':['H','HA','HB2', 'HB3','HD1', 'HD2','HE1', 'HE2'],
'PRO':['HA','HB2', 'HB3','HG2', 'HG3','HD2', 'HD3'],
'SER':['H','HA','HB2', 'HB3'],
'THR':['H','HA','HB','QG2'],
'VAL':['H','HA','HB','QG1','QG2'],
'TRP':['H','HA','HB2', 'HB3','HD1','HE1','HZ2','HH2','HZ3','HE3'],
'TYR':['H','HA','HB2', 'HB3','HD1', 'HD2','HE1', 'HE2']}
def getDistance(donor, acceptor,PDBdict):
	resn1, resi1, name1 = donor.split()
	resn2, resi2, name2 = acceptor.split()
	d = 0
	dist = []
	if resn1+name1 in Pseudo2Prot.keys():
		atoms1list = ['{:4} {:4} {:4}'.format(resn1.replace('HIST','HIS'),resi1,atom) for atom in Pseudo2Prot[resn1.replace('HIST','HIS')+name1]]
	if resn1+name1  not in Pseudo2Prot.keys():
		atoms1list = [donor]
	if resn2+name2 in Pseudo2Prot.keys():
		atoms2list = ['{:4} {:4} {:4}'.format(resn2.replace('HIST','HIS'),resi2,atom) for atom in Pseudo2Prot[resn2.replace('HIST','HIS')+name2]]
	if resn2+name2 not in Pseudo2Prot.keys():
		atoms2list = [acceptor]
	for a1 in atoms1list:
		for a2 in atoms2list:
			(x1,y1,z1) = PDBdict[a1]
			(x2,y2,z2) = PDBdict[a2]
			# d = d + np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**-6
			dist.append(np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2))
	reff = np.min(dist)
	# reff = np.round(d**(-1/6),2)
	# print('used %s %s d=%1.2f' %(donor, acceptor,reff))
	return reff

if 'HB2' in open(in_pdb).read():
	Methyl_groups = Methylpseudo
	Aromatic_groups = Aromaticpseudo
	bbatom = 'H'
if 'HB2' not in open(in_pdb).read():
	Methyl_groups = Methylcarbons
	Aromatic_groups = Aromaticcarbons
	bbatom = 'N'

Backbone, Aromatics, Methyls, Sequence, AllAtoms = [], [], [], [],[]
Coor = {}

if len(chains) == 1:
	with open(in_pdb) as In_pdb:
		for line in In_pdb:
			if line[0:4] == "ATOM" or line[0:4] == 'HETA':
				if line[21] == chains:
					index = '{:4} {:4} {:4}'.format(line[17:20].strip(),line[22:26].strip(),line[12:16].strip())
					Coor[index] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
					if line[17:20].strip() + line[22:26].strip() not in Sequence and line[17:20].strip() in AAA2A.keys():
						Sequence.append(line[17:20].strip() + line[22:26].strip())
						if line[17:20].strip() != 'PRO':
							Backbone.append('{:4} {:4} {:4}'.format(line[17:20].strip(),line[22:26].strip(),bbatom))
						if line[17:20].strip() in Methyl_groups.keys():
							for atom in Methyl_groups[line[17:20].strip()]:
								Methyls.append('{:4} {:4} {:4}'.format(line[17:20].strip(),line[22:26].strip(),atom))
						if line[17:20].strip() in Aromatic_groups.keys():
							for atom in Aromatic_groups[line[17:20].strip()]:
								Aromatics.append('{:4} {:4} {:4}'.format(line[17:20].strip(),line[22:26].strip(),atom))
						for atom in AllowedProtons[line[17:20].strip()]:
							AllAtoms.append('{:4} {:4} {:4}'.format(line[17:20].strip(),line[22:26].strip(),atom))
	droplist = []
	for atom in Backbone:
		if atom not in Coor.keys():
			droplist.append(atom)
	for atom in AllAtoms:
		if atom.split()[-1][0] != 'Q':
			if atom not in Coor.keys():
				droplist.append(atom)
	for atom in droplist:
		try:
			AllAtoms.remove(atom)
		except:
			pass
		try:
			Backbone.remove(atom)
		except:
			pass
	Doubledist = []
	D3,D5,DLshort, DLmed, DLlon,DLintra = 0,0,0,0,0,0
	for atom1 in AllAtoms:
		for atom2 in AllAtoms:
			if atom1 != atom2:
				dist = getDistance(atom1, atom2,Coor)
				if dist < max_dist:
					Doubledist.append(dist)
					if dist <= 3.0:
						D3+=0.5
					elif dist > 3.0 and dist <= 5.0:
						D5+=0.5
					diff = abs(float(atom1.split()[1]) - float(atom2.split()[1]))
					if diff == 0:DLintra+=0.5
					elif diff ==1:DLshort+=0.5
					elif diff >1 and diff < 5:DLmed+=0.5
					elif diff >5:DLlon+=0.5
	print('Double labeled sample in H2O')
	print('   Number of distances {:.0f}'.format(len(Doubledist)/2))
	print('   Median distance is {:2.1f}'.format(median(Doubledist)))
	print('   Average of {:3.1f} distance restraints per residue'.format((len(Doubledist)/len(Sequence))/2.0))
	print('   {:.0f} distance < 3A'.format(D3))
	print('   {:.0f} 3A > distance < 5A'.format(D5))
	print('   with intramolecular assignment |i-j|=0 : {:5.0f} ({:2.0f}%)'.format(DLintra,100*(DLintra/(len(Doubledist)/2))))
	print('   with short-range assignment    |i-j|=1 : {:5.0f} ({:2.0f}%)'.format(DLshort,100*(DLshort/(len(Doubledist)/2))))
	print('   with medium-range assignment 1<|i-j|<5 : {:5.0f} ({:2.0f}%)'.format(DLmed,100*(DLmed/(len(Doubledist)/2))))
	print('   with long-range assignment     |i-j|>=5: {:5.0f} ({:2.0f}%)'.format(DLlon,100*(DLlon/(len(Doubledist)/2))))
	print()
	MethyLabeled = []
	MethyLabeled.extend(Methyls)
	MethyLabeled.extend(Backbone)
	if len(Aromatics) != 0: MethyLabeled.extend(Aromatics)

	MElabeled = []
	D3,D5,MEshort, MEmed, MElon,MEintra = 0,0,0,0,0,0
	for atom1 in MethyLabeled:
		for atom2 in MethyLabeled:
			if atom1 != atom2:
				dist = getDistance(atom1, atom2,Coor)
				if dist < max_dist:
					if labeling_type == 'm' and atom1.split()[1] == atom2.split()[1] and atom1.split()[2][0] == atom2.split()[2][0]:
						pass
					else:
						MElabeled.append(dist)
						if dist <= 3.0:
							D3+=0.5
						elif dist > 3.0 and dist <= 5.0:
							D5+=0.5
						diff = abs(float(atom1.split()[1]) - float(atom2.split()[1]))
						if diff == 0:MEintra+=0.5
						elif diff ==1:MEshort+=0.5
						elif diff >1 and diff < 5:MEmed+=0.5
						elif diff >= 5:MElon+=0.5
	print('Deuterated ILVMAT-YF labeled sample in H2O')
	print('   Number of distances {:.0f}'.format(len(MElabeled)/2))
	print('   Median distance is {:2.1f}'.format(median(MElabeled)))
	print('   Average of {:3.1f} distance restraints per residue'.format((len(MElabeled)/len(Sequence))/2.0))
	print('   {:.0f} distance < 3A'.format(D3))
	print('   {:.0f} 3A > distance < 5A'.format(D5))
	print('   with intramolecular assignment |i-j|=0 : {:5.0f} ({:2.0f}%)'.format(MEintra,100*(MEintra/(len(MElabeled)/2))))
	print('   with short-range assignment    |i-j|=1 : {:5.0f} ({:2.0f}%)'.format(MEshort,100*(MEshort/(len(MElabeled)/2))))
	print('   with medium-range assignment 1<|i-j|<5 : {:5.0f} ({:2.0f}%)'.format(MEmed,100*(MEmed/(len(MElabeled)/2))))
	print('   with long-range assignment     |i-j|>=5: {:5.0f} ({:2.0f}%)'.format(MElon,100*(MElon/(len(MElabeled)/2))))
	print()

	print('Total of {:} labeled Methyls'.format(len(Methyls)))
	Me_Me_Dist = []
	for atom1 in Methyls:
		resn1, resi1, name1 = atom1.split()
		for atom2 in Methyls:
			resn2, resi2, name2 = atom2.split()
			if atom1 != atom2: 
				dist = getDistance(atom1, atom2,Coor)
				if dist < max_dist:
					if labeling_type == 'm' and resi1 == resi2 and name1[0] == name2[0]:
						pass
					else:
						Me_Me_Dist.append(dist)

	Hall_Me_Dist = []
	for atom1 in Methyls:
		resn1, resi1, name1 = atom1.split()
		for atom2 in MethyLabeled:
			resn2, resi2, name2 = atom2.split()
			if atom1 != atom2: 
				dist = getDistance(atom1, atom2,Coor)
				if dist < max_dist: 
					if labeling_type == 'm' and resi1 == resi2 and name1[0] == name2[0]:
						pass
					else:
						Hall_Me_Dist.append(dist)

	Ca_Me_Dist = []
	if len(Aromatics) != 0:
		print('Total of {:} YF Aromatics labeled'.format(len(Aromatics)))
		for atom2 in Methyls:
			resn1, resi1, name1 = atom1.split()
			for atom1 in Aromatics:
				resn2, resi2, name2 = atom2.split()
				if atom1 != atom2: 
					dist = getDistance(atom1, atom2,Coor)
					if dist < max_dist: 
						Ca_Me_Dist.append(dist)
		Hall_Ca_Dist = []
		for atom1 in Aromatics:
			resn1, resi1, name1 = atom1.split()
			for atom2 in MethyLabeled:
				resn2, resi2, name2 = atom2.split()
				if atom1 != atom2: 
					dist = getDistance(atom1, atom2,Coor)
					if dist < max_dist: 
						Hall_Ca_Dist.append(dist)

	print('Total of {:} amide groups'.format(len(Backbone)))
	Hn_NHn_Dist = []
	for atom1 in Backbone:
		resn1, resi1, name1 = atom1.split()
		for atom2 in Backbone:
			resn2, resi2, name2 = atom2.split()
			if atom1 != atom2: 
				dist = getDistance(atom1, atom2,Coor)
				if dist < max_dist: 
					Hn_NHn_Dist.append(dist)

	Hall_NHn_Dist = []
	for atom1 in Backbone:
		resn1, resi1, name1 = atom1.split()
		for atom2 in MethyLabeled:
			resn2, resi2, name2 = atom2.split()
			if atom1 != atom2: 
				dist = getDistance(atom1, atom2,Coor)
				if dist < max_dist: 
					Hall_NHn_Dist.append(dist)

	Cm_NHn_Dist = []
	for atom1 in Backbone:
		resn1, resi1, name1 = atom1.split()
		for atom2 in Methyls:
			resn2, resi2, name2 = atom2.split()
			if atom1 != atom2:
				dist = getDistance(atom1, atom2,Coor)
				if dist < max_dist: 
					Cm_NHn_Dist.append(dist)

	print('{:>30} {:^9} {:^9} {:^9} {:^9} {:^9} {:^9} {:^9}'.format('NOESY','Cm-CmHm','Hall-CmHm','Ca-CmHm','Hall-CaHa','Hn-NHn','Hall-NHn','Cm-NHn'))
	print('{:30} {:^9d} {:^9d} {:^9d} {:^9d} {:^9d} {:^9d} {:^9d}'.format('Number of crosspeaks',len(Me_Me_Dist),len(Hall_Me_Dist),len(Ca_Me_Dist),len(Hall_Ca_Dist),len(Hn_NHn_Dist),len(Hall_NHn_Dist),len(Cm_NHn_Dist)))
	print('{:30} {:^9.1f} {:^9.1f} {:^9.1f} {:^9.1f} {:^9.1f} {:^9.1f} {:^9.1f}'.format('Median Distance (Å)',median(Me_Me_Dist),median(Hall_Me_Dist),median(Ca_Me_Dist),median(Hall_Ca_Dist),median(Hn_NHn_Dist),median(Hall_NHn_Dist),median(Cm_NHn_Dist)))
	print('{:30} {:^9.1f} {:^9.1f} {:^9.1f} {:^9.1f} {:^9.1f} {:^9.1f} {:^9.1f}'.format('Average # crosspeaks per group',len(Me_Me_Dist)/len(Methyls),len(Hall_Me_Dist)/len(Methyls),len(Ca_Me_Dist)/len(Methyls),len(Hall_Ca_Dist)/len(Aromatics),len(Hn_NHn_Dist)/len(Backbone),len(Hall_NHn_Dist)/len(Backbone),len(Cm_NHn_Dist)/len(Backbone)))



