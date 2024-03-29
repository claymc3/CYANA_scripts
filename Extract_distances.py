'''
# ------------------------------------------------------------------------------
#
# Created by : Mary Clay PhD
# e-mail: mary.clay@stjude.org
# St Jude Children's Research Hospital 
# Department of Structural Biology Memphis, TN 
#
# ------------------------------------------------------------------------------

'''
import os
import sys
import pandas as pd
import numpy as np
'''
CNS noe.tbl fomat
 assign (resid 632 and name HN   )  (resid 633 and name HN   )  3.24  1.44  1.44 !peak 1 #plist 1 #SUP 0.91 #QF 0.91 
 line.split()
 0 assign
 1 (reside 
 2 residue 1 index 
 3 and
 4 name
 5 residue 1 atom name 
 6 )
 7 (resid
 8 residue 2 index
 9 and 
 10 name
 11 residue 2 atom name 
 12 )
 13 d (the distance)
 14 dminus 
 15 dplus 
 16 !peak
 17 peak number 
 18 #plist
 19 plist number 
 20 #SUP
 21 SUP value
 22 #QF (optional)
 23 QF values (optional)

d (the distance), and dminus, and dplus

CYANA upl/lol format 
632 LYS  H     633 ILE  H       4.46  #peak 1 #plist 1 #SUP 0.91 #QF 0.91
'''
DistancesDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv','upl','plist','peak'])


Pseudo2Prot = {'ALAQB':['HB1', 'HB2', 'HB3'], 'ARGQB':['HB2', 'HB3'], 'ARGQG':['HG2', 'HG3'], 'ARGQD':['HD2', 'HD3'], 'ARGQH1':['HH11', 'HH12'], 'ARGQH2':['HH21', 'HH22'], 'ASNQB':['HB2', 'HB3'], 'ASNQD2':['HD21', 'HD22'], 'ASPQB':['HB2', 'HB3'], 'CYSQB':['HB2', 'HB3'], 'CYSSQB':['HB2', 'HB3'], 'GLNQB':['HB2', 'HB3'], 'GLNQG':['HG2', 'HG3'], 'GLNQE2':['HE21', 'HE22'], 'GLUQB':['HB2', 'HB3'], 'GLUQG':['HG2', 'HG3'], 'GLYQA':['HA2', 'HA3'], 'HISQB':['HB2', 'HB3'], 'ILEQG1':['HG12', 'HG13'], 'ILEQG2':['HG21', 'HG22', 'HG23'], 'ILEQD1':['HD11', 'HD12', 'HD13'], 'LEUQB':['HB2', 'HB3'], 'LEUQD1':['HD11', 'HD12', 'HD13'], 'LEUQD2':['HD21', 'HD22', 'HD23'], 'LEUQQD':['HD11', 'HD12', 'HD13','HD21', 'HD22', 'HD23'], 'LYSQG':['HG2', 'HG3'], 'LYSQD':['HD2', 'HD3'], 'LYSQE':['HE2', 'HE3'], 'LYSQZ':['HZ2', 'HZ3'], 'METQB':['HB2', 'HB3'], 'METQG':['HG2', 'HG3'], 'METQE':['HE1', 'HE2', 'HE3'], 'PHEQB':['HB2', 'HB3'], 'PHEQD':['HD1', 'HD2'], 'PHEQE':['HE1', 'HE2'],'PROQB':['HB2', 'HB3'], 'PROQG':['HG2', 'HG3'], 'PROQD':['HD2', 'HD3'], 'SERQB':['HB2', 'HB3'], 'THRQG2':['HG21', 'HG22', 'HG23'], 'TRPQB':['HB2', 'HB3'], 'TYRQB':['HB2', 'HB3'], 'TYRQD':['HD1', 'HD2'], 'TYRQE':['HE1', 'HE2'], 'VALQB':['HB2', 'HB3'], 'VALQG1':['HG11', 'HG12', 'HG13'], 'VALQG2':['HG21', 'HG22', 'HG23'], 'VALQQG':['HG11', 'HG12', 'HG13','HG21', 'HG22', 'HG23']}
PseudoIDs = {'ALA':['ALAQB'],'ILE':['ILEQD1'],'LEU':['LEUQD1','LEUQD2'], 'VAL':['VALQG1','VALQG2'],'MET':['METQE'], 'THR':['THRQG2']}

Prot2Pseudo = {'ALAHB': 'QB','ALAHB1': 'QB', 'ALAHB2': 'QB', 'ALAHB3': 'QB', 'ARGHB2': 'QB', 'ARGHB3': 'QB', 'ARGHG2': 'QG', 'ARGHG3': 'QG', 'ARGHD2': 'QD', 'ARGHD3': 'QD', 'ARGHH11': 'QH1', 'ARGHH12': 'QH1', 'ARGHH21': 'QH2', 'ARGHH22': 'QH2', 'ASNHB2': 'QB', 'ASNHB3': 'QB', 'ASNHD21': 'QD2', 'ASNHD22': 'QD2', 'ASPHB2': 'QB', 'ASPHB3': 'QB', 'CYSHB2': 'QB', 'CYSHB3': 'QB', 'CYSSHB2': 'QB', 'CYSSHB3': 'QB', 'GLNHB2': 'QB', 'GLNHB3': 'QB', 'GLNHG2': 'QG', 'GLNHG3': 'QG', 'GLNHE21': 'QE2', 'GLNHE22': 'QE2', 'GLUHB2': 'QB', 'GLUHB3': 'QB', 'GLUHG2': 'QG', 'GLUHG3': 'QG', 'GLYHA2': 'QA', 'GLYHA3': 'QA', 'HISHB2': 'QB', 'HISHB3': 'QB', 'ILEHG1': 'QG1', 'ILEHG11': 'QG1', 'ILEHG12': 'QG1', 'ILEHG13': 'QG1', 'ILEHG21': 'QG2', 'ILEHG22': 'QG2', 'ILEHG23': 'QG2', 'ILEHD11': 'QD1', 'ILEHD1': 'QD1', 'ILEHD12': 'QD1', 'ILEHD13': 'QD1', 'LEUHB2': 'QB', 'LEUHB3': 'QB', 'LEUHD1': 'QD1', 'LEUHD11': 'QD1', 'LEUHD12': 'QD1', 'LEUHD13': 'QD1', 'LEUHD2': 'QD2', 'LEUHD21': 'QD2', 'LEUHD22': 'QD2', 'LEUHD23': 'QD2', 'LYSHG2': 'QG', 'LYSHG3': 'QG', 'LYSHD2': 'QD', 'LYSHD3': 'QD', 'LYSHE2': 'QE', 'LYSHE3': 'QE', 'LYSHZ2': 'QZ', 'LYSHZ3': 'QZ', 'METHB2': 'QB', 'METHB3': 'QB', 'METHG2': 'QG', 'METHG3': 'QG', 'METHE': 'QE', 'METHE1': 'QE', 'METHE2': 'QE', 'METHE3': 'QE', 'PHEHB2': 'QB', 'PHEHB3': 'QB', 'PHEHD1': 'QD', 'PHEHD2': 'QD', 'PHEHE1': 'QE', 'PHEHE2': 'QE', 'PROHB2': 'QB', 'PROHB3': 'QB', 'PROHG2': 'QG', 'PROHG3': 'QG', 'PROHD2': 'QD', 'PROHD3': 'QD', 'SERHB2': 'QB', 'SERHB3': 'QB', 'THRHG2': 'QG2', 'THRHG21': 'QG2', 'THRHG22': 'QG2', 'THRHG23': 'QG2', 'TRPHB2': 'QB', 'TRPHB3': 'QB', 'TYRHB2': 'QB', 'TYRHB3': 'QB', 'TYRHD1': 'QD', 'TYRHD2': 'QD', 'TYRHE1': 'QE', 'TYRHE2': 'QE', 'VALHB2': 'QB', 'VALHB3': 'QB', 'VALHG11': 'QG1', 'VALHG12': 'QG1', 'VALHG13': 'QG1', 'VALHG21': 'QG2', 'VALHG22': 'QG2', 'VALHG23': 'QG2'}
replacements ={
'ALAHA':'CA','ALAQB':'CB','ALAHB1':'CB','ALAHB2':'CB','ALAHB3':'CB',
'CYSHA':'CA','CYSHB2':'CB','CYSHB3':'CB','CYSQB':'CB',
'ASPHA':'CA','ASPHB2':'CB','ASPHB3':'CB','ASPQB':'CB',
'GLUHA':'CA','GLUHB2':'CB','GLUHB3':'CB','GLUQB':'CB','GLUHG2':'CG','GLUHG3':'CG','GLUQG':'CG',
'PHEHA':'CA','PHEHB2':'CB','PHEHB3':'CB','PHEQB':'CB','PHEQD':'CD1,CD2','PHEQE':'CE1,CE2','PHEHD1':'CD1','PHEHE1':'CE1','PHEHZ':'CZ','PHEHE2':'CE2','PHEHD2':'CD2',
'GLYHA2':'CA','GLYHA3':'CA','GLYQA':'CA',
'HISHA':'CA','HISHB2':'CB','HISHB3':'CB','HISQB':'CB','HISHD1':'ND1','HISHE2':'NE2','HISHD2':'CD2','HISHE1':'CE1',
'ILEHA':'CA','ILEHB':'CB','ILEQG2':'CG2','ILEHG21':'CG2','ILEHG22':'CG2','ILEHG23':'CG2','ILEHG12':'CG1','ILEHG13':'CG1','ILEQG1':'CG1','ILEQD1':'CD1','ILEHD11':'CD1','ILEHD12':'CD1','ILEHD13':'CD1',
'LYSHA':'CA','LYSHB2':'CB','LYSHB3':'CB','LYSQB':'CB','LYSHG2':'CG','LYSHG3':'CG','LYSHD2':'CD ','LYSHD3':'CD ','LYSQD':'CD ','LYSHE2':'CE','LYSHE3':'CE','LYSQE':'CE','LYSHZ1':'NZ','LYSHZ2':'NZ','LYSHZ3':'NZ','LYSQZ':'NZ',
'LEUHA':'CA','LEUHB2':'CB','LEUHB3':'CB','LEUQB':'CB','LEUHG':'CG','LEUHD11':'CD1','LEUHD12':'CD1','LEUHD13':'CD1','LEUQD1':'CD1','LEUHD21':'CD2','LEUHD22':'CD2','LEUHD23':'CD2','LEUQD2':'CD2','LEUQQD':'CD2,CD1',
'METHA':'CA','METHB2':'CB','METHB3':'CB','METQB':'CB','METHG2':'CG','METHG3':'CG','METQG':'CG','METQE':'CE','METHE1':'CE','METHE2':'CE','METHE3':'CE',
'ASNHA':'CA','ASNHB2':'CB','ASNHB3':'CB','ASNQB':'CB','ASNHD21':'ND2','ASNHD22':'ND2','ASNQD2':'ND2',
'PROHA':'CA','PROHB2':'CB','PROHB3':'CB','PROQB':'CB','PROHG2':'CG','PROHG3':'CG','PROQG':'CG','PROHD2':'CD','PROHD3':'CD','PROQD':'CD',
'GLNHA':'CA','GLNHB2':'CB','GLNHB3':'CB','GLNQB':'CB','GLNHG2':'CG','GLNHG3':'CG','GLNQG':'CG','GLNHE21':'NE2','GLNHE22':'NE2','GLNQE2':'NE2',
'ARGHA':'CA','ARGHB2':'CB','ARGHB3':'CB','ARGQB':'CB','ARGHG2':'CG','ARGHG3':'CG','ARGQG':'CG','ARGHD2':'CD','ARGHD3':'CD','ARGQD':'CD','ARGHE':'NE','ARGHH11':'NH1','ARGHH12':'NH1','ARGQH1':'NH1','ARGHH21':'NH2','ARGHH22':'NH2','ARGQH2':'NH2',
'SERHA':'CA','SERHB2':'CB','SERHB3':'CB','SERHG':'OG',
'THRHA':'CA','THRHB':'CB','THRHG1':'OG1','THRHG21':'CG2','THRHG22':'CG2','THRHG23':'CG2','THRQG2':'CG2',
'VALHA':'CA','VALHB':'CB','VALHG11':'CG1','VALHG12':'CG1','VALHG13':'CG1','VALQG1':'CG1','VALHG21':'CG2','VALHG22':'CG2','VALHG23':'CG2','VALQG2':'CG2','VALQQG':'CG1,CG2',
'TRPHA':'CA','TRPHB2':'CB','TRPHB3':'CB','TRPQB':'CB','TRPHD1':'CD1','TRPHE3':'CE3','TRPHE1':'NE1','TRPHZ3':'CZ3','TRPHZ2':'CZ2','TRPHH2':'CH2',
'TYRHA':'CA','TYRHB2':'CB','TYRHB3':'CB','TYRQB':'CB','TYRQD':'CD1,CD2','TYRQE':'CE1,CE2','TYRHD1':'CD1','TYRHE1':'CE1','TYRHE2':'CE2','TYRHD2':'CD2','TYRHH':'OH'}

AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "CYSS":"C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H","HIST": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P","cPRO":"P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S'}
A_dict = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS','I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

cwd = '/Users/mclay1/FGFR2/FGFR2_467-768_C491A/Structure_Calc/cyana_70/'
# cwd = '/Users/mclay1/FGFR3_structure/cyana_23/'

in_pdb = cwd + 'CNS/refinedPDB/FGFR2_r70_cya.pdb'
pdb_name = in_pdb.split('/')[-1].replace('.pdb','')
fupl = cwd + 'final.upl'
calc = cwd + 'CALC.cya'
init = cwd + 'init.cya'
manualongcons = [line.strip() for line in open(calc).readlines() if line.strip() and '.upl' in line][0].split()[2].split(',')
upls = [con for con in manualongcons if 'upl' in con and 'hbond' not in con]
fovw = cwd + 'final.ovw'
Starts,Ends = [], []
for mnum in range(1,21,1):
	start = open(in_pdb).readlines().index('MODEL' + '%9s\n' %str(mnum))
	exec('Coor' + str(mnum) + ' = {}')
	Starts.append(start)
	Ends.append(start-1)
Ends.append(len(open(in_pdb).readlines()))
Ends = Ends[1:]
pdb = open(in_pdb).readlines()
n = 0
seq = [line.strip().split() for line in open(cwd + open(init).readlines()[0].strip().split(':=')[-1] + '.seq').readlines() if '#' != line[0]]
Seqdict = {}
Sequence, ASequence = [],[]

for (start,end) in zip(Starts,Ends):
	n+=1
	print('Reading coordinates for model {:d}'.format(n))
	Coor = eval('Coor' + str(n))
	for x in range(start,end,1):
		line = pdb[x]
		if line[0:4] == "ATOM" or line[0:4] == 'HETA':
			index = '{:}{:}-{:}'.format(AAA_dict[line[17:20].strip()],line[22:26].strip(),line[12:16].strip())
			Coor[index] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
for mnum in range(1,21,1):
	Coor = eval('Coor' + str(mnum))
	for resn,resi in seq:
		if resn in PseudoIDs.keys():
			for pseudo in PseudoIDs[resn]:
				index = '{:}{:}-{:}'.format(AAA_dict[resn],resi,pseudo[3:])
				(x1,y1,z1) = Coor['{:}{:}-{:}'.format(AAA_dict[resn],resi,Pseudo2Prot[pseudo][0])]
				(x2,y2,z2) = Coor['{:}{:}-{:}'.format(AAA_dict[resn],resi,Pseudo2Prot[pseudo][1])]
				(x3,y3,z3) = Coor['{:}{:}-{:}'.format(AAA_dict[resn],resi,Pseudo2Prot[pseudo][2])]
				Coor[index] = [(x1+x2+x3)/3.0, (y1+y2+y3)/3.0, (z1+z2+z3)/3.0]

def getDistance(donor, acceptor,PDBdict):
	resn1, name1 = donor.split('-')
	resn2, name2 = acceptor.split('-')
	d = 0
	if 'Q' in name1: atoms1list = ['{:}-{:}'.format(resn1,atom) for atom in Pseudo2Prot[A_dict[resn1[0]]+name1]]
	if 'Q' in name2: atoms2list = ['{:}-{:}'.format(resn2,atom) for atom in Pseudo2Prot[A_dict[resn2[0]]+name2]]
	if 'Q' not in name1: atoms1list = [donor]
	if 'Q' not in name2: atoms2list = [acceptor]
	for a1 in atoms1list:
		for a2 in atoms2list:
			(x1,y1,z1) = PDBdict[a1]
			(x2,y2,z2) = PDBdict[a2]
			d = d + np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**-6
	reff = np.round(d**(-1/6),2)
	# print('used %s %s d=%1.2f' %(donor, acceptor,reff))
	return reff

violdict = {}
for line in open(fovw).readlines():
	if line[4:9] == 'Upper':
		dviol = line.split()
		cons1 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[dviol[2]],dviol[3],dviol[1],AAA_dict[dviol[6]],dviol[7],dviol[5])
		cons2 = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[dviol[6]],dviol[7],dviol[5],AAA_dict[dviol[2]],dviol[3],dviol[1])
		outline = ' #viol in {:} by {:}\n'.format(dviol[9],dviol[10])
		violdict[cons1] = outline
		violdict[cons2] = outline




DistancesDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv','upl','violation','num viol'])
DiffDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
print('Extracting Distances from %s' %fupl)
## Extract distances for the final upl no modification, and excluding ambiguous restraints
for line in open(fupl).readlines():
	if 'SUP' not in line:
		pass 
	else:
		cns = line.strip().split()
		atom1 = '{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2])
		atom2 = '{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],cns[5])
		connection = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
		DistancesDF.loc[connection,'upl'] = cns[6]
		diff = []
		for mnum in range(1,21,1):
			Coor = eval('Coor' + str(mnum))
			d = float(getDistance(atom1, atom2, Coor))
			DistancesDF.loc[connection,mnum] = d
			if np.round(d - float(cns[6]),2) > 0.1:
				DiffDF.loc[connection,mnum] = np.round(d - float(cns[6]),2)
				diff.append(np.round(d - float(cns[6]),2))
		if connection in violdict.keys():
			DistancesDF.loc[connection,'violation'] = violdict[connection]
		DistancesDF.loc[connection,'num viol'] = len(diff)
DistancesDF['mean'] = DistancesDF[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]].mean(axis=1)
DistancesDF['stdv'] = DistancesDF[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]].std(axis=1)
decimals = pd.Series([2, 2,2,2], index=['mean','stdv','mean diff','max diff'])
DistancesDF['mean diff'] = DiffDF.mean(axis=1)
DistancesDF['max diff'] = DiffDF.max(axis=1)
DistancesDF.round(2)
DistancesDF.to_csv(pdb_name + '_distances_proton.csv')
print(DistancesDF)
print(DistancesDF.shape)

## Now look at the other input upls 
# for upl in upls:
# 	print(upl)
# 	DiffDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
# 	DistancesDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv','upl'])
# 	for line in open(cwd+upl).readlines():
# 		if '#'in line[0:4]:
# 			pass 
# 		else:
# 			cns = line.strip().split()
# 			atom1 = '{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2])
# 			atom2 = '{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],cns[5])
# 			connection = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
# 			DistancesDF.loc[connection,'upl'] = float(cns[6])
# 			for mnum in range(1,21,1):
# 				Coor = eval('Coor' + str(mnum))
# 				d = getDistance(atom1, atom2,Coor)
# 				DistancesDF.loc[atom1 + ' ' + atom2,mnum] = d
# 				if np.round(d - float(cns[6]),2) > 0.1:
# 					DiffDF.loc[connection,mnum] = np.round(d - float(cns[6]),2)
# 			if connection in violdict.keys():
# 				DistancesDF.loc[connection,'violation'] = violdict[connection]
# 	DistancesDF['mean'] = DistancesDF[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]].mean()
# 	DistancesDF['stdv'] = DistancesDF[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]].std()
# 	DistancesDF['mean diff'] = DiffDF.mean(axis=1)
# 	DistancesDF['max diff'] = DiffDF.max(axis=1)
# 	DistancesDF.round({'mean':2,'std':2,'mean diff':2,'max diff':2})
# 	DistancesDF.to_csv(upl.replace('.upl','_distances_v1.csv'))
# 	print(DistancesDF)
# 	print(DistancesDF.shape)


### Now do this replacing protons with heavy atoms but drop the viaolation counts and so on 


pd.set_option("display.precision", 2)
DistancesDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv'])
DiffDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
print('Extracting Distances from %s' %fupl)
## Extract distances for the final upl no modification, and excluding ambiguous restraints
for line in open(fupl).readlines():
	if 'SUP' not in line:
		pass 
	else:
		cns = line.strip().split()
		atom1 = '{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2])
		atom2 = '{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],cns[5])
		if cns[1]+cns[2] in replacements.keys():
			atom1 = '{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],replacements[cns[1]+cns[2]])
		if cns[4]+cns[5] in replacements.keys():
			atom2 = '{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],replacements[cns[4]+cns[5]])
		connection = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
		DistancesDF.loc[connection,'upl'] = cns[6]
		diff = []
		for mnum in range(1,21,1):
			Coor = eval('Coor' + str(mnum))
			d = getDistance(atom1, atom2, Coor)
			DistancesDF.loc[connection,mnum] = d

DistancesDF['mean'] = DistancesDF[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]].mean(axis=1)
DistancesDF['stdv'] = DistancesDF[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]].std(axis=1)
decimals = pd.Series([2,2], index=['mean','stdv'])
DistancesDF.round(2)
DistancesDF.to_csv(pdb_name + '_distances_heavy_atom.csv')
print(DistancesDF)
print(DistancesDF.shape)

# ## Now look at the other input upls 
# for upl in upls:
# 	print(upl)
# 	DiffDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
# 	DistancesDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv'])
# 	for line in open(cwd+upl).readlines():
# 		if '#'in line[0:4]:
# 			pass 
# 		else:
# 			cns = line.strip().split()
# 			atom1 = '{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2])
# 			atom2 = '{:}{:}-{:}'.format(AAA_dict[cns[4]],cns[3],cns[5])
# 			connection = '{:}{:}-{:}-{:}{:}-{:}'.format(AAA_dict[cns[1]],cns[0],cns[2],AAA_dict[cns[4]],cns[3],cns[5])
# 			DistancesDF.loc[connection,'upl'] = float(cns[6])
# 			for mnum in range(1,21,1):
# 				Coor = eval('Coor' + str(mnum))
# 				d = getDistance(atom1, atom2,Coor)
# 				DistancesDF.loc[atom1 + ' ' + atom2,mnum] = d
# 				if np.round(d - float(cns[6]),2) > 0.1:
# 					DiffDF.loc[connection,mnum] = np.round(d - float(cns[6]),2)
# 			if connection in violdict.keys():
# 				DistancesDF.loc[connection,'violation'] = violdict[connection]
# 	DistancesDF['mean'] = np.round(DistancesDF[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]].mean(axis=1),2)
# 	DistancesDF['stdv'] = np.round(DistancesDF[[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]].std(axis=1),2)
# 	DistancesDF['mean diff'] = np.round(DiffDF.mean(axis=1),2)
# 	DistancesDF['max diff'] = np.round(DiffDF.max(axis=1),2)
# 	DistancesDF.to_csv(upl.replace('.upl','_distances_v1.csv'))
# 	print(DistancesDF)
# 	print(DistancesDF.shape)


