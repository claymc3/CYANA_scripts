from math import sqrt, cos, sin, acos, pi
import pandas as pd 
import os
import numpy as np
import sys
pdb_columns = ['name', 'resn', 'resid', 'Chain', 'X', 'Y', 'Z', 'Element']
# Read in PDB file one line at a time, if the first four letter ar ATOM or HETA then it will parse the data into the 
# data frame, using the atome index int he PDB as the row index in the data frame. 
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V' }
A_dict = { "I": "ILE", "A": "ALA", "V": "VAL", "L": "LEU", "M": "MET", "T" : "THR","F" : "PHE", "Y": "TYR" }
Methyl_groups = {"LEU" : ['CD1','CD2'], "ILE" : ['CD1'], "VAL" : ['CG1','CG2'], "THR" : ['CG2'], "MET" : ['CE'], "ALA" : ['CB']}
Methyls = ['ILE', 'LEU','VAL','MET','ALA', 'THR']
Aromatics = ['PHE', 'TYR']
Aromatic_groups = {"PHE" : ['CE1', 'CE2'], "TYR": ['CE1','CHE2']}
Pass_Atoms = ['N   ALA', 'H   ALA', 'N   ARG', 'H   ARG', 'N   ASN', 'H   ASN', 'N   ASP', 'H   ASP', 'N   CYS', 'H   CYS', 'N   GLU', 'H   GLU', 'N   GLN', 'H   GLN', 'N   GLY', 'H   GLY', 'N   HIS', 'H   HIS', 'N   ILE', 'H   ILE', 'N   LEU', 'H   LEU', 'N   LYS', 'H   LYS', 'N   MET', 'H   MET', 'N   PHE', 'H   PHE', 'N   SER', 'H   SER', 'N   THR', 'H   THR', 'N   TRP', 'H   TRP', 'N   TYR', 'H   TYR', 'N   VAL', 'H   VAL']


if len(sys.argv) == 1:
	print('Usage: GetDieh pdb chain')
	exit()

def crossProduct(u,v):
	## Calculates the cross product of two 3d vectors (as 1-d arrays).
	prod = [0.,0.,0.]
	prod[0] = u[1]*v[2]-u[2]*v[1]
	prod[1] = u[2]*v[0]-u[0]*v[2]
	prod[2] = u[0]*v[1]-u[1]*v[0]
	return prod

def norm(v):
	mag = sqrt(v[0]**2 + v[1]**2 + v[2]**2)
	normv = [v[i]/mag for i in range(3)]
	return normv

def getOrthoNorm(v,ref):
	temp = [v[i] - ref[i]* dotProduct(v,ref) for i in range(3)]
	orth = norm(temp)
	return orth

def dotProduct(u,v):
	## Calculates the dot product between two vectors.
	return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]


def calcDihedrals(A,B,C,D):
	##Calculates dihedral angle side or main chain

	U2 = norm([C[i] - B[i] for i in range(3)])
	U1 = getOrthoNorm([A[i] - B[i] for i in range(3)], U2)
	U3 = getOrthoNorm([D[i] - C[i] for i in range(3)], U2)

	angle = 180/pi * acos(dotProduct(U1,U3))
	sign = dotProduct(crossProduct(U1,U3),U2)

	if sign < 0:
		angle = -angle
	return angle


# phi main chain torsion angle for atoms C-1,N,CA,C
# psi main chain torsion angle for atoms N,CA,C,N+1
# chi1 side chain torsion angle for atoms N,CA,CB,*G
# chi2 side chain torsion angle for atoms CA,CB,*G,*D
# chi3 side chain torsion angle for atoms CB,*G,*D,*E
# chi4 side chain torsion angle for atoms *G,*D,*E,*Z
# chi5 side chain torsion angle for atoms *D,*E,*Z, NH1
# return phi, psi
#MODEL        1
#

SideDihe = {'R':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD' ], ['chi3', 'CB', 'CG', 'CD', 'NE'], ['chi4', 'CG', 'CD', 'NE', 'CZ'], ['chi5', 'CD', 'NE', 'CZ' ,'NH1']],
 'N':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'OD1']],
 'D':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'OD1']],
 'C':[['chi1', 'N', 'CA', 'CB','SG']],
 'Q':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD'], ['chi3','CB', 'CG', 'CD', 'OE1']],
 'E':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD'], ['chi3', 'CB', 'CG', 'CD', 'OE1']],
 'H':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2','CA', 'CB', 'CG', 'ND1']],
 'I':[['chi1', 'N', 'CA', 'CB','CG1'], ['chi2', 'CA', 'CB', 'CG1', 'CD1']],
 'L':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
 'K':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD'],['chi3', 'CB', 'CG', 'CD', 'CE'],['chi4', 'CG', 'CD', 'CE', 'NZ']],
 'M':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'SD'],['chi3', 'CB', 'CG', 'SD', 'CE']],
 'F':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
 'P':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD']],
 'S':[['chi1', 'N', 'CA', 'CB','OG']],
 'T':[['chi1', 'N', 'CA', 'CB','OG1']],
 'W':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
 'Y':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
 'V':[['chi1', 'N', 'CA', 'CB','CG1']]}



Dihedrals = pd.DataFrame(columns=['phi','psi','chi1','chi2'])
Phidihed = pd.DataFrame()
Psidided = pd.DataFrame()
Chi1dihed = pd.DataFrame()
Chi2dihed = pd.DataFrame()

### XXX#-atom:[x,y,z]
in_pdb = sys.argv[1]
chains = sys.argv[2]
# in_pdb = '4k33.pdb'
# chains = 'A'
Coords = {}
Sequence = []
PDB_df = pd.DataFrame(columns=pdb_columns)
Starts = []

if 'MODEL        1' not in open(in_pdb).read():
	for line in open(in_pdb):
		if line[0:4] == "ATOM" or line[0:4] == 'HETA':
			if  line[17:20].strip() in AAA_dict.keys():
				index =  AAA_dict[line[17:20].strip()] + line[22:26].strip() + '-' + line[12:16].strip()
				if AAA_dict[line[17:20].strip()] + line[22:26].strip() not in Sequence:
					Sequence.append(AAA_dict[line[17:20].strip()] + line[22:26].strip())
				Coords[index] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

	for i in range(1,len(Sequence)-1,1):
		phi = calcDihedrals(Coords[Sequence[i-1]+ '-C'],Coords[Sequence[i]+ '-N'],Coords[Sequence[i]+ '-CA'],Coords[Sequence[i]+ '-C'])
		psi = calcDihedrals(Coords[Sequence[i]+ '-N'],Coords[Sequence[i]+ '-CA'],Coords[Sequence[i]+ '-C'],Coords[Sequence[i+1]+ '-N'])
		Dihedrals.loc[Sequence[i].split('-')[0],'phi'] = np.round(phi,1)
		Dihedrals.loc[Sequence[i].split('-')[0],'psi'] = np.round(psi,1)
	for res in Sequence: 
		if res[0] in SideDihe.keys():
			for dihe in SideDihe[res[0]]:
				if res+ '-' + dihe[1] in Coords.keys() and res+ '-' + dihe[2] in Coords.keys() and res+ '-' + dihe[3] in Coords.keys() and res+ '-' + dihe[4] in Coords.keys():
					ang = calcDihedrals(Coords[res+ '-' + dihe[1]],Coords[res+ '-' + dihe[2]],Coords[res+ '-' + dihe[3]],Coords[res+ '-' + dihe[4]])
					if ang < 0: ang = ang + 360.0
					Dihedrals.loc[res.split('-')[0],dihe[0]] = np.round(ang,1)
	Dihedrals.to_csv(in_pdb.replace('.pdb', '_dihedrals.csv'))

# if 'MODEL        1' not in open(in_pdb).read():
# 	for mnum in range(1,21,1):
# 		start = open(in_pdb).readlines().index('MODEL' + '%9s\n' %str(mnum))
# 		exec('Coor' + str(mnum) + ' = {}')
# 		Coor = eval('Coor' + str(mnum))
# 		Starts.append(start)
# 		for line in open(in_pdb).readlines()[start:]:
# 			if line == 'ENDMDL\n': break
# 			if line[0:4] == "ATOM" or line[0:4] == 'HETA':
# 				index =  line[17:20].strip() + line[22:26].strip() + '-' + line[12:16].strip()
# 				if line[17:20].strip() + line[22:26].strip() not in Sequence:
# 					Sequence.append(line[17:20].strip() + line[22:26].strip())
# 				Coords[index] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]

# 	print(Coords)
# 	for i in range(1,len(Sequence)-1,1):
# 	    phi = calcDihedrals(Coords[Sequence[i-1]+ '-C'],Coords[Sequence[i]+ '-N'],Coords[Sequence[i]+ '-CA'],Coords[Sequence[i]+ '-C'])
# 	    psi = calcDihedrals(Coords[Sequence[i]+ '-N'],Coords[Sequence[i]+ '-CA'],Coords[Sequence[i]+ '-C'],Coords[Sequence[i+1]+ '-N'])
# 	    Dihedrals.loc[Sequence[i].split('-')[0],'phi'] = np.round(phi,1)
# 	    Dihedrals.loc[Sequence[i].split('-')[0],'psi'] = np.round(psi,1)
# 	for res in Sequence: 
# 	    if res[0] in SideDihe.keys():
# 	        for dihe in SideDihe[res[0]]:
# 	            if res+ '-' + dihe[1] in Coords.keys() and res+ '-' + dihe[2] in Coords.keys() and res+ '-' + dihe[3] in Coords.keys() and res+ '-' + dihe[4] in Coords.keys():
# 	                ang = calcDihedrals(Coords[res+ '-' + dihe[1]],Coords[res+ '-' + dihe[2]],Coords[res+ '-' + dihe[3]],Coords[res+ '-' + dihe[4]])
# 	                Dihedrals.loc[res.split('-')[0],dihe[0]] = np.round(ang,1)

print(Dihedrals)
