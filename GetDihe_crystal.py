from math import sqrt, cos, sin, acos, pi
import pandas as pd 
import os
import re
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
from math import ceil, floor
from rama_config import RAMA_PREFERENCES
from rama_config import ROTA_PREFERENCES
import requests
mpl.rcParams['pdf.fonttype'] = 42
if sys.platform == "linux":
	mpl.rcParams['font.sans-serif'] = 'Dejavu Sans'
else:
	mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.linewidth'] = 1.0
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['xtick.labelsize'] = mpl.rcParams['ytick.labelsize']=8
plt.rcParams['mathtext.default'] = 'regular'

pdb_columns = ['name', 'resn', 'resid', 'Chain', 'X', 'Y', 'Z', 'Element']
# Read in PDB file one line at a time, if the first four letter ar ATOM or HETA then it will parse the data into the 
# data frame, using the atome index int he PDB as the row index in the data frame. 
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",
 "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S','CYSS':'C', 'HIST':'H'}
A_dict = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS', 'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN','G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP','A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

RAMA_PREF_VALUES = None
ROTA_PREF_VALUES = None

def _cache_RAMA_PREF_VALUES():
	## Values in data span -179 to 179 in incrments of 2, so values are translated to 0 to 360 in increments of 1
	f_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
	RAMA_PREF_VALUES = {}
	for key, val in RAMA_PREFERENCES.items():
		RAMA_PREF_VALUES[key] = np.full((361, 361), 0, dtype=np.float64)
		with open(os.path.join(f_path, val["file"])) as fn:
			for line in fn:
				if not re.match('^\s*#', line):
					x = int(float(line.split()[1]))
					y = int(float(line.split()[0]))
					for nx in np.arange(x-1,x+1,1):
						for ny in np.arange(y-1,y+1,1):
							RAMA_PREF_VALUES[key][nx + 180][ny + 180] = float(line.split()[-1])
	return RAMA_PREF_VALUES

def _cache_ROTA_PREF_VALUES():
	f_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
	ROTA_PREF_VALUES = {}
	for key, val in ROTA_PREFERENCES.items():
		ROTA_PREF_VALUES[key] = np.full((361, 361), 0, dtype=np.float64)
		with open(os.path.join(f_path, val["file"])) as fn:
			for line in fn:
				if not re.match('^\s*#', line):
					x = int(float(line.split()[1]))
					y = int(float(line.split()[0]))
					for nx in np.arange(x-1,x+1,1):
						for ny in np.arange(y-1,y+1,1):
							ROTA_PREF_VALUES[key][nx][ny] = float(line.split()[-1])
							if key in ['F','Y','D']:
								ROTA_PREF_VALUES[key][nx+180][ny] = float(line.split()[-1])
	return ROTA_PREF_VALUES

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
	angle = 180/pi * acos(np.round(dotProduct(U1,U3),3))
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

#SideDihe = {
# 'R':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD' ], ['chi3', 'CB', 'CG', 'CD', 'NE'], ['chi4', 'CG', 'CD', 'NE', 'CZ'], ['chi5', 'CD', 'NE', 'CZ' ,'NH1']],
# 'N':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'OD1']],
# 'D':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'OD1']],
# 'C':[['chi1', 'N', 'CA', 'CB','SG']],
# 'Q':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD'], ['chi3','CB', 'CG', 'CD', 'OE1']],
# 'E':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD'], ['chi3', 'CB', 'CG', 'CD', 'OE1']],
# 'H':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2','CA', 'CB', 'CG', 'ND1']],
# 'I':[['chi1', 'N', 'CA', 'CB','CG1'], ['chi2', 'CA', 'CB', 'CG1', 'CD1']],
# 'L':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
# 'K':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD'],['chi3', 'CB', 'CG', 'CD', 'CE'],['chi4', 'CG', 'CD', 'CE', 'NZ']],
# 'M':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'SD'],['chi3', 'CB', 'CG', 'SD', 'CE']],
# 'F':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
# 'P':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD']],
# 'S':[['chi1', 'N', 'CA', 'CB','OG']],
# 'T':[['chi1', 'N', 'CA', 'CB','OG1']],
# 'W':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
# 'Y':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
# 'V':[['chi1', 'N', 'CA', 'CB','CG1']]}
SideDihe = {
 'R':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD' ]],
 'N':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'OD1']],
 'D':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'OD1']],
 'Q':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD']],
 'E':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD']],
 'H':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2','CA', 'CB', 'CG', 'ND1']],
 'I':[['chi1', 'N', 'CA', 'CB','CG1'], ['chi2', 'CA', 'CB', 'CG1', 'CD1']],
 'K':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD']],
 'L':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
 'M':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'SD']],
 'F':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
 'Y':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']],
 'W':[['chi1', 'N', 'CA', 'CB','CG'], ['chi2', 'CA', 'CB', 'CG', 'CD1']]}
def plot_phi_psi_ramachandran(res, ax, PhiDF, PsiDF,axtext,ypos,colnames,plotcolors,logfile):

	global RAMA_PREF_VALUES

	if RAMA_PREF_VALUES is None:
		RAMA_PREF_VALUES = _cache_RAMA_PREF_VALUES()
	outtext,plotc = [],[]
	normals = {"phi":[],"psi":[]}
	outliers = {"phi":[],"psi":[]}
	aa_type = PhiDF.loc[res,'type']
	if pd.isna(PhiDF.loc[res,'type']): aa_type = 'General'
	outline = "   "
	outcount = 0
	outtext.append([r"$\phi$ = {:} $\pm$ {:}".format(PhiDF.loc[res,'mean'],PhiDF.loc[res,'stdv']),'black'])
	outtext.append([r"$\psi$ = {:} $\pm$ {:}".format(PsiDF.loc[res,'mean'],PsiDF.loc[res,'stdv']),'black'])
	for mnum in colnames:
		if not np.isnan(PsiDF.loc[res,mnum]) and not np.isnan(PhiDF.loc[res,mnum]):
			# ax.scatter(normals["phi"], normals["psi"],marker='o',s= 30,facecolors=plotcolors[mnum], edgecolors= 'none', linewidth=1.0)
			if RAMA_PREF_VALUES[aa_type][int(PsiDF.loc[res,mnum])+ 180][int(PhiDF.loc[res,mnum]) + 180] < RAMA_PREFERENCES[aa_type]["bounds"][1]:
				outliers["phi"].append(PhiDF.loc[res,mnum])
				outliers["psi"].append(PsiDF.loc[res,mnum])
				outline = outline  + "{:} ".format(mnum, PhiDF.loc[res,mnum],PsiDF.loc[res,mnum])
				outcount+=1
			# else:
			normals["phi"].append(PhiDF.loc[res,mnum])
			normals["psi"].append(PsiDF.loc[res,mnum])
			plotc.append(plotcolors[mnum])
	if outcount != 0:
		logfile.write("{:} Phi, Psi disallowed in:\n".format(res))
		logfile.write('  {:}\n'.format(outline))
	ax.imshow(RAMA_PREF_VALUES[aa_type], cmap=RAMA_PREFERENCES[aa_type]["cmap"],
			norm=colors.BoundaryNorm(RAMA_PREFERENCES[aa_type]["bounds"], RAMA_PREFERENCES[aa_type]["cmap"].N),
			extent=(-180, 180, 180, -180))
	ax.scatter(normals["phi"], normals["psi"],marker='o',s= 30,c=plotc, edgecolors= 'none', linewidth=1.0)
	ax.set_xlabel(r'$\mathrm{\phi}$')
	ax.set_ylabel(r'$\mathrm{\psi}$')
	tcolor = 'black'
	title = res
	ax.set_title(title, color = tcolor)
	if len(outtext) > 0:
		for text, col in outtext:
			axtext.text(-0.3,ypos, text, color = col, fontsize = 7)
			ypos = ypos - 0.06
	ax.set_xlim([-180,180])
	ax.set_xticks([-180,-120,-60,0,60,120,180])
	ax.set_yticks([-180,-120,-60,0,60,120,180])
	ax.set_ylim([-180,180])
	ax.grid(visible=True, which='major', axis='both',linestyle='--')
	# plt.tight_layout(w_pad = 0.0001)

def plot_chi1_chi2_ramachandran(res, ax, chi1DF, chi2DF, axtext, ypos,colnames,plotcolors,logfile):

	global ROTA_PREF_VALUES

	if ROTA_PREF_VALUES is None:
		ROTA_PREF_VALUES = _cache_ROTA_PREF_VALUES()
	outtext,plotc = [],[]
	normals = {"chi1":[],"chi2":[]}
	outliers = {"chi1":[],"chi2":[]}
	outline = "   "
	aa_type = res[0]
	outcount = 0
	outtext.append([r"$\chi1$ = {:} $\pm$ {:}".format(chi1DF.loc[res,'mean'],chi1DF.loc[res,'stdv']),'black'])
	outtext.append([r"$\chi2$ = {:} $\pm$ {:}".format(chi2DF.loc[res,'mean'],chi2DF.loc[res,'stdv']),'black'])
	for mnum in colnames:
		if not np.isnan(chi1DF.loc[res,mnum]) and not np.isnan(chi2DF.loc[res,mnum]):
			if ROTA_PREF_VALUES[aa_type][int(chi2DF.loc[res,mnum])][int(chi1DF.loc[res,mnum])] < ROTA_PREFERENCES[aa_type]["bounds"][1]:
				outliers["chi1"].append(chi1DF.loc[res,mnum])
				outliers["chi2"].append(chi2DF.loc[res,mnum])
				outcount+=1
				outline = outline  + "{:} ".format(mnum)
			# else:
			normals["chi1"].append(chi1DF.loc[res,mnum])
			normals["chi2"].append(chi2DF.loc[res,mnum])
			plotc.append(plotcolors[mnum])
	if outcount != 0: 
		logfile.write("{:} chi1, chi2 disallowed in:\n".format(res))
		logfile.write('  {:}\n'.format(outline))
	ax.imshow(ROTA_PREF_VALUES[aa_type], cmap=ROTA_PREFERENCES[aa_type]["cmap"],
			norm=colors.BoundaryNorm(ROTA_PREFERENCES[aa_type]["bounds"], ROTA_PREFERENCES[aa_type]["cmap"].N),
			extent=(0, 360, 360, 0))
	ax.scatter(normals["chi1"], normals["chi2"],marker='o',s= 30,c=plotc, edgecolors= 'none', linewidth=1.0)
	ax.set_xlabel(r'$\mathrm{\chi}1$')
	ax.set_ylabel(r'$\mathrm{\chi}2$')
	tcolor = 'black'
	title = res
	ax.set_title(title, color = tcolor)
	if len(outtext) > 0:
		for text, col in outtext:
			axtext.text(-0.3,ypos, text, color = col, fontsize = 7)
			ypos = ypos - 0.06
	ax.set_xlim([0,360])
	ax.set_xticks([0,60,120,180,240,300,360])
	ax.set_yticks([0,60,120,180,240,300,360])
	ax.set_ylim([0,360])
	ax.grid(visible=True, which='major', axis='both',linestyle='--')
	# plt.tight_layout(w_pad = 0.0001)


###----------------------------------------------------------------------------###
##	Extract the starting and ending lines for the 20 NMR models in the file.	##
##	Build a separate dictionary for each set of coordinates with the format:	##
##								A#-atom:[x,y,z] 								##
##	where A# is the residue single letter code and index 						##
###----------------------------------------------------------------------------###
if len(sys.argv)==1:
	print('''
Usage:
	getdihe [input_file]

	input_file: text file containing the following on separate lines

	Example :	getdihe FGFR1_active_input.txt

	Example input file (explanation below):
		FGFR1_Active_pdbs
		P11362
		464-762
		3GQI green
		5FLF blue
		3KXX orange

	outname:	name to give output files, no spaces allowed

	UniProt_id: UniProt accession id number for protein, the script 
				will pull that fasts file from https://www.uniprot.org
				and use only chains which match this UniProt ID 

	residues:	start-end indices of the residues you are interested in 468-768

	PDB ID color: the PDB id space and color to use for this PDB on the plots

	The script will extract all the possible phi,psi and chi1/chi2
	dihedral angles. It also checks the protein sequence from the 
	PDB and identifies mutation and phosphorylation sites, and identifies
	missing residues or side chains in the PDB. These are recorded in the 
	outname_log.txt file, along with the identity of any disallowed dihderals

	Dependencies:pandas, numpy, matplotlib, requests
		''')
	exit()

infile = open(sys.argv[1]).readlines()
# infile = open('FGFR1_structures_input.txt').readlines()

outname = infile[0].strip()
uniprotid = infile[1].strip()
seqbounds = infile[2].strip()
logfile =  open("{:}_log.txt".format(outname),'w')
inpdbs,pdbnames, colorsd = [], [], {}
x=-1
for line in infile[3:]:
	if line.strip():
		inpdbs.append(line.split()[0])
		pdbname=line.split()[0].split('/')[-1].split('.')[0].replace('-','_')
		pdbnames.append(pdbname)
		colorsd[pdbname]= line.split()[1].strip()
# print(inpdbs)
# print(seqbounds)
# exit()
UNPdict,UNPnseq = {},[]
UNPseq = ''
dbsb = int(seqbounds.split('-')[0])
dbse = int(seqbounds.split('-')[1])+1
for sline in requests.get('https://www.uniprot.org/uniprot/'+uniprotid+'.fasta', allow_redirects=True).iter_lines(decode_unicode=True):
	if sline.rstrip() and ">" not in sline and "#" not in sline:
		UNPseq = UNPseq + sline.strip()
UNPseq = UNPseq[dbsb - 1:int(seqbounds.split('-')[1])]
index = dbsb
for res in UNPseq:
	UNPdict[str(index)] = res
	UNPnseq.append(res+str(index))
	index+=1
plotcolors,colnames = {},[]
MasterDict ={}
for in_pdb in inpdbs:
	pdbname = in_pdb.split('/')[-1].split('.')[0].replace('-','_')
	print('Extracting dihedrals for: {:}'.format(pdbname))
	chains = []
	pSEQ_dict,pSEQ = {},[]
	MasterDict['Seq_{:}'.format(pdbname)] = {}
	SEQdict = MasterDict['Seq_{:}'.format(pdbname)]
	if not re.search('.pdb',in_pdb):
		pdblines = list(requests.get('https://files.rcsb.org/view/'+in_pdb+'.pdb', allow_redirects=True).iter_lines(decode_unicode=True))
		for line in pdblines:
			if line[0:5] == 'DBREF':
				name = line[33:41].strip()
				if name == uniprotid:
					chains.append(line[12])
			if line[0:5] == 'ATOM ':
				break
	if re.search('.pdb',in_pdb):
		pdblines = open(in_pdb).readlines()
		# if 'DBREF' not in open(in_pdb).read():
		for line in pdblines:
			if line[0:5] == 'ATOM ' and line[21] not in chains:
				chains.append(line[21])
	for chain in chains:
		MasterDict['Res_{:}_{:}'.format(pdbname,chain)] = {}
		MasterDict['Coor_{:}_{:}'.format(pdbname,chain)] = {}
		colnames.append('{:}_{:}'.format(pdbname,chain))
		plotcolors['{:}_{:}'.format(pdbname,chain)] = colorsd[pdbname]
	for line in pdblines:
		if line[0:4] == "ATOM" or line[0:4] == 'HETA':
			if line[21] in chains:
				chain = line[21]
				Coor = MasterDict['Coor_{:}_{:}'.format(pdbname,chain)]
				Resdict = MasterDict['Res_{:}_{:}'.format(pdbname,chain)]
				if line[17:20].strip() in AAA_dict.keys():
					if line[22:26].strip() not in pSEQ: 
						pSEQ.append(line[22:26].strip())
						pSEQ_dict[line[22:26].strip()] = line[17:20].strip()
						SEQdict[int(line[22:26].strip())] = '{:}{:}'.format(AAA_dict[line[17:20].strip()],line[22:26].strip())
					resid = '{:}{:}-{:}'.format(AAA_dict[line[17:20].strip()], line[22:26].strip(),line[12:16].strip())
					Coor[resid]= [float(line[30:38]), float(line[38:46]), float(line[46:54])]
					if line[22:26].strip() in Resdict.keys() and line[12:16].strip()[0] != 'H':
						Resdict[line[22:26].strip()].append(line[12:16].strip())
					elif line[22:26].strip() not in Resdict.keys() and line[12:16].strip()[0] != 'H':
						Resdict[line[22:26].strip()]=[line[12:16].strip()]

	logfile.write('{:}\n'.format(in_pdb.upper()))
	for x in range(dbsb,dbse):
		if str(x) in pSEQ_dict.keys():
			if pSEQ_dict[str(x)] in ["PTR", "TPO", "SEP"]: logfile.write('  Phosphorlated p{:}{:}\n'.format(UNPdict[str(x)],x))
			elif A_dict[UNPdict[str(x)]] != pSEQ_dict[str(x)]:
				logfile.write('  mutation {:}{:}{:}\n'.format(UNPdict[str(x)],x,AAA_dict[pSEQ_dict[str(x)]]))
	for chain in chains:
		missingSide,missingRes ='',''
		Resdict = MasterDict['Res_{:}_{:}'.format(pdbname,chain)]
		for x in range(dbsb,dbse):
			if str(x) in Resdict.keys():
				if pSEQ_dict[str(x)] not in ['ALA','GLY'] and len(Resdict[str(x)]) < 6:
					missingSide = missingSide + '{:}{:} '.format(AAA_dict[pSEQ_dict[str(x)]],x)
			if str(x) not in Resdict.keys():
				missingRes= missingRes + '{:}{:} '.format(UNPdict[str(x)],x)
		if len(missingSide) != 0:
			logfile.write('  {:} {:} Missing side chains:\n    {:}\n'.format(in_pdb.upper(),chain,missingSide))
		if len(missingRes) != 0:
			logfile.write('  {:} {:} Missing Residues:\n    {:}\n'.format(in_pdb.upper(),chain,missingRes))

# colnames.extend(['mean','stdv'])
PhiDF =  pd.DataFrame(columns=colnames)
MasterDict['PhiDF'] = PhiDF
PsiDF =  pd.DataFrame(columns=colnames)
MasterDict['PsiDF'] = PsiDF
chi1DF = pd.DataFrame(columns=colnames)
MasterDict['Chi1DF'] = Chi1DF
chi2DF = pd.DataFrame(columns=colnames)
MasterDict['Chi2DF'] = Chi2DF
dihedDF = pd.DataFrame(columns=colnames)
for entry in colnames:
	Coords = MasterDict['Coor_' + entry]
	Sequence = MasterDict['Seq_{:}'.format(entry[0:-2])]
	for i in range(dbsb+1,dbse,1):
		if i-1 in Sequence.keys() and i in Sequence.keys() and i+1 in Sequence.keys():
			if Sequence[i]+ '-N' in Coords.keys() and Sequence[i]+ '-C' in Coords.keys() and Sequence[i]+ '-CA' in Coords.keys() and Sequence[i+1]+ '-N' in Coords.keys() and Sequence[i-1]+ '-C' in Coords.keys():
				phi = calcDihedrals(Coords[Sequence[i-1]+ '-C'],Coords[Sequence[i]+ '-N'],Coords[Sequence[i]+ '-CA'],Coords[Sequence[i]+ '-C'])
				PhiDF.loc[Sequence[i],entry] = np.round(phi,1)
				psi = calcDihedrals(Coords[Sequence[i]+ '-N'],Coords[Sequence[i]+ '-CA'],Coords[Sequence[i]+ '-C'],Coords[Sequence[i+1]+ '-N'])
				PsiDF.loc[Sequence[i],entry] = np.round(psi,1)
	for res in UNPnseq: 
		if res[0] in SideDihe.keys():
			for dihe in SideDihe[res[0]]:
				diheDF = MasterDict[dihe[0] + 'DF']
				if res+ '-' + dihe[1] in Coords.keys() and res+ '-' + dihe[2] in Coords.keys() and res+ '-' + dihe[3] in Coords.keys() and res+ '-' + dihe[4] in Coords.keys():
					ang = calcDihedrals(Coords[res+ '-' + dihe[1]],Coords[res+ '-' + dihe[2]],Coords[res+ '-' + dihe[3]],Coords[res+ '-' + dihe[4]])
					if ang < 0: ang = ang + 360.0
					diheDF.loc[res,entry] = np.round(ang,1)
for i in range(1,len(UNPnseq)-1,1):
	if UNPnseq[i] in PhiDF.index.to_list():
		if UNPnseq[i+1][0] == 'P' and UNPnseq[i][0] != 'G': PhiDF.loc[UNPnseq[i],'type'] = "PRE-PRO"
		elif UNPnseq[i][0] == 'P':PhiDF.loc[UNPnseq[i],'type'] = "PRO"
		elif UNPnseq[i][0] == 'G':PhiDF.loc[UNPnseq[i],'type'] = "GLY"
		else: PhiDF.loc[UNPnseq[i],'type'] = "General"
for res in UNPnseq: 
	if res[0] in SideDihe.keys() and res in chi1DF.index.to_list():
		chi1DF.loc[res,'type'] = res[0]
PhiDF['mean'] = PhiDF[colnames].mean(axis=1).astype(float).round(2)
PhiDF['stdv'] = PhiDF[colnames].std(axis=1).astype(float).round(2)
PsiDF['mean'] = PsiDF[colnames].mean(axis=1).astype(float).round(2)
PsiDF['stdv'] = PsiDF[colnames].std(axis=1).astype(float).round(2)
chi1DF['mean'] = chi1DF[colnames].mean(axis=1).astype(float).round(2)
chi1DF['stdv'] = chi1DF[colnames].std(axis=1).astype(float).round(2)
chi2DF['mean'] = chi2DF[colnames].mean(axis=1).astype(float).round(2)
chi2DF['stdv'] = chi2DF[colnames].std(axis=1).astype(float).round(2)

for res in UNPnseq:
	if res in PhiDF.index.to_list():
		dihedDF.loc['{:}_Phi'.format(res)] = PhiDF.loc[res].copy()
	if res in  PsiDF.index.to_list():
		dihedDF.loc['{:}_Psi'.format(res)] = PsiDF.loc[res].copy()
	if res in chi1DF.index.to_list():
		dihedDF.loc['{:}_Chi1'.format(res)] = chi1DF.loc[res].copy()
	if res in chi1DF.index.to_list():
		dihedDF.loc['{:}_Chi2'.format(res)] = chi2DF.loc[res].copy()
dihedDF['mean'] = dihedDF[colnames].mean(axis=1).astype(float).round(2)
dihedDF['stdv'] = dihedDF[colnames].std(axis=1).astype(float).round(2)
# print(dihedDF)
dihedDF.to_csv(outname+'_dihed.csv')
logfile.write('\n\n')

import time
start_time = time.time()
pdf = PdfPages('{:}_overview.pdf'.format(outname))
count = 0
from matplotlib.gridspec import GridSpec
for res in UNPnseq:
	count+=1
	if count in np.arange(1,len(UNPnseq),25):
		print('Plotting results for {:} ({:} of {:})'.format(res, count, len(UNPnseq)))
	if res in PhiDF.index.to_list() or res in chi1DF.index.to_list():
		if res in PhiDF.index.to_list() and res in chi1DF.index.to_list():
			fig = plt.figure(figsize=(9,3))
			gs = GridSpec(1,3,width_ratios = (1,1,1))
			ax1 = fig.add_subplot(gs[0])
			ax2 = fig.add_subplot(gs[1])
			ax0 = fig.add_subplot(gs[2])
			# fig, (ax1,ax2,ax3,ax0) =plt.subplots(1,4,figsize=(9,3), width_ratios = [6,6,3,3])
			plot_phi_psi_ramachandran(res, ax1, PhiDF, PsiDF,ax0, 0.90,colnames,plotcolors,logfile)
			plot_chi1_chi2_ramachandran(res, ax2, chi1DF, chi2DF,ax0, 0.78,colnames,plotcolors,logfile)
		if res in PhiDF.index.to_list() and res not in chi1DF.index.to_list():
			fig = plt.figure(figsize=(6,3))
			gs = GridSpec(1,2,width_ratios = (1,1))
			ax1 = fig.add_subplot(gs[0])
			ax0 = fig.add_subplot(gs[1])
			# fig, (ax1,ax2,ax0) =plt.subplots(1,3,figsize=(6,3), width_ratios = [2,1,1])
			plot_phi_psi_ramachandran(res, ax1, PhiDF, PsiDF,ax0, 0.90, colnames,plotcolors,logfile)
		if res not in PhiDF.index.to_list() and res in chi1DF.index.to_list():
			fig = plt.figure(figsize=(6,3))
			gs = GridSpec(1,2,width_ratios = (1,1))
			ax1 = fig.add_subplot(gs[0])
			ax0 = fig.add_subplot(gs[1])
			# fig, (ax1,ax2,ax0) =plt.subplots(1,3,figsize=(6,3),width_ratios = [2,1,1])
			plot_chi1_chi2_ramachandran(res, ax1, chi1DF, chi2DF,ax0, 0.90, colnames,plotcolors,logfile)
		ax0.axis('off')
		xpos = -0.5
		for x in range(0,len(inpdbs),10):
			i = x
			y = 0.66
			xpos = xpos + 0.2
			for j in range(10):
				ax0.text(xpos, y, pdbnames[i].upper(), color = colorsd[pdbnames[i]], fontsize = 8)
				i+=1
				y = y - 0.06
				if i== len(inpdbs): break
		plt.tight_layout()
		pdf.savefig(transparent=True)
		plt.close()
pdf.close()
logfile.close()

