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
DistancesDF = pd.DataFrame()


in_pdb = '/Users/mclay1/FGFR3_structure/cyana_23/final.pdb'
in_upl = '/Users/mclay1/FGFR3_structure/cyana_23/final.upl'
chains = 'A'
Coords = {}
Sequence = []
# PDB_df = pd.DataFrame(columns=pdb_columns)
Starts = []
Ends = []
for mnum in range(1,21,1):
	start = open(in_pdb).readlines().index('MODEL' + '%9s\n' %str(mnum))
	exec('Coor' + str(mnum) + ' = {}')
	Starts.append(start)
	Ends.append(start-1)
Ends.append(open(in_pdb).readlines().index('END\n'))
Ends = Ends[1:]
pdb = open(in_pdb).readlines()
n = 0
for (start,end) in zip(Starts,Ends):
	n+=1
	Coor = eval('Coor' + str(n))
	for x in range(start,end,1):
		line = pdb[x]
		if line[0:4] == "ATOM" or line[0:4] == 'HETA':
			index =  line[17:20].strip() + '-' + line[22:26].strip() + '-' + line[12:16].strip()
			Coor[index] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
			if line[17:20].strip() in ['ILE','LEU','ALA','VAL','MET','THR'] and line[12:16].strip() in ['HB1','HE1','HD11','HD21','HG11','HG21']:
				index =  line[17:20].strip() + '-' + line[22:26].strip() + '-' + line[12:16].strip()[:-1].replace('H','Q')
				print(line[12:16].strip())
				line2 = pdb[x+1]
				line3 = pdb[x+2]
				print(line2[12:16].strip())
				print(line3[12:16].strip())
				xcoor = (float(line[30:38]) + float(line2[30:38]) + float(line3[30:38]))/3.0
				ycoor = (float(line[38:46]) + float(line2[38:46]) + float(line3[38:46]))/3.0
				zcoor = (float(line[46:54]) + float(line2[46:54]) + float(line3[46:54]))/3.0
				Coor[index] = [xcoor,ycoor,zcoor]
			if line[17:20].strip() not in ['ILE','LEU','ALA','VAL','MET','THR'] and len(line[12:16].strip()) == 4:
				index =  line[17:20].strip() + '-' + line[22:26].strip() + '-' + 'Q' + line[13:16].strip()[:-1]
				line2 = pdb[x+1]
				xcoor = (float(line[30:38]) + float(line2[30:38]))/2.0
				ycoor = (float(line[38:46]) + float(line2[38:46]))/2.0
				zcoor = (float(line[46:54]) + float(line2[46:54]))/2.0
				Coor[index] = [xcoor,ycoor,zcoor]
			if line[17:20].strip() in ['PHE','TYR'] and line[12:16].strip() == 'HB2':
				index =  line[17:20].strip() + '-' + line[22:26].strip() + '-' + 'QB'
				line2 = pdb[x+1]
				xcoor = (float(line[30:38]) + float(line2[30:38]))/2.0
				ycoor = (float(line[38:46]) + float(line2[38:46]))/2.0
				zcoor = (float(line[46:54]) + float(line2[46:54]))/2.0
				Coor[index] = [xcoor,ycoor,zcoor]
for line in open(in_upl).readlines():
	if 'SUP' not in line:
		pass 
	else:
		cns = line.strip().split()
		atom1 = cns[1].replace('HIST','HIS') + '-' + cns[0] + '-' + cns[2]
		atom2 = cns[4].replace('HIST','HIS') + '-' + cns[3] + '-' + cns[5]
		for mnum in range(1,21,1):
			Coor = eval('Coor' + str(mnum))
			(x1,y1,z1) = Coor[atom1]
			(x2,y2,z2) = Coor[atom2]
			d = np.round(np.sqrt((x1-x2)**2 + (y1-y2)**2+ (z1-z2)**2),2)
			DistancesDF.loc[atom1 + '-' + atom2,mnum] = d
	DistancesDF['mean'] = np.round(DistancesDF.mean(axis=1),2)
for line in open(in_upl).readlines():
	if 'SUP' not in line:
		pass 
	else:
		cns = line.strip().split()
		atom1 = cns[1].replace('HIST','HIS') + '-' + cns[0] + '-' + cns[2]
		atom2 = cns[4].replace('HIST','HIS') + '-' + cns[3] + '-' + cns[5]
	DistancesDF.loc[atom1 + '-' + atom2,'upl'] = cns[6]
DistancesDF.to_csv('test.csf')
print(DistancesDF)
print(DistancesDF.shape)
