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


in_pdb = '/Volumes/common/Kinases/FGFR3/FGFR3_459-755_C482A_C582S/Structure_Calc/cyana_15/CNS/refinedPDB/FGFR3_r15_cya.pdb'
in_upl = '/Volumes/common/Kinases/FGFR3/FGFR3_459-755_C482A_C582S/Structure_Calc/cyana_15/CNS/'
chains = 'A'
Coords = {}
Sequence = []
PDB_df = pd.DataFrame(columns=pdb_columns)
Starts = []
for mnum in range(1,21,1):
	start = open(in_pdb).readlines().index('MODEL' + '%9s\n' %str(mnum))
	exec('Coor' + str(mnum) + ' = {}')
	Coor = eval('Coor' + str(mnum))
	Starts.append(start)
	for line in open(in_pdb).readlines()[start:]:
		if line == 'ENDMDL\n': break
		if line[0:4] == "ATOM" or line[0:4] == 'HETA':
			index =  line[17:20].strip() + line[22:26].strip() + '-' + line[12:16].strip()
			if line[17:20].strip() + line[22:26].strip() not in Sequence:
				Sequence.append(line[17:20].strip() + line[22:26].strip())
			Coords[index] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
DistancesDF = pd.DataFrame(columns=['m1','m2','m3','m4','m5','m6','m7','m8','m9','m10','m11','m12','m13','m14','m15','m16','m17','m18','m19','m20'])



