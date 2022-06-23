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


PseudoAtoms = {'ALAQB':['HB1', 'HB2', 'HB3'], 'ARGQB':['HB2', 'HB3'], 'ARGQG':['HG2', 'HG3'], 'ARGQD':['HD2', 'HD3'], 'ARGQH1':['HH11', 'HH12'], 'ARGQH2':['HH21', 'HH22'], 'ASNQB':['HB2', 'HB3'], 'ASNQD2':['HD21', 'HD22'], 'ASPQB':['HB2', 'HB3'], 'CYSQB':['HB2', 'HB3'], 'CYSSQB':['HB2', 'HB3'], 'GLNQB':['HB2', 'HB3'], 'GLNQG':['HG2', 'HG3'], 'GLNQE2':['HE21', 'HE22'], 'GLUQB':['HB2', 'HB3'], 'GLUQG':['HG2', 'HG3'], 'GLYQA':['HA2', 'HA3'], 'HISQB':['HB2', 'HB3'], 'ILEQG1':['HG11', 'HG12', 'HG13'], 'ILEQG2':['HG21', 'HG22', 'HG23'], 'ILEQD1':['HD11', 'HD12', 'HD13'], 'LEUQB':['HB2', 'HB3'], 'LEUQD1':['HD11', 'HD12', 'HD13'], 'LEUQD2':['HD21', 'HD22', 'HD23'], 'LEUQQD':['HD11', 'HD12', 'HD13','HD21', 'HD22', 'HD23'], 'LYSQG':['HG2', 'HG3'], 'LYSQD':['HD2', 'HD3'], 'LYSQE':['HE2', 'HE3'], 'LYSQZ':['HZ2', 'HZ3'], 'METQB':['HB2', 'HB3'], 'METQG':['HG2', 'HG3'], 'METQE':['HE1', 'HE2', 'HE3'], 'PHEQB':['HB2', 'HB3'], 'PHEQD':['HD1', 'HD2'], 'PHEQE':['HE1', 'HE2'],'PROQB':['HB2', 'HB3'], 'PROQG':['HG2', 'HG3'], 'PROQD':['HD2', 'HD3'], 'SERQB':['HB2', 'HB3'], 'THRQG2':['HG21', 'HG22', 'HG23'], 'TRPQB':['HB2', 'HB3'], 'TYRQB':['HB2', 'HB3'], 'TYRQD':['HD1', 'HD2'], 'TYRQE':['HE1', 'HE2'], 'VALQB':['HB2', 'HB3'], 'VALQG1':['HG11', 'HG12', 'HG13'], 'VALQG2':['HG21', 'HG22', 'HG23'], 'VALQQG':['HG11', 'HG12', 'HG13','HG21', 'HG22', 'HG23']}
cwd = '/Volumes/common/Kinases/FGFR3/FGFR3_459-755_C482A_C582S/Structure_Calc/cyana_25/'
# cwd = '/Users/mclay1/FGFR3_structure/cyana_23/'

in_pdb = cwd + 'final.pdb'
fupl = cwd + 'final.upl'
calc = cwd + 'CALC.cya'
manualongcons = [line.strip() for line in open(calc).readlines() if line.strip() and '.upl' in line][0].split()[2].split(',')
upls = [con for con in manualongcons if 'upl' in con and 'hbond' not in con]
fovw = cwd + 'final.ovw'
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
# Fixit = {}
for (start,end) in zip(Starts,Ends):
	n+=1
	print('Reading coordinates for model %d' %n)
	Coor = eval('Coor' + str(n))
	for x in range(start,end,1):
		line = pdb[x]
		if line[0:4] == "ATOM" or line[0:4] == 'HETA':
			index = '%3s %4s %-4s'%(line[17:20].strip(),line[22:26].strip(),line[12:16].strip())
			Coor[index] = [float(line[30:38]),float(line[38:46]),float(line[46:54])]
			## Generate pseudo atoms for methyl groups
			# print(index)
			if line[17:20].strip() in ['ILE','LEU','ALA','VAL','MET','THR'] and line[12:16].strip() in ['HB1','HE1','HD11','HD21','HG11','HG21']:
				index = '%3s %4s %-4s'%(line[17:20].strip(),line[22:26].strip(),line[12:16].strip()[:-1].replace('H','Q'))
				# print('methyl ' + index)
				line2 = pdb[x+1]
				line3 = pdb[x+2]
				# Fixit[line[17:20].strip()+line[12:16].strip()[:-1].replace('H','Q')] = [line[12:16].strip(),line2[12:16].strip(),line3[12:16].strip()]
				xcoor = (float(line[30:38]) + float(line2[30:38]) + float(line3[30:38]))/3.0
				ycoor = (float(line[38:46]) + float(line2[38:46]) + float(line3[38:46]))/3.0
				zcoor = (float(line[46:54]) + float(line2[46:54]) + float(line3[46:54]))/3.0
				Coor[index] = [xcoor,ycoor,zcoor]
			## Generate pseudo atoms for side chain exchangeable atoms 
			if line[17:20].strip() not in ['ILE','LEU','ALA','VAL','MET','THR'] and line[12:16].strip()[0] == 'H':
				if line[17:20].strip() + 'Q' + line[12:16].strip()[1:-1] in PseudoAtoms.keys():
					if line[12:16].strip() == PseudoAtoms[line[17:20].strip() + 'Q' + line[12:16].strip()[1:-1]][0]:
						index = '%3s %4s %-4s'%(line[17:20].strip(),line[22:26].strip(),'Q' + line[12:16].strip()[1:-1])
						print('side chain '+ index)
						line2 = pdb[x+1]
						# Fixit[line[17:20].strip()+'Q' + line[12:16].strip()[1:-1]] = [line[12:16].strip(),line2[12:16].strip()]
						xcoor = (float(line[30:38]) + float(line2[30:38]))/2.0
						ycoor = (float(line[38:46]) + float(line2[38:46]))/2.0
						zcoor = (float(line[46:54]) + float(line2[46:54]))/2.0
						Coor[index] = [xcoor,ycoor,zcoor]

def getdistance(doner, acceptor,PDBdict):
	resn1, resi1, name1 = doner.split()
	resn2, resi2, name2 = acceptor.split()
	distances = []
	if 'Q' in name1 and 'Q' not in name2:
		atoms1list = ['%3s %4s %-4s'%(resn1.replace('HIST','HIS'),resi1,atom) for atom in PseudoAtoms[resn1.replace('HIST','HIS')+name1]]
		atoms1list.append(doner)
		for a1 in atoms1list:
			(x1,y1,z1) = PDBdict[a1]
			(x2,y2,z2) = PDBdict[acceptor]
			distances.append(np.round(np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2),2))
			print('%s %s d=%1.2f' %(doner,acceptor, np.round(np.sqrt((x1-x2)**2 + (y1-y2)**2+ (z1-z2)**2),2)))
		d = min(distances)
	if 'Q' not in name1 and 'Q' in name2:
		atoms2list = ['%3s %4s %-4s'%(resn2.replace('HIST','HIS'),resi2,atom) for atom in PseudoAtoms[resn2.replace('HIST','HIS')+name2]]
		atoms2list.append(acceptor)
		for a2 in atoms2list:
			(x1,y1,z1) = PDBdict[doner]
			(x2,y2,z2) = PDBdict[a2]
			distances.append(np.round(np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2),2))
			print('%s %s d=%1.2f' %(doner, a2,np.round(np.sqrt((x1-x2)**2 + (y1-y2)**2+ (z1-z2)**2),2)))
		d = min(distances)
	if 'Q' in name1 and 'Q' in name2:
		atoms1list = ['%3s %4s %-4s'%(resn1.replace('HIST','HIS'),resi1,atom) for atom in PseudoAtoms[resn1.replace('HIST','HIS')+name1]]
		atoms2list = ['%3s %4s %-4s'%(resn2.replace('HIST','HIS'),resi2,atom) for atom in PseudoAtoms[resn2.replace('HIST','HIS')+name2]]
		atoms1list.append(doner)
		atoms2list.append(acceptor)
		for a1 in atoms1list:
			for a2 in atoms2list:
				(x1,y1,z1) = PDBdict[a1]
				(x2,y2,z2) = PDBdict[a2]
				distances.append(np.round(np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2),2))
				print('%s %s d=%1.2f' %(a1,a2,np.round(np.sqrt((x1-x2)**2 + (y1-y2)**2+ (z1-z2)**2),2)))
		d = min(distances)
	if 'Q' not in name1 and 'Q' not in name2:
		(x1,y1,z1) = PDBdict[atom1]
		(x2,y2,z2) = PDBdict[atom2]
		d = np.round(np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2),2)
	print('used %s %s d=%1.2f' %(doner, acceptor,d))
	return d



pd.set_option("display.precision", 2)
DistancesDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv','upl','plist','peak'])
for line in open(fupl).readlines():
	if 'SUP' not in line:
		pass 
	else:
		cns = line.strip().split()
		atom1 = '%3s %4s %-4s'%(cns[1].replace('HIST','HIS'),cns[0],cns[2])
		atom2 = '%3s %4s %-4s'%(cns[4].replace('HIST','HIS'),cns[3],cns[5])
		for mnum in range(1,21,1):
			Coor = eval('Coor' + str(mnum))
			d = getdistance(atom1, atom2,Coor)
			DistancesDF.loc[atom1 + ' ' + atom2,mnum] = d

DistancesDF['mean'] = np.round(DistancesDF.mean(axis=1),2)
DistancesDF['stdv'] = np.round(DistancesDF.std(axis=1),2)
for line in open(fupl).readlines():
	if 'SUP' not in line:
		pass 
	else:
		cns = line.strip().split()
		atom1 = '%3s %4s %-4s'%(cns[1].replace('HIST','HIS'),cns[0],cns[2])
		atom2 = '%3s %4s %-4s'%(cns[4].replace('HIST','HIS'),cns[3],cns[5])
	DistancesDF.loc[atom1 + ' ' + atom2,'upl'] = float(cns[6])
	DistancesDF.loc[atom1 + ' ' + atom2,'peak'] = cns[8]
	DistancesDF.loc[atom1 + ' ' + atom2,'plist'] = cns[10]
	ovwent = '%-4s %4s %4s - %-4s %4s %4s' %(cns[2], cns[1],cns[0],cns[5],cns[4],cns[3])
	# if 'peak %s list %s' %(cns[8], cns[10]) in open(fovw).read():
	if ovwent in open(fovw).read():
		DistancesDF.loc[atom1 + ' ' + atom2,'violation'] = 'yes'

DistancesDF['diff'] = DistancesDF['mean'] - DistancesDF['upl']
DistancesDF.to_csv('final_distances.csv')
print(DistancesDF)
print(DistancesDF.shape)

for upl in upls:
	print(upl)
	DistancesDF = pd.DataFrame(columns=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,'mean','stdv','upl'])
	for line in open(cwd+upl).readlines():
		if '#'in line[0:4]:
			pass 
		else:
			cns = line.strip().split()
			atom1 = '%3s %4s %-4s'%(cns[1].replace('HIST','HIS'),cns[0],cns[2])
			atom2 = '%3s %4s %-4s'%(cns[4].replace('HIST','HIS'),cns[3],cns[5])
			for mnum in range(1,21,1):
				Coor = eval('Coor' + str(mnum))
				d = getdistance(atom1, atom2,Coor)
				DistancesDF.loc[atom1 + ' ' + atom2,mnum] = d
	DistancesDF['mean'] = np.round(DistancesDF.mean(axis=1),2)
	DistancesDF['stdv'] = np.round(DistancesDF.std(axis=1),2)
	for line in open(cwd+upl).readlines():
		if '#'in line[0:4]:
			pass 
		else:
			cns = line.strip().split()
			atom1 = '%3s %4s %-4s'%(cns[1].replace('HIST','HIS'),cns[0],cns[2])
			atom2 = '%3s %4s %-4s'%(cns[4].replace('HIST','HIS'),cns[3],cns[5])
		DistancesDF.loc[atom1 + ' ' + atom2,'upl'] = float(cns[6])
		ovwent = '%-4s %4s %4s - %-4s %4s %4s' %(cns[2], cns[1],cns[0],cns[5],cns[4],cns[3])
		print(ovwent)
		if ovwent in open(fovw).read():
			DistancesDF.loc[atom1 + ' ' + atom2,'violation'] = 'yes'
	DistancesDF['diff'] = DistancesDF['mean'] - DistancesDF['upl']
	DistancesDF.to_csv(upl.replace('.upl','_distances.csv'))
	print(DistancesDF)
	print(DistancesDF.shape)




