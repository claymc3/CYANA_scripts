import os
import sys
import re
import matplotlib as mpl
import numpy as np

AAA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
	'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
	'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 
	'ADE': 'A', 'GUA': 'G', 'CYT': 'C', 'URA': 'U', 'THY': 'T','HIST':'H', 'HISE':'H'}

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
#'ALAH':'N','CYSH':'N','ASPH':'N','GLUH':'N','PHEH':'N','GLYH':'N','HISH':'N','ILEH':'N','LYSH':'N','LEUH':'N','METH':'N','ASNH':'N','GLNH':'N','ARGH':'N','SERH':'N','THRH':'N','VALH':'N','TRPH':'N','TYRH':'N',

# cwd = os.getcwd() + '/'
# calc = cwd + 'CALC.cya'
cwd = '/Volumes/common/Kinases/FGFR3/FGFR3_459-755_C482A_C582S/Structure_Calc/cyana_23/'
calc = cwd+'CALC.cya'
seqfile = cwd+'FGFR3_AS_KD.seq'
fupl = cwd+'final.upl'

manualongcons = [line.strip() for line in open(calc).readlines() if line.strip() and '.upl' in line][0].split()[2].split(',')
upls = [con for con in manualongcons if 'upl' in con and 'hbond' not in con]

## Read in the peak list files, and generate a dictionary connecting index to peaks list file. 
cya_plists = [line.strip().replace('.peaks','-cycle7.peaks') for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
plists = {}
x = 0
for plist in cya_plists:
	x+=1
	plists[str(x)]=plist

print(plists)
### Read Log file and extract the calibration constants for each peaks list file. They will be connected by index not list file name
Calibration_Cnsts = {}
print(cya_plists)
log = '/Volumes/common/Kinases/FGFR3/FGFR3_459-755_C482A_C582S/Structure_Calc/cyana_23/log'
for line in open(log).readlines():
	if 'Calibration constant for peak list' in line:
		print(line)
		Calibration_Cnsts[line.replace(':','').strip().split()[-2]] = float(line[-9:-1])
print(Calibration_Cnsts)

cya_protlist = [line.strip().replace('.prot','-final.prot') for line in open(calc).readlines() if line.strip() and 'prot' in line][0].split()[2].split(',')
print(cya_protlist)

num2AAA = {} #key = index, value = group
pszLines =[line.rstrip() for line in open(seqfile).readlines() if line.rstrip() and ">" not in line and "#" not in line]
for line in pszLines:
	num2AAA[line.split()[1]] = line.split()[0]

    # cya_plist = [line.strip() for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
    # for x in range(len(cya_plist)):
    #   if cya_plist[x].replace('.peaks','-cycle7.peaks') == fin.split('/')[-1]:
    #     lupl = [line.strip() for line in open(fupl).readlines() if line.strip() and 'plist '+ str(x+1) in line]
    #     for upl in lupl:
    #       d = upl.split()[6]
    #       pkn = upl.split()[8]
    #       dist[pkn] = d 
### Now lest take each peak file, translate the assignment to a upl fromat restrain and calculate the distance based on the calibraton constant
for x in range(len(cya_plists)):
	dist = {}
	Ical = Calibration_Cnsts[str(x+1)]
	peaklist = cya_plists[x]
	print(peaklist)
	lupl = [line.strip() for line in open(fupl).readlines() if line.strip() and 'plist '+ str(x+1) in line]
	if x in range(len(cya_protlist)):
		prot = cya_protlist[x]
	if not x in range(len(cya_protlist)):
		prot = cya_protlist[-1]
	protlines = [line.strip() for line in open(cwd+prot).readlines() if line.strip() and line[0] != "#"]
	cyntrans = {} #key = cyana #code; value = [sparky group-atom, chain]
	for line in protlines:
		cyntrans[line.split()[0]]= '%4s %-4s %-4s' %(line.split()[4],num2AAA[line.split()[4]],line.split()[3])
	print(open(cwd+peaklist).readlines()[0][-2])
	if open(cwd+peaklist).readlines()[0][-2] == '3':
		print(peaklist)
		for line in open(cwd+peaklist).readlines():
			if '#' not in line[0:10] and '#VC' not in line:
				ents = line.strip().split()
				d = np.round((Ical/float(ents[6]))**(1/6),2)
				if ents[10] != ents[12]:
					outupl = '%s   %s   %6.2f  %3.2E #peak %s #plist %s' %(cyntrans[ents[10]],cyntrans[ents[12]],d,float(ents[6]),ents[0],str(x+1))
					print(outupl)
	if open(cwd+peaklist).readlines()[0][-2] == '4':
		print(peaklist)
		for line in open(cwd+peaklist).readlines():
			if '#' not in line[0:10] and '#VC' not in line:
				ents = line.strip().split()
				d = np.round((Ical/float(ents[7]))**(1/6),2)
				if ents[11] != ents[13]:
					outupl = '%s   %s   %6.2f  %3.2E #peak %s #plist %s' %(cyntrans[ents[11]],cyntrans[ents[13]],d,float(ents[7]),ents[0],str(x+1))
					print(outupl)
# for line in open(fupl).readlines():


#%4s %s  %-3s   %4s %s  %-3s
