import os
import sys
import re
import matplotlib as mpl

AAA_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
	'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
	'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 
	'ADE': 'A', 'GUA': 'G', 'CYT': 'C', 'URA': 'U', 'THY': 'T','HIST':'H', 'HISE':'H'}
# cwd = os.getcwd() + '/'
# calc = cwd + 'CALC.cya'

calc = '/Users/mclay1/FGFR3_structure/cyana_23/CALC.cya'
seqfile = '/Users/mclay1/FGFR3_structure/cyana_23/FGFR3_AS_KD.seq'


seq_dict = {} #key = index, value = group
pszLines =[line.rstrip() for line in open(seqfile).readlines() if line.rstrip() and ">" not in line and "#" not in line]
for line in pszLines:
	if line.split()[0] in AAA_dict.keys():
		seq_dict[line.split()[1]] = AAA_dict[line.split()[0]] + str(line.split()[1])

protlines= [line.strip() for line in open(fprot).readlines() if line.strip() and line[0] != "#"]
cyntrans = {} #key = cyana #code; value = [sparky group-atom, chain]
    for line in protlines:
      cyntrans[line.split()[0]]= seq_dict[line.split()[4]] + line.split()[3]


calc = '/Users/mclay1/FGFR3_structure/cyana_23/CALC.cya'
cya_plists = [line.strip().replace('.peaks','-cycle7.peaks') for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
plists = {}
x = 0
for plist in cya_plists:
	x+=1
	plists[str(x)]=plist

print(plists)

Calibration_Cnsts = {}
print(cya_plists)
log = '/Users/mclay1/FGFR3_structure/cyana_23/log'
for line in open(log).readlines():
	if 'Calibration constant for peak list' in line:

		Calibration_Cnsts[line.replace(':','').strip().split()[-2]] = float(line[-9:-1])
print(Calibration_Cnsts)
