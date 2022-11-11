'''
# ------------------------------------------------------------------------------
#
# Created by : Mary Clay PhD
# e-mail: mary.clay@stjude.org
# St Jude Children's Research Hospital 
# Department of Structural Biology Memphis, TN 
# 06/28/2022
# ------------------------------------------------------------------------------

'''
import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from scipy.stats import trim_mean
from scipy.stats import tstd
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker 
import seaborn as sb
import re
import mplcursors
####----------------------------------------------------------------------------------------------####
##																									##
##			 	Setting controlling plot appearance: Font, line widths, ticks, and legend  			##
##																									##
####----------------------------------------------------------------------------------------------####
AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H","HIST": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S'}


# cwd = '/Volumes/common/Kinases/FGFR3/FGFR3_459-755_C482A_C582S/Structure_Calc/cyana_25/'
cwd = '/Volumes/common/Kinases/FGFR2/FGFR2_467-768_C491A/Structure_Calc/cyana_46/'
calc = cwd + 'CALC.cya'
cya_plists = [line.strip() for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
prots = [line.strip() for line in open(calc).readlines() if line.strip() and 'prot' in line][0].split()[2].split(',')
init = cwd + 'init.cya'
# print(open(init).readlines()[0].strip().split(':=')[-1])
seq = [line.strip().split() for line in open(cwd + open(init).readlines()[0].strip().split(':=')[-1] + '.seq').readlines() if '#' != line[0]]
Seqdict = {}
Sequence = []
for resn,resi in seq:
	Seqdict[resi] = AAA_dict[resn] + resi
	Sequence.append(AAA_dict[resn] + resi)
Assignments,exceptions = [], []
for prot in prots:
	for line in open(cwd + prot.replace('.prot','-final.prot')).readlines():
		if line.strip():
			if '#' != line.strip()[0]:
				resi = line.strip().split()[4]
				atom = line.strip().split()[3]
				assign = Seqdict[resi] + '-' + atom
				if assign not in Assignments and atom[0] in ['H','Q']:
					Assignments.append(assign)
				# if assign[0] in ['T','G','M'] and atom == 'H':
				# 	exceptions.append(assign)
				if assign[0] in ['F','Y'] and atom in ['HB2','HB3']:
					exceptions.append(assign)
print(exceptions)
assigndict = {}
from itertools import combinations
ADpairs = [ '{:} - {:}'.format(comb[0],comb[1]) for comb in combinations(Assignments,2)]
for comb in ADpairs:
	assigndict[comb] = []

noa7 = cwd + 'cycle7.noa'
noalines = open(noa7).readlines()
noassignment = open('no_assignment.noa','w')
unused_assignment = open('unused_assignment.noa','w')
Used_aupls = []
for x in range(len(noalines)):
	line = noalines[x]
	if 'Peak' in noalines[x] and 'out of' in noalines[x+1]:
		if '0 out of' not in noalines[x+1] and 'diagonal' not in noalines[x]:
			peak = int(noalines[x].strip().split()[1])
			plist = noalines[x].strip().split()[3]
			for y in range(2,int(noalines[x+1].split()[0])+2,1):
				cns = noalines[x+y].strip().split()
				if cns[4] == '+':
					atom1,resn1, resi1, atom2, resn2, resi2, pshift,drange= cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, cns[13]
				if cns[3] == '+':
					atom1,resn1, resi1, atom2, resn2, resi2,pshift,drange = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, cns[12]
				group1 = '{:}{:}-{:}'.format(AAA_dict[resn1],resi1, atom1)
				group2 = '{:}{:}-{:}'.format(AAA_dict[resn2],resi2, atom2)
				outline = '    {:>10} - {:<10}   {:<9}  Peak {:} from {:}  pshift {:0.2f}\n'.format(group1, group2,drange,peak,plist,pshift)
				if '{:} - {:}'.format(group1,group2)in ADpairs:
					assigndict['{:} - {:}'.format(group1,group2)].append(outline)
				if '{:} - {:}'.format(group2,group1)in ADpairs:
					assigndict['{:} - {:}'.format(group2,group1)].append(outline)
print('finised good')
ADpairs2 = []
for con in ADpairs:
	if len(assigndict[con]) > 1:
		ADpairs2.append(con)
for x in range(len(noalines)):
	line = noalines[x]
	if 'Peak' in noalines[x] and 'out of' in noalines[x+1]:
		if '0 out of' in noalines[x+1] and '0 out of 0' not in noalines[x+1]:
			peak = int(noalines[x].strip().split()[1])
			plist = noalines[x].strip().split()[3]
			unused_assignment.write(noalines[x])
			for y in range(2,int(noalines[x+1].split()[3])+2,1):
				cns = noalines[x+y].strip().split()
				if noalines[x+y].strip()[0] in ['!','*']:
					atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, cns[13]
				if noalines[x+y].strip()[0] not in ['!','*']:
					atom1,resn1, resi1, atom2, resn2, resi2,pshift, drange = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, cns[12]
				group1 = '{:}{:}-{:}'.format(AAA_dict[resn1],resi1, atom1)
				group2 = '{:}{:}-{:}'.format(AAA_dict[resn2],resi2, atom2)
				outline = '#    {:>10} - {:<10}   {:<9}  Peak {:} from {:}  pshift {:0.2f} unused\n'.format(group1, group2,drange,peak,plist,pshift)
				if '{:} - {:}'.format(group1,group2)in ADpairs2:
					assigndict['{:} - {:}'.format(group1,group2)].append(outline)
				if '{:} - {:}'.format(group2,group1)in ADpairs2:
					assigndict['{:} - {:}'.format(group2,group1)].append(outline)
				if pshift > 0.75 :
					unused_assignment.write(noalines[x+y])
			if 'Violated' in noalines[x+y+1]:
				unused_assignment.write(noalines[x+y+1])
			unused_assignment.write('\n')
		if '0 out of 0' in noalines[x+1]:
			noassignment.write(noalines[x])
noassignment.close()
unused_assignment.close()
assigned = open('assinged.list','w')
questionable = open('Q_assinged.list','w')
qpeaks = []
for con in ADpairs:
	if len(assigndict[con]) > 1:
		assigned.write('{:}  {:}:\n'.format(con,len(assigndict[con])))
		assigned.writelines(assigndict[con])
		assigned.write('\n')
	if len(assigndict[con]) == 1:
		print(con.split('-'))
		if (con.split(' - ')[0] in exceptions) and (con.split(' - ')[1] in exceptions):
			assigned.write('{:}  {:}:\n'.format(con,len(assigndict[con])))
			assigned.writelines(assigndict[con])
			assigned.write('\n')
		if (con.split(' - ')[0] not in exceptions) and (con.split(' - ')[1] not in exceptions):
			qpeaks.append(assigndict[con][0])
qpeaks = sorted(qpeaks,key = lambda x: (x.split()[7],x.split()[5]))
assigned.close()
questionable.writelines(qpeaks)
questionable.close()

