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

linewidths = 1.5
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.linewidth'] = linewidths
mpl.rcParams['xtick.direction'] = mpl.rcParams['ytick.direction']='out'
mpl.rcParams['xtick.labelsize'] = mpl.rcParams['ytick.labelsize']=10
mpl.rcParams['xtick.major.size'] = mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = mpl.rcParams['ytick.major.width']=linewidths
mpl.rcParams['xtick.minor.size'] = 2
mpl.rcParams['xtick.minor.width'] = linewidths
# mpl.rcParams['axes.spines.right'] = False
# mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['legend.loc'] = 'best'
mpl.rcParams['legend.borderpad'] = 0.01
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['legend.handlelength'] = 0
mpl.rcParams['legend.scatterpoints'] = 1
mpl.rcParams['xtick.major.bottom'] = mpl.rcParams['ytick.major.left'] = True
mpl.rcParams['xtick.major.top'] = mpl.rcParams['ytick.major.right'] = True
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.sf'] = 'sans\\-serif'
plt.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['xtick.major.pad']=mpl.rcParams['ytick.major.pad']= 2
mpl.rcParams['axes.labelpad'] = 3
# mpl.mathtext.FontConstantsBase.sup1 = 0.25

# cwd = '/Volumes/common/Kinases/FGFR3/FGFR3_459-755_C482A_C582S/Structure_Calc/cyana_25/'
cwd = '/Volumes/common/Kinases/FGFR2/FGFR2_467-768_C491A/Structure_Calc/cyana_46/'
calc = cwd + 'CALC.cya'
prots = [line.strip() for line in open(calc).readlines() if line.strip() and 'prot' in line][0].split()[2].split(',')
init = cwd + 'init.cya'
# print(open(init).readlines()[0].strip().split(':=')[-1])
seq = [line.strip().split() for line in open(cwd + open(init).readlines()[0].strip().split(':=')[-1] + '.seq').readlines() if '#' != line[0]]
Seqdict = {}
Sequence = []
for resn,resi in seq:
	Seqdict[resi] = AAA_dict[resn] + resi
	Sequence.append(AAA_dict[resn] + resi)
Assignments = []
for prot in prots:
	exec("{:}_list = []".format(prot.replace('.prot','')))
	protlist = eval("{:}_list".format(prot.replace('.prot','')))
	Assignlist = eval("{:}_list".format(prot.replace('.prot','')))
	for line in open(cwd + prot.replace('.prot','-final.prot')).readlines():
		if '#' != line[0] and line.strip():
			resi = line.strip().split()[4]
			atom = line.strip().split()[3]
			assign = Seqdict[resi] + '-' + atom
			protlist.append(assign)
			if assign not in Assignments and atom[0] in ['H','Q']:
				Assignments.append(assign)
assigndict = {}
from itertools import combinations
ADpairs = [ '{:}-{:}'.format(comb[0],comb[1]) for comb in combinations(Assignments,2)]
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
			peak = noalines[x].strip().split()[1]
			plist = noalines[x].strip().split()[3]
			for y in range(2,int(noalines[x+1].split()[0])+2,1):
				cns = noalines[x+y].strip().split()
				if cns[4] == '+':
					atom1,resn1, resi1, atom2, resn2, resi2, pshift,drange= cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, cns[13]
				if cns[3] == '+':
					atom1,resn1, resi1, atom2, resn2, resi2,pshift,drange = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, cns[12]
				group1 = '{:}{:}-{:}'.format(AAA_dict[resn1],resi1, atom1)
				group2 = '{:}{:}-{:}'.format(AAA_dict[resn2],resi2, atom2)
				outline = '    {:<10} - {:<10}   {:<9}  #peak {:} from {:}  pshift {:0.2f}\n'.format(group1, group2,drange,peak,plist,pshift)
				if '{:}-{:}'.format(group1,group2)in ADpairs:
					assigndict['{:}-{:}'.format(group1,group2)].append(outline)
				if '{:}-{:}'.format(group2,group1)in ADpairs:
					assigndict['{:}-{:}'.format(group2,group1)].append(outline)
print('finised good')
ADpairs2 = []
for con in ADpairs:
	if len(assigndict[con]) > 1:
		ADpairs2.append(con)
for x in range(len(noalines)):
	line = noalines[x]
	if 'Peak' in noalines[x] and 'out of' in noalines[x+1]:
		if '0 out of' in noalines[x+1] and '0 out of 0' not in noalines[x+1]:
			peak = noalines[x].strip().split()[1]
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
				outline = '#    {:<10} - {:<10}   {:<9}  #peak {:} from {:}  pshift {:0.2f} unused\n'.format(group1, group2,drange,peak,plist,pshift)
				if '{:}-{:}'.format(group1,group2)in ADpairs2:
					assigndict['{:}-{:}'.format(group1,group2)].append(outline)
				if '{:}-{:}'.format(group2,group1)in ADpairs2:
					assigndict['{:}-{:}'.format(group2,group1)].append(outline)
				if pshift > 0.75 :
					unused_assignment.write(noalines[x+y])
			if 'Violated' in noalines[x+y+1]:
				unused_assignment.write(noalines[x+y+1])
			unused_assignment.write('\n')
		if '0 out of 0' in noalines[x+1]:
			noassignment.write(noalines[x])
			noassignment.write(noalines[x+1])
			noassignment.write('\n')
noassignment.close()
unused_assignment.close()
assigned = open('assinged.list','w')
questionable = open('Q_assinged.list','w')
for con in ADpairs:
	if len(assigndict[con]) > 1:
		assigned.write('{:}  {:}:\n'.format(con,len(assigndict[con])))
		assigned.writelines(assigndict[con])
		assigned.write('\n')
	if len(assigndict[con]) == 1:
		questionable.write('{:}  {:}:\n'.format(con,len(assigndict[con])))
		questionable.writelines(assigndict[con])
		questionable.write('\n')
assigned.close()
questionable.close()

