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

def analize_noa(cwd, outdir, calc, noa7, Seqdict, violdict, qupldict):
	cya_plists = [line.strip() for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
	prots = [line.strip() for line in open(calc).readlines() if line.strip() and 'prot' in line][0].split()[2].split(',')

	plist_dict = {}
	for x in range(len(cya_plists)):
		plist = cya_plists[x].replace('.peaks','')
		plist_dict[plist] = str(x)
		exec("unused{:} = open('{:}{:}_unused.list','w')".format(str(x), outdir, plist))
		exec("no_assign{:} = open('{:}{:}_no_assign.list','w')".format(str(x), outdir, plist))
		#exec("questionable{:} = open('{:}{:}_questionable.list','w')".format(str(x), outdir, plist))
		exec("questlist{:} = []".format(str(x)))
		exec("peaks{:} = {{}}".format(str(x)))
	for x in range(len(cya_plists)):
		pdict = eval('peaks' + str(x))
		for line in open(cwd + cya_plists[x]):
			if line.strip():
				if line.strip()[0] != '#':
					pdict[int(line.split()[0])] = [float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]
	Assignments= []
	for prot in prots:
		for line in open(cwd + prot.replace('.prot','-final.prot')).readlines():
			if line.strip():
				if '#' != line.strip()[0]:
					resi = line.strip().split()[4]
					atom = line.strip().split()[3]
					assign = Seqdict[resi] + '-' + atom
					if assign not in Assignments and atom[0] in ['H','Q']:
						Assignments.append(assign)
	assigndict = {}
	from itertools import combinations
	ADpairs = [ '{:}-{:}'.format(comb[0],comb[1]) for comb in combinations(Assignments,2)]
	for comb in ADpairs:
		assigndict[comb] = []

	noalines = open(noa7).readlines()
	for x in range(len(noalines)):
		line = noalines[x]
		if 'Peak' in noalines[x] and 'out of' in noalines[x+1]:
			if '0 out of' not in noalines[x+1] and 'diagonal' not in noalines[x]:
				peak = int(noalines[x].strip().split()[1])
				plist = noalines[x].strip().split()[3].replace('.peaks','')
				pdict = eval('peaks' + plist_dict[plist])
				questionable = eval('questlist' + plist_dict[plist])
				QF = noalines[x+1].split()[-1].replace(':','')
				for y in range(2,int(noalines[x+1].split()[0])+2,1):
					cns = noalines[x+y].strip().split()
					if cns[4] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[1],cns[2],int(cns[3]), cns[5], cns[6], int(cns[7]), float(cns[10])/100, cns[13]
					if cns[3] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[0],cns[1],int(cns[2]), cns[4], cns[5], int(cns[6]), float(cns[9])/100, cns[12]
					group1 = '{:}{:}-{:}'.format(AAA_dict[resn1],resi1, atom1)
					group2 = '{:}{:}-{:}'.format(AAA_dict[resn2],resi2, atom2)
					conect = '{:}-{:}'.format(group1,group2)
					outline = '{:^28}   {:<9}  Peak {:} from {:}  pshift {:3.2f}\n'.format(conect,drange,peak,plist,pshift)
					if '{:}-{:}'.format(group1,group2)in ADpairs:
						assigndict['{:}-{:}'.format(group1,group2)].append(outline)
					if '{:}-{:}'.format(group2,group1)in ADpairs:
						assigndict['{:}-{:}'.format(group2,group1)].append(outline)
					if float(QF) <= 0.6:
						questionable.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:>9}A  {:<3.2f}   #poor/low support\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,drange,float(QF),'poor constraint'))
					if conect in violdict.keys():
						print(conect)
						print(violdict[conect])
						print(pdict[peak])
						outline = "{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:>9}A  {:<3.2f}  {:}".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,drange,pshift,violdict[conect])
						questionable.append(outline)
					if conect in qupldict.keys():
						print(conect)
						outline = "{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:>9}A  {:<3.2f}  {:}".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,drange,pshift,qupldict[conect])
						questionable.append(outline)
	print('finished assinged')

	ADpairs2 = []
	for con in ADpairs:
		if len(assigndict[con]) >= 1:
			ADpairs2.append(con)
	for x in range(len(noalines)):
		line = noalines[x]
		if 'Peak' in noalines[x] and '0 out of' in noalines[x+1]:
			peak = int(noalines[x].strip().split()[1])
			plist = noalines[x].strip().split()[3].replace('.peaks','')
			pdict = eval('peaks' + plist_dict[plist])
			noassingmnet = eval('no_assign' + plist_dict[plist])
			unused = eval('unused' + plist_dict[plist])
			if '0 out of 0' not in noalines[x+1]:
				nopt = int(noalines[x+1].split()[3])
				for y in range(2,int(noalines[x+1].split()[3])+2,1):
					cns = noalines[x+y].strip().split()
					if noalines[x+y].strip()[0] in ['!','*']:
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, cns[13]
					if noalines[x+y].strip()[0] not in ['!','*']:
						atom1,resn1, resi1, atom2, resn2, resi2,pshift, drange = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, cns[12]
					group1 = '{:}{:}-{:}'.format(AAA_dict[resn1],resi1, atom1)
					group2 = '{:}{:}-{:}'.format(AAA_dict[resn2],resi2, atom2)
					conect = '{:}-{:}'.format(group1,group2)
					outline = '#{:^28}   {:<9}  Peak {:} from {:}  pshift {:0.2f} unused\n'.format(conect,drange,peak,plist,pshift)
					if '{:}-{:}'.format(group1,group2)in ADpairs2:
						assigndict['{:}-{:}'.format(group1,group2)].append(outline)
					if '{:}-{:}'.format(group2,group1)in ADpairs2:
						assigndict['{:}-{:}'.format(group2,group1)].append(outline)
					if pshift > 0.75:
						unused.write("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:>9}A  {:<3.2f}  #{:}\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,drange,pshift,noalines[x+nopt+2].strip()[:-1].replace('Violated','Viol').replace('structures','')))
					if pshift < 0.75 and int(noalines[x+1].split()[3]) == 1:
						noassingmnet.write("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:>9}A  {:<3.2f}  #{:}\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,drange,pshift,noalines[x+nopt+2].strip()[:-1].replace('Violated','Viol').replace('structures','')))
			if '0 out of 0' in noalines[x+1]:
				pdict = eval('peaks' + plist_dict[plist])
				noassingmnet.write("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:>9}A  {:}\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],'none',noalines[x].strip().split()[-2],'na'))
	assigned = open(outdir + 'Assigned.list','w')
	for con in ADpairs:
		if len(assigndict[con]) > 1:
			assigned.write('{:}  {:}:\n'.format(con,len(assigndict[con])))
			assigned.writelines(assigndict[con])
			assigned.write('\n')
		if len(assigndict[con]) == 1:
			pline = assigndict[con][0].split()
			peak = int(pline[3])
			plist = pline[5]
			pshift = pline[7]
			pdict = eval('peaks' + plist_dict[plist])
			questout =  eval('questlist' + plist_dict[plist])
			questout.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:>9}A  {:}\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],pline[0],pline[1],pline[7]))

	for x in range(len(cya_plists)):
		eval("unused{:}.close()".format(str(x)))
		eval("no_assign{:}.close()".format(str(x)))
		exec("questionable{:} = open('{:}{:}_questionable.list','w')".format(str(x), outdir, cya_plists[x].replace('.peaks','')))
		questlist = eval('questlist' + plist_dict[cya_plists[x].replace('.peaks','')])
		qlist = sorted(questlist, key = lambda x: float(x.strip().split()[0]))
		eval("questionable{:}.writelines(qlist)".format(str(x)))
		eval("questionable{:}.close()".format(str(x)))
	assigned.close()

	print("Finished cycle7.noa analysis")


