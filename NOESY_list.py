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
from scipy.stats import trim_mean
from scipy.stats import tstd
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker 
import seaborn as sb
import re
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
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
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
mpl.mathtext.FontConstantsBase.sup1 = 0.25

# cwd = '/Volumes/common/Kinases/FGFR3/FGFR3_459-755_C482A_C582S/Structure_Calc/cyana_25/'
cwd = '/Volumes/common/Kinases/FGFR2/FGFR2_467-768_C491A/Structure_Calc/cyana_21/'
calc = cwd + 'CALC.cya'
prots = [line.strip() for line in open(calc).readlines() if line.strip() and 'prot' in line][0].split()[2].split(',')
print(prots)
init = cwd + 'init.cya'
print(open(init).readlines()[0].strip().split(':=')[-1])
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
			if assign not in Assignments:
				Assignments.append(assign)

cya_plists = [line.strip() for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
noa = cwd + 'cycle7.noa'

cya_plists = [line.strip().replace('.peaks','') for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')

protdict = {}
for x in range(len(cya_plists)):
	exec("plist{:}_unused_w1 = []".format(str(x+1)))
	exec("plist{:}_unused_w2 = []".format(str(x+1)))
	exec("plist{:}_unused = []".format(str(x+1)))
	if x < len(prots):
		protlist = eval("{:}_list".format(prots[x].replace('.prot','')))
	if x > len(prots):
		protlist = eval("{:}_list".format(prots[-1].replace('.prot','')))
	protdict[x+1] = [eval("plist{:}_unused_w1".format(str(x+1))), eval("plist{:}_unused_w2".format(str(x+1))),eval("plist{:}_unused".format(str(x+1))),protlist]

noalines = open(noa).readlines()
pl = 0
unusedupls = []
unusedupl = open('final_unused_3.upl','w')
usedupl = open('final_used.upl','w')
Used_upls = []
Used_aupls = []
for noelist in cya_plists:
	pl+=1
	nota,single,amb,notused,peak,incr,dia,ndia  = 0, 0, 0, 0, 0, 0, 0,0
	for x in range(len(noalines)):
		line = noalines[x]
		if noelist in line and 'out of' in noalines[x+1]:
			peak+=1
			peakn = line.strip().split()[1]
			dist = line.strip().split()[-2]
			if '0 out of 0' in noalines[x+1]:
				nota+=1
			if '0 out of' in noalines[x+1] and '0 out of 0' not in noalines[x+1]:
				notused+=1
			if '1 out of' in noalines[x+1] and 'diagonal' in line:
				single+= 1
			if '0 out of' not in noalines[x+1] and 'diagonal' not in line:
				single+= 1
				# print(noalines[x+ int(noalines[x+1].split()[3])+2])
				for y in range(2,int(noalines[x+1].split()[3])+2,1):
					cns = noalines[x+y].strip().split()
					if noalines[x+y].strip()[0] in ['!','*'] and cns[4] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, maxd = cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, float(cns[13].split('-')[0])
					if noalines[x+y].strip()[0] not in ['!','*'] and cns[3] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2,pshift, maxd = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, float(cns[12].split('-')[0])
					upl = '%4d %-4s %-4s %4d %-4s %-4s' %(resi1, resn1, atom1, resi2, resn2, atom2)
					upl2 = '%4d %-4s %-4s %4d %-4s %-4s' %(resi2, resn2, atom2, resi1, resn1, atom1)
					if upl not in Used_upls and pshift > 0.9 and maxd < 8.0:
						# print(upl)
						Used_aupls.extend([upl,upl2])
						usedupl.write('%4d %-4s %-4s %4d %-4s %-4s    %s  #peak %s #plist %d #Pshift %0.2f\n' %(resi1, resn1, atom1, resi2, resn2, atom2, dist, peakn, pl, pshift ))
			if noalines[x+1].strip().split()[0] > '1' and 'diagonal' not in line:
				amb+=1
			if 'increased' in line:
				incr+=1 
			if 'diagonal' in line and '0 out of' not in noalines[x+1]:
				dia+=1
			if 'diagonal' in line and '0 out of' in noalines[x+1]:
				ndia+=1
			if '0 out of' in noalines[x+1] and '0 out of 0' not in noalines[x+1]:
				for y in range(2,int(noalines[x+1].split()[3])+2,1):
					cns = noalines[x+y].strip().split()
					if noalines[x+y].strip()[0] in ['!','*']:
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, maxd = cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, float(cns[13].split('-')[0])
					if noalines[x+y].strip()[0] not in ['!','*']:
						atom1,resn1, resi1, atom2, resn2, resi2,pshift, maxd = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, float(cns[12].split('-')[0])
					upl = '%4d %-4s %-4s %4d %-4s %-4s' %(resi1, resn1, atom1, resi2, resn2, atom2)
					upl2 = '%4d %-4s %-4s %4d %-4s %-4s' %(resi2, resn2, atom2, resi1, resn1, atom1)
					if upl not in Used_upls and pshift > 0.9 and maxd < 9.0:
						# print(upl)
						Used_upls.extend(upl)
						unusedupl.write('{:<4d} {:^4s} {:^4s} {:<4d} {:^4s} {:^4s}    {:}  #peak {:} #plist {:} #Pshift {:0.2f} out of {:}\n'.format(resi1, resn1, atom1, resi2, resn2, atom2, dist, peakn, pl, pshift,noalines[x+1].split()[3]))
						unusedupls.append('{:<4d} {:^4s} {:^4s} {:<4d} {:^4s} {:^4s}   #plist {:} #options {:}'.format(resi1, resn1, atom1, resi2, resn2, atom2, pl,noalines[x+1].split()[3]))
						w1unplist, w2unplist, unplist, protlist =protdict[pl]
						unplist.append('{:}{:d}-{:} {:}{:d}-{:}   #plist {:} #options {:}'.format(AAA_dict[resn1],resi1 , atom1, AAA_dict[resn2],resi2, atom2, pl,noalines[x+1].split()[3]))
						w1unplist.append('{:}{:d}-{:}'.format(AAA_dict[resn1],resi1 , atom1))
						w2unplist.append('{:}{:d}-{:}'.format(AAA_dict[resn2],resi2, atom2))


pdf = PdfPages('plist_test_plots.pdf')
for x in range(len(cya_plists)):
	w1list,w2list, w1w2list,prot = protdict[x+1]
	df = pd.DataFrame(index=prot)
	print(cya_plists[x])
	print(len(w1list))
	for atom in prot:
		if w1list.count(atom) != 0: df.loc[atom,'atom1'] = w1list.count(atom)
		if w2list.count(atom) != 0: df.loc[atom,'atom2'] = w2list.count(atom)
	df['total'] =df['atom1']+df['atom2']
	df1 = df[(df['atom1'] > 1.0) ].copy(deep=True)
	if len(df1.index.tolist()) > 0:
		nsubplots = round(len(df1.index.tolist())/30,0)
		if round(len(df1.index.tolist())/30,1) - nsubplots > 0.0:
			nsubplots = nsubplots + 1
		if nsubplots == 0 : nsubplots = 1
		print(nsubplots)
		fig_height = 3.0 * nsubplots
		entry_width = 5.0/30
		fig_width = 0.78 + entry_width * 30
		if fig_width < 3.0: 
			fig_width = 3.0
		if fig_height <= 2.0: 
			fig_height = 3.0
		fig=plt.figure(figsize=(fig_width,fig_height))
		spi = 0
		for i in range(0,len(df1.index.tolist()),30):
			spi = spi + 1
			z = i
			temp = []
			for y in range(30):
				temp.append(df1.index.tolist()[z])
				z = z +1 
				if z == len(df1.index.tolist()):break
			dfp = df1.reindex(temp)
			ax = fig.add_subplot(int(nsubplots),1,spi)
			ax.bar(dfp.index.tolist(), dfp['atom1'],0.9, color = 'blue', edgecolor='none', label = 'atom1')
			ax.tick_params(axis='x', labelrotation = 90)
		ax.set_title(cya_plists[x] + ' atom1')
		ax.set_xlabel('Residue Number')
		plt.tight_layout(pad = 0.4, w_pad = 0.4, h_pad = 0.4)
		pdf.savefig(transparent=True)
		plt.close()
	df2 = df[(df['atom2'] > 1.0) ].copy(deep=True)
	if len(df2.index.tolist()) > 0:
		nsubplots = round(len(df2.index.tolist())/30,0)
		if round(len(df2.index.tolist())/30,1) - nsubplots > 0.0:
			nsubplots = nsubplots + 1
		if nsubplots == 0: nsubplots = 1
		fig_height = 3.0 * nsubplots
		if fig_height <= 2.0: 
			fig_height = 3.0
		entry_width = 5.0/30
		fig_width = 0.78 + entry_width * 30
		if fig_width < 3.0: 
			fig_width = 3.0
		fig=plt.figure(figsize=(fig_width,fig_height))
		spi = 0
		for i in range(0,len(df2.index.tolist()),30):
			spi = spi + 1
			temp = []
			z = i
			for y in range(30):
				temp.append(df2.index.tolist()[z])
				z = z + 1 
				if z == len(df2.index.tolist()):break
			dfp = df2.reindex(temp)
			ax = fig.add_subplot(int(nsubplots),1,spi)
			ax.bar(dfp.index.tolist(), dfp['atom2'],0.9,color = 'purple', edgecolor='none', label = 'atom2')
			ax.tick_params(axis='x', labelrotation = 90)
		ax.set_title(cya_plists[x] + ' atom2')
		ax.set_xlabel('Residue Number')
		plt.tight_layout(pad = 0.4, w_pad = 0.4, h_pad = 0.4)
		pdf.savefig(transparent=True)
		plt.close()
	df3 = df[(df['total'] > 1.0) ].copy(deep=True)
	if len(df3.index.tolist()) > 0:
		nsubplots = round(len(df3.index.tolist())/30,0)
		if round(len(df3.index.tolist())/30,1) - nsubplots > 0.0:
			nsubplots = nsubplots + 1
		if nsubplots == 0: nsubplots = 1
		fig_height = 3.0 * nsubplots
		if fig_height <= 2.0: 
			fig_height = 3.0
		entry_width = 5.0/30
		fig_width = 0.78 + entry_width * 30
		if fig_width < 3.0: 
			fig_width = 3.0
		fig=plt.figure(figsize=(fig_width,fig_height))
		spi = 0
		for i in range(0,len(df3.index.tolist()),30):
			spi = spi + 1
			temp = []
			z = i
			for y in range(30):
				temp.append(df3.index.tolist()[z])
				z = z + 1 
				if z == len(df3.index.tolist()):break
			dfp = df3.reindex(temp)
			ax = fig.add_subplot(int(nsubplots),1,spi)
			ax.bar(dfp.index.tolist(), dfp['atom2'],0.9,color='orange', edgecolor='none', label = 'atom2')
			ax.tick_params(axis='x', labelrotation = 90)
		ax.set_title(cya_plists[x] + ' total')
		ax.set_xlabel('Residue Number')
		plt.tight_layout(pad = 0.4, w_pad = 0.4, h_pad = 0.4)
		pdf.savefig(transparent=True)
	plt.close()
pdf.close()

	# print(df2)
	# print(df.dropna(inplace=True).shape)
# sx = range(0,len(Isolated2),4)[-1]
# PDBiso = '   '
# for x in range(0, len(Isolated2),4)[:-1]:
# 	self.summary_list.append('   %s %s %s %s' %(Isolated2[x],Isolated2[x+1],Isolated2[x+2],Isolated2[x+3]))
# lx = len(Isolated2) - sx

	# print(noelist)
	# print("Total number of peaks %d" %peak)
	# print('Number with single assignment %d'%single)
	# print('Number ambigious assignment %d'%amb)
	# print('Number assigned but unused %d'%notused)
	# print('Number never assinged %d'%nota)
	# print('Number of diagonal assignments %d' %dia)
	# print('Number of diagonal not assignmened %d' %ndia)
	# print('Number of peaks with increase distacne cut off %d' %incr)
	# print('total %d' %(nota+notused+single+amb+dia+ndia))
	# print()
unusedupl.close()
# print(plist1_unused)
# print(len(plist1_unused))


