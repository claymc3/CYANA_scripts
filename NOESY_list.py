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
cwd = '/Volumes/common/Kinases/FGFR2/FGFR2_467-768_C491A/Structure_Calc/cyana_21/'
calc = cwd + 'CALC.cya'
prots = [line.strip() for line in open(calc).readlines() if line.strip() and 'prot' in line][0].split()[2].split(',')
print(prots)
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
			if assign not in Assignments:
				Assignments.append(assign)

cya_plists = [line.strip() for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')
noa = cwd + 'cycle7.noa'

cya_plists = [line.strip().replace('.peaks','') for line in open(calc).readlines() if line.strip() and 'peaks' in line][0].split()[2].split(',')

protdict = {}
for x in range(len(cya_plists)):
	exec("plist{:}_unused = []".format(str(x+1)))
	exec("plist{:}_used = []".format(str(x+1)))
	if x < len(prots):
		protlist = eval("{:}_list".format(prots[x].replace('.prot','')))
	if x > len(prots):
		protlist = eval("{:}_list".format(prots[-1].replace('.prot','')))
	protdict[x+1] = [eval("plist{:}_unused".format(str(x+1))),protlist]

noalines = open(noa).readlines()
pl = 0
unusedupls = []
unusedupl = open('final_unused_3.upl','w')
usedupl = open('final_used.upl','w')
Used_aupls = []
for noelist in cya_plists:
	pl+=1
	for x in range(len(noalines)):
		line = noalines[x]
		if noelist in line and 'out of' in noalines[x+1]:
			peakn = line.strip().split()[1]
			dist = line.strip().split()[-2]
			if '0 out of' not in noalines[x+1] and 'diagonal' not in line:
				# print(noalines[x+ int(noalines[x+1].split()[3])+2])
				for y in range(2,int(noalines[x+1].split()[0])+2,1):
					cns = noalines[x+y].strip().split()
					if cns[4] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, maxd = cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, float(cns[13].split('-')[0])
					if cns[3] == '+':
						atom1,resn1, resi1, atom2, resn2, resi2,pshift, maxd = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, float(cns[12].split('-')[0])
					#if pshift > 0.9 and maxd < 8.0:
						# print(upl)
					Used_aupls.append('{:}{:d} {:^4s} {:}{:d} {:^4s}   #plist {:} #options {:}'.format(AAA_dict[resn1],resi1, atom1, AAA_dict[resn2],resi2, atom2, pl,noalines[x+1].split()[3]))
					usedupl.write('%4d %-4s %-4s %4d %-4s %-4s    %s  #peak %s #plist %d #Pshift %0.2f\n' %(resi1, resn1, atom1, resi2, resn2, atom2, dist, peakn, pl, pshift ))
			if '0 out of' in noalines[x+1] and '0 out of 0' not in noalines[x+1]:
				for y in range(2,int(noalines[x+1].split()[3])+2,1):
					cns = noalines[x+y].strip().split()
					if noalines[x+y].strip()[0] in ['!','*']:
						atom1,resn1, resi1, atom2, resn2, resi2, pshift, maxd = cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, float(cns[13].split('-')[0])
					if noalines[x+y].strip()[0] not in ['!','*']:
						atom1,resn1, resi1, atom2, resn2, resi2,pshift, maxd = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, float(cns[12].split('-')[0])
					upl = '%4d %-4s %-4s %4d %-4s %-4s' %(resi1, resn1, atom1, resi2, resn2, atom2)
					# upl2 = '%4d %-4s %-4s %4d %-4s %-4s' %(resi2, resn2, atom2, resi1, resn1, atom1)
					if pshift > 0.9 and maxd < 9.0:
						# print(upl)
						unusedupl.write('{:<4d} {:^4s} {:^4s} {:<4d} {:^4s} {:^4s}    {:}  #peak {:} #plist {:} #Pshift {:0.2f} out of {:}\n'.format(resi1, resn1, atom1, resi2, resn2, atom2, dist, peakn, pl, pshift,noalines[x+1].split()[3]))
						unusedupls.append('{:}{:d} {:^4s} {:}{:d} {:^4s}   #plist {:} #options {:}'.format(AAA_dict[resn1],resi1, atom1, AAA_dict[resn2],resi2, atom2, pl,noalines[x+1].split()[3]))
						unplist, protlist =protdict[pl]
						unplist.append('{:}{:d} {:} {:}{:d} {:}   #plist {:} #options {:}'.format(AAA_dict[resn1],resi1 , atom1, AAA_dict[resn2],resi2, atom2, pl,noalines[x+1].split()[3]))

# pdf = PdfPages('plist_test_plots.pdf')
for x in range(len(cya_plists)):
	print(cya_plists[x])
	upllist,prot = protdict[x+1]
	df = pd.DataFrame(index=Sequence, columns = Sequence)
	for upl in upllist:
		tvalue = df.loc[upl.split()[0],upl.split()[2]]
		if pd.isna(tvalue): tvalue = 0.0
		df.loc[upl.split()[0],upl.split()[2]] = tvalue + 1
	df = df[df > 1]
	df2 = df.dropna(axis=0, how= 'all')
	df3 = df2.dropna(axis=1, how= 'all')
	fig=plt.figure(figsize=(5.3,5))
	ax = fig.add_subplot(111)
	cmap = sb.color_palette("inferno_r", as_cmap=True)
	cmap.set_under(color='white')
	im = ax.imshow(df3.fillna(0.0), cmap=cmap, vmin = 1.0)
	ax.set_xticks(np.arange(len(df3.columns.to_list())), labels=df3.columns.to_list())
	ax.set_yticks(np.arange(len(df3.index.to_list())), labels=df3.index.to_list())
	ax.tick_params(axis='x', labelrotation = 90)
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right", size="5%", pad=0.05)
	cbar = ax.figure.colorbar(im, cax=cax, shrink=0.5,pad=0.01)
	cbar.set_label(label='Number of Observations',labelpad=1)
	# xlabels = df3.columns.to_list()
	# ylabels = df3.index.to_list()
	cursor = mplcursors.cursor(hover=True).connect("add", lambda sel: sel.annotation.set_text(df3.index.to_list()[sel.index[0]]+'-'+df3.columns.to_list()[sel.index[1]]))	
	ax.set_ylabel('Atom 1')
	ax.set_xlabel('Atom 2')
	ax.set_title(cya_plists[x] + ' Unused Peaks')
	plt.draw()



df = pd.DataFrame(index=Sequence, columns = Sequence)
for upl in unusedupls:
	tvalue = df.loc[upl.split()[0],upl.split()[2]]
	if pd.isna(tvalue): tvalue = 0.0
	df.loc[upl.split()[0],upl.split()[2]] = tvalue + 1
df = df[df > 2]
df2 = df.dropna(axis=0, how= 'all')
df3 = df2.dropna(axis=1, how= 'all')
fig1=plt.figure()
ax = fig1.add_subplot(111)
cmap = sb.color_palette("inferno_r", as_cmap=True)
cmap.set_under(color='white')
im = ax.imshow(df3.fillna(0.0), cmap=cmap, vmin = 1.0)
ax.set_xticks(np.arange(len(df3.columns.to_list())), labels=df3.columns.to_list())
ax.set_yticks(np.arange(len(df3.index.to_list())), labels=df3.index.to_list())
ax.tick_params(axis='x', labelrotation = 90)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = ax.figure.colorbar(im, cax=cax, shrink=0.5,pad=0.01)
cbar.set_label(label='Number of Observations',labelpad=1)
cursor = mplcursors.cursor(hover=True).connect("add", lambda sel: sel.annotation.set_text(df3.index.to_list()[sel.index[0]]+'-'+df3.columns.to_list()[sel.index[1]]))
ax.set_ylabel('Atom 1')
ax.set_xlabel('Atom 2')
ax.set_title('Total Unused Peaks')
plt.draw()

df = pd.DataFrame(index=Sequence, columns = Sequence)
for upl in Used_aupls:
	tvalue = df.loc[upl.split()[0],upl.split()[2]]
	if pd.isna(tvalue): tvalue = 0.0
	df.loc[upl.split()[0],upl.split()[2]] = tvalue + 1
df = df[df > 2]
df2 = df.dropna(axis=0, how= 'all')
df3 = df2.dropna(axis=1, how= 'all')
fig1, ax = plt.subplots(figsize=(22,20))
cmap = sb.color_palette("inferno_r", as_cmap=True)
cmap.set_under(color='white')
im = ax.imshow(df3.fillna(0.0), cmap=cmap, vmin = 1.0)
ax.set_xticks(np.arange(len(df3.columns.to_list())), labels=df3.columns.to_list())
ax.set_yticks(np.arange(len(df3.index.to_list())), labels=df3.index.to_list())
ax.tick_params(axis='x', labelrotation = 90)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.05)
cbar = ax.figure.colorbar(im, cax=cax, shrink=0.2,pad=0.01)
cbar.set_label(label='Number of Observations',labelpad=1)
ax.set_ylabel('Atom 1')
ax.set_xlabel('Atom 2')
ax.set_title('Total Used peaks')
cursor = mplcursors.cursor(hover=True).connect("add", lambda sel: sel.annotation.set_text(df3.index.to_list()[sel.index[0]]+'-'+df3.columns.to_list()[sel.index[1]]))
plt.subplots_adjust(bottom=0.2, right=0.85, top= 0.85)
plt.tight_layout()
plt.draw()
plt.show()
# pdf.close()

unusedupl.close()


