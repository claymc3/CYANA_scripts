import os
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
import matplotlib.colors as mplcolors
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.linewidth'] = 1.0
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['xtick.labelsize'] = mpl.rcParams['ytick.labelsize']=8
# mpl.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.default'] = 'regular'
# mpl.rcParams['mathtext.sf']= 'sans\\-serif'
# for x in np.arange(2.5,362.5,5):
# 	print(x)
bcmap = mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF'])
gcmap = mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C'])
pcmap = mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF'])
ycmap = mplcolors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F'])
# exit()
boxes = {
'Arginine (ARG,R)':[(37.0,91.0,145.0,215.0),(151.0,220.0,138.0,226.0),(258.0,323.0,138.0,235.0),(267.0,325.0,250.0,325.0),(155.0,210.0,41.0,99.0)],
'Asparagine (ASN,N)':[(42.0,83.0,0.0,105),(42.0,83.0,263.0,359.0),(161.0,223.0,0.0,95.0),(161.0,223.0,200.0,359.0),(262.0,323.0,0,190.0),(262.0,323.0,243.0,359.0)],
'Glutamine (GLN,Q)':[(40.0,87.0,152.0,212.0),(155.0,215.0,141.0,216.0),(261.0,323.0,138.0,223.0),(265.0,323.0,255.0,338.0),(265.0,310.0,65.0,100.0),(154.0,211.0,40.0,96.0)],
'Glutamic Acid (GLU,E)':[(41.0,88.0,153.0,211.0),(153.0,215.0,141.0,215.0),(260.0,322.0,141.0,226.0),(266.0,322.0,259.0,324.0),(266.0,322.0,57.0,107.0),(154.0,212.0,38.0,96.0)],
'Histadine (HIS,H)':[(36.0,95.0,229.0,323.0),(36.0,95.0,49.0,126.0),(150.0,214.0,14.0,325.0),(260.0,330.0,33.0,359.0)],
'Isoleucine (ILE,I)':[(40.0,80.0,145.0,195.0),(168.0,212.0,145.0,185.0),(275.0,325.0,138.0,198.0), (278.0,325.0,274.0,322.0),(172.0,210.0,48.0,82.0)],
'Leucine (LEU,L)':[(155.0,212.0,40.0,88.0), (264.0,324.0,143.0,205.0)],
'Lysine (LYS,K)':[(42.0,88.0,154.0,210.0),(153.0,215.0,139.0,215.0),(261.0,321.0,139.0,225.0),(269.0,325.0,253.0,325.0),(156.0,205.0,41.0,96.0)],
'Methionine (MET, M)':[(45.0,84.0,157.0,210.0),(157.0,208.0,150.0,208.0),(265.0,322.0,146.0,212.0),(266.0,320.0,270.0,325.0),(161.0,208.0,42.0,86.0)],
'Tryptophan (TRP,W)':[(31.0,90.0,59.0,113.0),(31.0,90.0,244.0,315.0),(151.0,212.0,0.0,110.0),(151.0,212.0,216.0,360.0),(255.0,325.0,0.0,142.0),(255.0,325.0,246.0,360.0)],
'Phenylalanine/Tyrosine':[(40.0,85.0,70.0,112.0),(40.0,85.0,249.0,290.0),(155.0,205.0,21.0,114.0),(155.0,205.0,200.0,291.0),(267.0,320.0,60.0,180.0),(267.0,320.0,240.0,358.0)],
'Aspartate (ASP,D)':[(41.0,86.0,0.0,360.0),(160.0,220.0,0.0,360.0),(260.0,318.0,0.0,360.0)]}

pdf = PdfPages('Rotamer_Chi1-Chi2_templates.pdf')
files = [
['rota2018-arg.data', gcmap, [0, 0.02, 0.2, 1], 'Arginine (ARG,R)'],
['rota2018-asn.data', gcmap, [0, 0.02, 0.2, 1], 'Asparagine (ASN,N)'],
['rota2018-gln.data', gcmap, [0, 0.02, 0.2, 1], 'Glutamine (GLN,Q)'],
['rota2018-glu.data', ycmap, [0, 0.02, 0.2, 1], 'Glutamic Acid (GLU,E)'],
['rota2018-his.data', pcmap, [0, 0.02, 0.2, 1], 'Histadine (HIS,H)'],
['rota2018-ile.data', bcmap, [0, 0.02, 0.3, 1], 'Isoleucine (ILE,I)'],
['rota2018-leu.data', bcmap, [0, 0.02, 0.2, 1], 'Leucine (LEU,L)'],
['rota2018-lys.data', gcmap, [0, 0.02, 0.2, 1], 'Lysine (LYS,K)'],
['rota2018-met.data', bcmap, [0, 0.04, 0.2, 1], 'Methionine (MET, M)'],
['rota2018-trp.data', pcmap, [0, 0.02, 0.2, 1], 'Tryptophan (TRP,W)']]
for file,cmap, bounds, label in files:
	print(file)
	data = np.full((360, 360), 0, dtype=np.float64)
	linedict = {}
	with open(file) as fn:
		for line in fn:
			if line[0] != "#":
				x = int(float(line.split()[1]))
				y = int(float(line.split()[0]))
				for nx in np.arange(x-1,x+1,1):
					for ny in np.arange(y-1,y+1,1):
						data[nx][ny] = float(line.split()[-1])
	fig, ax =plt.subplots(figsize=(3,3))
	plt.imshow(data, cmap=cmap, norm=colors.BoundaryNorm(bounds, cmap.N),extent=(0, 360, 360, 0))
	ax.set_xlabel(r'$\chi$1')
	ax.set_ylabel(r'$\chi$2')
	ax.set_title(label)
	ax.set_xlim([0,360])
	print(label)
	if label in boxes.keys():
		for (x1,x2,y1,y2) in boxes[label]:
			ax.plot([x1, x1], [y1, y2], color="black",linewidth = 1.0)
			ax.plot([x2, x2], [y1, y2], color="black",linewidth = 1.0)
			ax.plot([x1, x2], [y1, y1], color="black",linewidth = 1.0)
			ax.plot([x1, x2], [y2, y2], color="black",linewidth = 1.0)
	ax.set_xticks([0,60,120,180,240,300,360])
	ax.set_yticks([0,60,120,180,240,300,360])
	ax.set_ylim([0,360])
	ax.grid(visible=True, which='major', axis='both',linestyle='--')
	plt.tight_layout()
	pdf.savefig()
	#plt.close()

files = [
['rota2018-asp.data', ycmap, [0, 0.02, 0.2, 1],'Aspartate (ASP,D)'],
['rota2018-phetyr.data', pcmap, [0, 0.05, 0.2, 1], 'Phenylalanine/Tyrosine']]
for file,cmap, bounds, label in files:
	print(file)
	data = np.full((361, 361), 0, dtype=np.float64)
	linedict = {}
	with open(file) as fn:
		for line in fn:
			if line[0] != "#":
				x = int(float(line.split()[1]))
				y = int(float(line.split()[0]))
				for nx in np.arange(x-1,x+1,1):
					for ny in np.arange(y-1,y+1,1):
						data[nx][ny] = line.split()[-1]
						data[nx+180][ny] = line.split()[-1]
	fig, ax =plt.subplots(figsize=(3,3))
	plt.imshow(data, cmap=cmap, norm=colors.BoundaryNorm(bounds, cmap.N),extent=(0, 360, 360, 0))
	ax.set_xlabel(r'$\chi1$')
	ax.set_ylabel(r'$\chi2$')
	ax.set_title(label)
	ax.set_xlim([0,360])
	if label in boxes.keys():
		for (x1,x2,y1,y2) in boxes[label]:
			ax.plot([x1, x1], [y1, y2], color="black",linewidth = 1.0)
			ax.plot([x2, x2], [y1, y2], color="black",linewidth = 1.0)
			ax.plot([x1, x2], [y1, y1], color="black",linewidth = 1.0)
			ax.plot([x1, x2], [y2, y2], color="black",linewidth = 1.0)
	ax.set_xticks([0,60,120,180,240,300,360])
	ax.set_yticks([0,60,120,180,240,300,360])
	ax.set_ylim([0,360])
	ax.grid(visible=True, which='major', axis='both',linestyle='--')
	plt.tight_layout()
	pdf.savefig()
plt.show()
	#plt.close()
pdf.close()



 








