import os
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
import matplotlib.colors as mplcolors
mpl.rcParams['font.sans-serif'] = 'arial'
mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.linewidth'] = 1.0
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['xtick.labelsize'] = mpl.rcParams['ytick.labelsize']=8
# for x in np.arange(2.5,362.5,5):
# 	print(x)
bcmap = mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF'])
gcmap = mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C'])
pcmap = mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF'])
ycmap = mplcolors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F'])
# exit()

pdf = PdfPages('Rotamer_Chi1-Chi2_templates.pdf')
files = [
['rota2018-arg.data', gcmap, [0, 0.02, 0.2, 1], 'Arginine (ARG,R)'],
['rota2018-asn.data', gcmap, [0, 0.02, 0.2, 1], 'Asparagine (ASN,N)'],
['rota2018-gln.data', gcmap, [0, 0.02, 0.2, 1], 'Glutamine (GLN,Q)'],
['rota2018-glu.data', ycmap, [0, 0.02, 0.2, 1], 'Glutamic Acid (GLU,E)'],
['rota2018-his.data', pcmap, [0, 0.02, 0.2, 1], 'Histadine (HIS,H)'],
['rota2018-ile.data', bcmap, [0, 0.03, 0.3, 1], 'Isoleucine (ILE,I)'],
['rota2018-leu.data', bcmap, [0, 0.04, 0.2, 1], 'Leucine (LEU,L)'],
['rota2018-lys.data', gcmap, [0, 0.02, 0.2, 1], 'Lysince (LYS,K)'],
['rota2018-met.data', bcmap, [0, 0.04, 0.2, 1], 'Methionine (MET, M)'],
['rota2018-trp.data', pcmap, [0, 0.02, 0.2, 1], 'Tryptophan (TRP,W)']]
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
						data[nx][ny] = float(line.split()[-1])
	fig, ax =plt.subplots(figsize=(3,3))
	plt.imshow(data, cmap=cmap, norm=colors.BoundaryNorm(bounds, cmap.N),extent=(0, 360, 360, 0))
	ax.set_xlabel(r'$\chi$1')
	ax.set_ylabel(r'$\chi$2')
	ax.set_title(label)
	ax.set_xlim([0,360])
	ax.set_xticks([0,60,120,180,240,300,360])
	ax.set_yticks([0,60,120,180,240,300,360])
	ax.set_ylim([0,360])
	ax.grid(visible=True, which='major', axis='both',linestyle='--')
	plt.tight_layout()
	pdf.savefig()
	plt.close()

files = [
['rota2018-asp.data', ycmap, [0, 0.02, 0.2, 1],'Aspartic Acid (ASP,D) '],
['rota2018-phetyr.data', pcmap, [0, 0.02, 0.2, 1], 'Phenylalanine/Tyrosine (PHE,F/TRY,Y)']]
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
	ax.set_xlabel(r'$\chi$1')
	ax.set_ylabel(r'$\chi$2')
	ax.set_title(label)
	ax.set_xlim([0,360])
	ax.set_xticks([0,60,120,180,240,300,360])
	ax.set_yticks([0,60,120,180,240,300,360])
	ax.set_ylim([0,360])
	ax.grid(visible=True, which='major', axis='both',linestyle='--')
	plt.tight_layout()
	pdf.savefig()
	plt.close()
pdf.close()


