import os
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
import matplotlib.colors as mplcolors
# for x in np.arange(2.5,362.5,5):
# 	print(x)
bcmap = mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF'])
gcmap = mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C'])
pcmap = mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF'])
ycmap = mplcolors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F'])
# exit()

pdf = PdfPages('Ramachandran_Phi-Psi_templates.pdf')
files = [
['rama8000-general-noGPIVpreP.data', bcmap, [0, 0.001, 0.02, 1],'General'],
['rama8000-ileval-nopreP.data', bcmap, [0, 0.002, 0.02, 1],'Ile Val General'],
['rama8000-prepro-noGP.data', pcmap, [0, 0.002, 0.02, 1],'Pre-Proline'],
['rama8000-gly-sym.data', ycmap,[0, 0.002, 0.02, 1],'Glycine'],
['rama8000-cispro.data', gcmap, [0, 0.002, 0.02, 1],'cis-Proline'],
['rama8000-transpro.data', gcmap, [0, 0.002, 0.02, 1],'trans-Proline']
]


for file,cmap,bounds,label  in files:
	fig, ax =plt.subplots(figsize=(3,3))
	data = np.full((360, 360), 0, dtype=np.float64)
	with open(file) as fn:
		for line in fn:
			if line[0] != "#":
				x = int(float(line.split()[1]))
				y = int(float(line.split()[0]))
				for nx in np.arange(x-1,x+1,1):
					for ny in np.arange(y-1,y+1,1):
						data[nx+180][ny+180] = float(line.split()[-1])
	fig, ax =plt.subplots(figsize=(3,3))
	plt.imshow(data, cmap=cmap, norm=colors.BoundaryNorm(bounds, cmap.N),extent=(-180, 180, 180, -180))
	ax.set_xlabel(r'$\phi$1')
	ax.set_ylabel(r'$\psi$2')
	ax.set_title(label)
	ax.set_xlim([-180,180])
	ax.set_xticks([-180,-120,-60,0,60,120,180])
	ax.set_yticks([-180,-120,-60,0,60,120,180])
	ax.plot([-180, 180], [0, 0], color="black")
	ax.plot([0, 0], [-180, 180], color="black")
	ax.set_ylim([-180,180])
	ax.grid(visible=True, which='major', axis='both',linestyle='--')
	plt.tight_layout()
	pdf.savefig()
	plt.close()

pdf.close()
