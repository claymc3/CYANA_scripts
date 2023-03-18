import matplotlib.colors as mplcolors
import os
import numpy as np

RAMA_PREFERENCES = {
    "General": {
        "file": os.path.join('top500-angles/pct/rama/', 'rama500-general.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.01, 0.02, 1],
    },
    "GLY": {
        "file": os.path.join('top500-angles/pct/rama/', 'rama500-gly-sym.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F']),
        "bounds": [0, 0.002, 0.02, 1],
    },
    "PRO": {
        "file": os.path.join('top500-angles/pct/rama/', 'rama500-pro.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
    },
    "PRE-PRO": {
        "file": os.path.join('top500-angles/pct/rama/', 'rama500-prepro.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.002, 0.02, 1],
    }
}
ILWNDH = []
count = 5 
trans_dict = {}
for x in np.arange(0,360,5):
	for y in np.arange(0,360,5):
		count+=1
		for nx in np.arange(x,x+5,1):
			for ny in np.arange(y,y+5,1):
				ILWNDH.append((count, nx, ny))
YFD = []
count = 5 
for x in np.arange(0,360,5):
	for y in np.arange(0,180,5):
		count+=1
		for nix in np.arange(x,x+5,1):
			for ny in np.arange(y,y+5,1):
				YFD.append((count, nx,ny))
				YFD.append((count, nx,ny+180))
MQ = []
count = 6 
trans_dict = {}
for x in np.arange(0,360,8):
	for y in np.arange(0,360,8):
		for z in np.arange(0,360,8):
			count+=1
			for nx in np.arange(x,x+8,1):
				for ny in np.arange(y,y+8,1):
					MQ.append((count, nx,ny))

ROTA_PREFERENCES = {
    "I": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-ile.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.002, 0.02, 1],
        "incr": ILWNDH,
    },
    "L": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-leu.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": ILWNDH,
    },
    "M": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-met.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": MQ,
    },
    "F": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-phetyr.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": YFD,
    },
    "Y": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-phetyr.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": YFD,
    },
    "W": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-trp.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": ILWNDH,
    },
    "H": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-his.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": [[ 360.0, 5.0], [ 360.0, 5.0]],
    },
    # "R": {
    #     "file": os.path.join('top500-angles/pct/rota/', 'rota500-arg.data'),
    #     "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
    #     "bounds": [0, 0.002, 0.02, 1],
    #      "incr": [[360.0, 10.0], [360.0, 10.0], [360.0, 10.0], [360.0, 10.0]],
    # },
    "D": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-asp.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr":YFD,
    },
    # "E": {
    #     "file": os.path.join('top500-angles/pct/rota/', 'rota500-glu.data'),
    #     "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
    #     "bounds": [0, 0.002, 0.02, 1],
    #      "incr": [[360.0, 8.0], [360.0, 8.0], [180.0, 7.826086956521739]],
    # },
    "N": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-asn.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": ILWNDH,
    },
    "Q": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-gln.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": MQ,
    },
    "K": {
        "file": os.path.join('top500-angles/pct/rama/', 'rota500-lys.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": [[360.0, 10.0], [360.0, 10.0], [ 60.0, 10.0], [360.0, 10.0]],
    }
}


