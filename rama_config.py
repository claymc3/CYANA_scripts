import matplotlib.colors as mplcolors
import os

RAMA_PREFERENCES = {
    "General": {
        "file": os.path.join('top500-angles/pct/rama/', 'rama500-general.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.0005, 0.02, 1],
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

ROTA_PREFERENCES = {
    "I": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-ile.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.002, 0.02, 1],
        "incr": 2.5,
    },
    "L": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-leu.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 2.5,
    },
    "M": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-met.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 4.0,
    },
    "F": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-phetyr.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 2.5,
    },
    "Y": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-phetyr.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 2.5,
    },
    "W": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-trp.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 2.5,
    },
    "H": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-his.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 2.5,
    },
    "R": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-arg.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 5.0,
    },
    "D": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-asp.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 2.5,
    },
    "E": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-glu.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 4.0,
    },
    "N": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-asn.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 2.5,
    },
    "Q": {
        "file": os.path.join('top500-angles/pct/rota/', 'rota500-gln.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 4.0,
    },
    "K": {
        "file": os.path.join('top500-angles/pct/rama/', 'rota500-lys.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
         "incr": 5.0,
    }
}



