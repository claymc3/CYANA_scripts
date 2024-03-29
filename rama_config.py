import matplotlib.colors as mplcolors
import os
import numpy as np

RAMA_PREFERENCES = {
    "General": {
        "file": os.path.join('top2018_DataSet/rama/', 'rama2018-general-noGPpreP.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.001, 0.02, 1],
    },
    "GLY": {
        "file": os.path.join('top2018_DataSet/rama/', 'rama2018-gly-sym.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F']),
        "bounds": [0, 0.002, 0.02, 1],
    },
    "PRO": {
        "file": os.path.join('top2018_DataSet/rama/', 'rama2018-pro.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
    },
    "PRE-PRO": {
        "file": os.path.join('top2018_DataSet/rama/', 'rama2018-prepro-noGP.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.002, 0.02, 1],
    }
}

ROTA_PREFERENCES = {
    "R": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-arg.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.02, 0.2, 1],
    },
    "N": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-asn.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.02, 0.2, 1],
    },
    "D": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-asp.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F']),
        "bounds": [0, 0.02, 0.2, 1],
    },
    "Q": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-gln.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.02, 0.2, 1],
    },
    "E": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-glu.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F']),
        "bounds": [0, 0.02, 0.2, 1],
    },
    "H": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-his.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.02, 0.2, 1],
    },
    "I": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-ile.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.02, 0.2, 1],
    },
    "L": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-leu.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.02, 0.2, 1],
    },
    "K": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-lys.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.02, 0.2, 1],
    },
    "M": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-met.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.04, 0.2, 1],
    },
    "F": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-phetyr.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.05, 0.2, 1],
    },
    "Y": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-phetyr.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.05, 0.2, 1],
    },
    "W": {
        "file": os.path.join('top2018_DataSet/rota/', 'rota2018-trp.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.02, 0.2, 1],
    }
}


