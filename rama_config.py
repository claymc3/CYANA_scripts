import matplotlib.colors as mplcolors
import os
import numpy as np

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
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-ile_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.03, 0.1, 1],
    },
    "L": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-leu_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.05, 0.1, 1],
    },
    "M": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-met_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.2, 0.3, 1],
    },
    "F": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-phetyr_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.05, 0.1, 1],
    },
    "Y": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-phetyr_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.05, 0.1, 1],
    },
    "W": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-trp_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.05, 0.1, 1],
    },
    "H": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-his_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#E0D1FF', '#BE9EFF']),
        "bounds": [0, 0.05, 0.1, 1],
    },
    "R": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-arg_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.05, 0.1, 1],
    },
    "D": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-asp_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.05, 0.1, 1],
    },
    "E": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-glu_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.05, 0.1, 1],
    },
    "N": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-asn_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.05, 0.1, 1],
    },
    "Q": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-gln_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.05, 0.1, 1],
    },
    "K": {
        "file": os.path.join('top500-angles/pct/rota_2/', 'rota500-lys_chi1chi2.data'),
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.05, 0.1, 1],
    }
}


