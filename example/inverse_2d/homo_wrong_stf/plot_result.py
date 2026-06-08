
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

from python.libpy_io import *
from python.libpy_utility import *
from matplotlib.ticker import (MultipleLocator)

set_font()

n1 = 201
n2 = 201

method = ['waveform', 'envelope', 'adaptive', 'local-adaptive', 'dtw_1', 'dtw_2']
label = ['L2', 'Envelope', 'AWI', 'LAWI', 'GWI-1', 'GWI-2']

nm = len(label)

fig, ax = plt.subplots(2, nm, figsize=(11, 3.5), constrained_layout=False)

for i in range(nm):
    g = read_array('./test_low_' + method[i] + '/iteration_1/model/grad_vp.bin', (n1, n2))  
    gmax = np.max(np.abs(g))*0.75
    ax[0, i].imshow(g, cmap='bwr', vmin=-gmax, vmax=gmax)

    g = read_array('./test_high_' + method[i] + '/iteration_1/model/grad_vp.bin', (n1, n2))  
    gmax = np.max(np.abs(g))*0.75
    ax[1, i].imshow(g, cmap='bwr', vmin=-gmax, vmax=gmax)

sgnp = ['-', '+']

for i in range(2):
    for j in range(nm):

        ax[i, j].tick_params(bottom=True)

        if i == 0:
            ax[i, j].xaxis.set_ticks_position('top')
            ax[i, j].xaxis.set_label_position('top')

        if j >= 1 and j < nm - 1:
            ax[i, j].set_yticklabels([])
            ax[i, j].tick_params(axis='y', which='both', right=True, labelright=True)
        if j == nm - 1:
            ax[i, j].yaxis.set_ticks_position('both')
            ax[i, j].yaxis.set_label_position('right')
            ax[i, j].tick_params(axis='y', which='both', right=True, labelleft=False, labelright=True)
            ax[i, j].set_ylabel("$" + sgnp[i] + "\\Delta V_p$", fontsize=14, fontweight="bold", rotation=0)
            ax[i, j].yaxis.labelpad = 5
            ax[i, j].yaxis.label.set_verticalalignment('center')
            ax[i, j].yaxis.label.set_horizontalalignment('left')

        ax[i, j].tick_params(axis='x', which='both', bottom=True, top=True)
        
        ax[i, j].xaxis.set_major_locator(MultipleLocator(100))
        ax[i, j].xaxis.set_minor_locator(MultipleLocator(10))
        ax[i, j].yaxis.set_major_locator(MultipleLocator(100))
        ax[i, j].yaxis.set_minor_locator(MultipleLocator(10))

        if i < 1 and i > 0:
            ax[i, j].set_xticklabels([])
        if i == 1:
            ax[i, j].set_title(label[j], y=-0.4, fontweight='bold', fontsize=14)

        ax[i, j].plot(30, 100, '*', markersize=7, color='k')
        ax[i, j].plot(170, 100, 'v', markersize=3, color='k')

# plt.tight_layout()
plt.savefig('./grad.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
