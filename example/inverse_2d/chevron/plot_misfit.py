
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from python.libpy_io import *
from python.libpy_utility import *
from matplotlib.ticker import (MultipleLocator)

set_font()

n1 = 30
n2 = 50
n3 = 70
n4 = 100

sa = read_array('./test_a/data_misfit.txt', (n1, 3), ascii=True)
sb = read_array('./test_b/data_misfit.txt', (n2, 3), ascii=True)
sc = read_array('./test_c/data_misfit.txt', (n3, 3), ascii=True)
sd = read_array('./test_d/data_misfit.txt', (n4, 3), ascii=True)

fig, ax = plt.subplots(1, 1, figsize=(10, 3))

ax.plot(sa[:n1 + 1, 0], sa[:n1 + 1, 2], 'v', linestyle='-', linewidth=1, markersize=3, color='k', label='Stage 1: 1-2-5-6 Hz')
ax.plot(n1 + sb[:n2 + 1, 0], sb[:n2 + 1, 2], 'v', linestyle='-', linewidth=1, markersize=3, color='b', label='Stage 2: 1-2-7-8 Hz')
ax.plot(n1 + n2 + sc[:n3 + 1, 0], sc[:n3 + 1, 2], 'v', linestyle='-', linewidth=1, markersize=3, color='green', label='Stage 3: 1-2-9-10 Hz')
ax.plot(n1 + n2 + n3 + sd[:n4 + 1, 0], sd[:n4 + 1, 2], 'v', linestyle='-', linewidth=1, markersize=3, color='r', label='Stage 4: 1-2-13-15 Hz')

plt.xlim(-2, 250)
plt.ylim(0, 1.1)
plt.legend(ncols=2, fontsize=12)
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_minor_locator(MultipleLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(0.25))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))

ax.set_ylabel('Normalized Data Misfit', fontsize=14)
ax.set_xlabel('Iteration Number', fontsize=14)

ax.tick_params(axis='both', which='major', labelsize=12)

plt.savefig('./plot/misfit.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
