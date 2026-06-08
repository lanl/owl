import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mplt
from matplotlib import rcParams
from matplotlib.ticker import (FormatStrFormatter, MultipleLocator)
from matplotlib.colors import LightSource

# Fonts
basefamily = 'sans-serif'
basefont = 'Arial'
fontset = 'custom'
rcParams['font.family'] = basefamily
rcParams['font.' + basefamily] = basefont
mplt.rcParams['mathtext.fontset'] = fontset
mplt.rcParams['mathtext.rm'] = basefont
mplt.rcParams['mathtext.sf'] = basefont
mplt.rcParams['mathtext.it'] = basefont + ':italic'
mplt.rcParams['mathtext.bf'] = basefont + ':bold'

def read_array(filename, shape, dtype=np.float32, ascii=False, totorch=False):

    if ascii:
        w = np.zeros(shape, dtype=dtype)
        f = open(filename, 'r')
        l = 1
        while l <= shape[0]:
            w[l - 1, :] = f.readline().strip().split()
            l = l + 1
        f.close()

    else:
        if type(shape) is tuple:
            w = np.fromfile(filename, count=np.prod(shape), dtype=dtype)
            w = np.reshape(w, shape[::-1])
            w = np.transpose(w)
        else:
            w = np.fromfile(filename, count=shape, dtype=dtype)

    if totorch:
        w = torch.from_numpy(w).type(torch.FloatTensor)

    return w


################################################################################
nx = 818
ny = 640

topo = read_array('./model/topo.bin', (ny, nx))
# topo = np.transpose(topo)
topo = np.flipud(topo)

from matplotlib.colors import LightSource

ls = LightSource(azdeg=360, altdeg=45)
rgb = ls.hillshade(topo, vert_exag=0.2)

xmin = 342855.56048806757 - (13 - 1) * 50
xmax = xmin + (nx - 1) * 50
ymin = 3956991.4307338847 - (14 - 1) * 50
ymax = ymin + (ny - 1) * 50

print(xmin, ymin)

fig, ax = plt.subplots(figsize=(5, 7.5))

plt.imshow(rgb, cmap='gray', extent=[xmin, xmax, ymin, ymax])
plt.xlim((xmin, xmax))
plt.ylim((ymin, ymax))

plt.xlabel("UTM NAD83 13N Easting (m)")
plt.ylabel("UTM NAD83 13N Northing (m)")

ax.grid(which='major', axis='both', linestyle='dashed', color='gray', linewidth=0.5)
ax.xaxis.set_major_formatter(FormatStrFormatter('%6d'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%6d'))
ax.xaxis.set_major_locator(MultipleLocator(10000))
ax.xaxis.set_minor_locator(MultipleLocator(1000))
ax.yaxis.set_major_locator(MultipleLocator(5000))
ax.yaxis.set_minor_locator(MultipleLocator(1000))
ax.set_axisbelow(True)

plt.scatter(150 * 50.0 + xmin, 550 * 50.0 + ymin, s=50, edgecolor=None, marker='*', facecolor='r', label='Source')
ax.legend(loc='lower right')

rx = np.asarray([(8 - 1) * 4000.0 + 2000.0, (6 - 1) * 4000.0 + 2000.0, (10 - 1) * 4000.0 + 2000.0])
ry = np.asarray([(9 - 1) * 3200.0 + 1600.0, (4 - 1) * 3200.0 + 1600.0, (4 - 1) * 3200.0 + 1600.0])

plt.scatter(rx + xmin, ry + ymin, s=10, edgecolor=None, marker='v', facecolor='lime', label='Receiver')
ax.legend(loc='lower right')

plt.tight_layout()
plt.savefig('./topo.pdf', dpi=300, bbox_inches='tight', pad_inches=2.0 / 72.0)
