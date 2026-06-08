
import math
import numpy as np

def regtick(start, end, interval, mtick):

    tick = np.arange(start, end + 0.1 * abs(interval), interval)
    minor_tick_interval = interval / (mtick + 1.0)
    minor_tick = np.arange(start, end + 0.1 * np.abs(minor_tick_interval), minor_tick_interval)

    return tick, minor_tick


def regspace(start, stop, step=1):
    dir = 1 if (step > 0) else -1
    return np.arange(start, stop + dir, step)


def next_power_of_2(x):

    if x == 0:
        return 1
    else:
        return 2**math.ceil(math.log2(x))


def rescale(w, range=(0, 1)):

    if np.max(w) != np.min(w):

        dr = range[1] - range[0]
        w = w - np.min(w)
        w = w / np.max(w)
        w = w * dr + range[0]

    else:
        w = np.ones_like(w) * range[0]

    return w


## Convert string to bool
def str2bool(v):

    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'y', 'true', 't', 'on', '1'):
        return True
    elif v.lower() in ('no', 'n', 'false', 'f', 'off', '0'):
        return False
    else:
        print(' Error: Argument must be one of yes/no, y/n, true/false, t/f, on/off, 1/0. ')
        exit()


## convert strings to bool
def strs2bool(v):

    v = v.split(',')
    n = len(v)
    r = np.zeros(n, dtype=bool)

    for i in range(n):

        if isinstance(v[i], bool):
            r[i] = v[i]
        if v[i].lower() in ('yes', 'y', 'true', 't', 'on', '1'):
            r[i] = True
        elif v[i].lower() in ('no', 'n', 'false', 'f', 'off', '0'):
            r[i] = False
        else:
            print(' Error: Argument must be one of yes/no, y/n, true/false, t/f, on/off, 1/0. ')
            exit()

    return r


## forward integer indices
def forward_range(start, n, step=1):

    r = np.zeros(n, dtype=np.int32)
    for k in range(n):
        r[k] = start + k * step

    return r


## backward integer indices
def backward_range(end, n, step=1):

    r = np.zeros(n, dtype=np.int32)
    for k in range(n):
        r[k] = end - k * step

    return r


## strict start-end indices
def strict_range(start, end, step=1):

    if end < start:
        s = -np.abs(step)
    else:
        s = np.abs(step)
    r = np.zeros(np.int32(np.floor((end - start) / s)) + 1, dtype=np.int32)
    for k in range(np.size(r)):
        r[k] = start + k * s

    return r


## date and time
from datetime import datetime


def date_time():

    now = datetime.now()
    return now.strftime(" %Y/%m/%d %H:%M:%S ")


## convert number to string
def num2str(x, f=''):

    # v=10.4
    #
    # print('% 6.2f' % v)
    #   10.40
    #
    # print('% 12.1f' % v)
    #         10.4
    #
    # print('%012.1f' % v)
    # 0000000010.4

    if type(x) is int and f == '':
        ff = "{:d}"
    elif type(x) is float and f == 'd':
        ff = "{:.0f}"
    else:
        ff = "{:" + f + "}"

    return ff.format(x)


## Get numpy array and transfer to CPU
def get_numpy(w):

    return w.squeeze().data.cpu().numpy()


## Find indices
def find_indices(list_to_check, item_to_find):

    indices = []
    for idx, value in enumerate(list_to_check):
        if value == item_to_find:
            indices.append(idx)

    return indices


## Set font for maplotlib plotting
import matplotlib as mplt
from matplotlib import rcParams

def set_font():

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
