
import numpy as np

def count_nonempty_lines(filename):

    file = open(filename, "r")
    line_count = 0
    for line in file:
        if line != "\n":
            line_count += 1
    file.close()

    return line_count


def read_array(filename, shape, dtype=np.float32, ascii=False):

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

    return w


def write_array(w, filename, dtype=np.float32, ascii=False):

    w = np.asarray(w, dtype=dtype)
    w = np.transpose(w)

    if ascii:
        np.savetxt(filename, w.transpose(), fmt='%e', delimiter=' ')
    else:
        w.tofile(filename)
