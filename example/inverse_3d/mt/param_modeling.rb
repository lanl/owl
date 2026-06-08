
nx = 101
ny = 101
nz = 101
dx = 10
dy = 10
dz = 10

dt = 1.0e-3
tmax = 1

ns = 1
file_geometry = ./geometry/geometry.txt

which_medium = elastic-iso
model_name = vp, vs

file_vp = model/vp.bin
file_vs = model/vs.bin

snaps = 0.0, 0.1, 0.5

dir_synthetic = data

verbose = y

ngroup = 1
rankx = 2
ranky = 2
rankz = 2

