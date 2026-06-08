
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

niter_max = 1

which_medium = elastic-tti

model_aux = vp, vs
file_vp = model/vp.bin
file_vs = model/vs.bin

model_update = mt
file_mt = model/mt_init.bin

min_mt = -2
max_mt = 2

dir_record = data

niter_max = 5
dir_working = test_mt

verbose = y

ngroup = 1
rankx = 2
ranky = 2
rankz = 2

