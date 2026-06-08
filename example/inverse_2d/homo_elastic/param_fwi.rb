
nx = 201
nz = 201
dx = 10
dz = 10

dt = 1.0e-3
tmax = 2

ns = 1
file_geometry = ./geometry/geometry.txt

niter_max = 1

which_medium = elastic-iso

process_grad = smooth
grad_smooth_x = 30
grad_smooth_z = 30

model_update = vp, vs
file_vp = model/vp_init.bin
file_vs = model/vs_init.bin

data_name = z
trace_discard_threshold = 1.0e-2
