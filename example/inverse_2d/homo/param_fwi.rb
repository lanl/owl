
nx = 201
nz = 201
dx = 10
dz = 10

dt = 1.0e-3
tmax = 1.5

ns = 1
file_geometry = ./geometry/geometry.txt

niter_max = 1

process_grad = smooth
grad_smooth_x = 30
grad_smooth_z = 30

model_update = vp
file_vp = model/vp_init.bin
