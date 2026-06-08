
nx = 221
ny = 201
nz = 81

dx = 20
dy = 20
dz = 20

dt = 1.5e-3
tmax = 3

ns = 110
file_geometry = ./geometry/geometry.txt

which_medium = acoustic-iso
file_vp = model/vp_init.bin

dir_record = data

model_update = vp

process_grad = mask, smooth, mask
grad_smooth_x = 1:60, 100:30
grad_smooth_z = 1:40, 100:20
grad_mask = model/mask.bin

min_vp = 1900
max_vp = 4800

yn_energy_precond = y

niter_max = 200

ngroup = 55
rankx = 4
ranky = 4
rankz = 2
