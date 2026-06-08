
nx = 401
nz = 111
dx = 20
dz = 20

dt = 1.5e-3
tmax = 6

ns = 80
file_geometry = geometry/geometry.txt

file_vp = model/vp_init.bin

dir_record = data

model_update = vp

process_grad = mask, smooth, mask
grad_smooth_x = 1:60, 100:30
grad_smooth_z = 1:40, 100:20
grad_mask = model/mask.bin

min_vp = 1000
max_vp = 5000

yn_energy_precond = y

niter_max = 200
