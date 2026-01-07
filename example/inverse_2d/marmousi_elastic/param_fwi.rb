
nx = 401
nz = 111
dx = 20
dz = 20

dt = 1.5e-3
tmax = 10

ns = 80
file_geometry = geometry/geometry.txt

dir_record = data

which_medium = elastic-iso
model_update = vp, vs
file_vp = model/vp_init.bin
file_vs = model/vs_init.bin

process_grad = mask, smooth, mask
grad_smooth_x = 1:60, 100:30
grad_smooth_z = 1:40, 100:20
grad_mask = model/mask.bin

min_vp = 950
max_vp = 5000
min_vs = 600
max_vs = 3000

step_max_vp = 100
step_max_vs = 100

yn_energy_precond = y

niter_max = 100

verbose = y

