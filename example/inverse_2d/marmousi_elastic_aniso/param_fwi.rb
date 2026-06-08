
nx = 401
nz = 111
dx = 20
dz = 20

dt = 1.5e-3
tmax = 10

ns = 80
file_geometry = geometry/geometry.txt

dir_record = data

which_medium = elastic-tti
anisotropy_type = a-t

model_update = vp, vs, epsilon, eta
file_vp = model/vp_init.bin
file_vs = model/vs_init.bin
file_epsilon = model/eps_init.bin
file_eta = model/eta_init.bin

model_aux = theta
file_theta = model/the.bin

process_grad = mask, smooth, mask
grad_smooth_x = 1:60, 100:30
grad_smooth_z = 1:40, 100:20
grad_mask = model/mask.bin

min_vp = 950
max_vp = 5000
min_vs = 600
max_vs = 3000
min_epsilon = 0.15
max_epsilon = 0.35
min_eta = 0.05
max_eta = 0.2

step_max_vp = 100
step_max_vs = 100
step_max_epsilon = 0.01
step_max_eta = 0.01

yn_energy_precond = y

niter_max = 100

verbose = y


