
nx = 301
nz = 101

dx = 10
dz = 10

dt = 1.25e-4
data_dt = 1.0e-3
tmax = 2.5

ns = 60
file_geometry = ./geometry/geometry.txt

which_medium = elastic-tti
anisotropy_type = iso
model_update = vp, vs
model_aux = rho
file_vp = model/vp_init.bin
file_vs = model/vs_init.bin
file_rho = model/rho.bin

min_vp = 2600
max_vp = 4800
min_vs = 1500
max_vs = 2800

step_max_vp = 100
step_max_vs = 50

yn_free_surface = y
measure_source_depth_from_surface = y
measure_receiver_depth_from_surface = y
free_surface_dz_refine = 2
file_topo = ./model/ftopo.txt

process_grad = smooth, mask
grad_smooth_x = 20
grad_smooth_z = 10
grad_mask = model/mask.bin

dir_record = data

niter_max = 200

yn_energy_precond = y

verbose = y
