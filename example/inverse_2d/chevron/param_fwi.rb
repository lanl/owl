
nx = 3820
nz = 480
dx = 12.5
dz = 12.5
dt = 1.25e-3
tmax = 8

ns = 1600
file_geometry = ./geometry/chevron_geometry.txt

dir_record = ./field_data/csg

module_update = vp

yn_adpx = y
yn_free_surface = y

min_vp = 1500
max_vp = 4500

step_max_vp = 1:100, 20:25

process_grad = mask, smooth, mask
grad_mask = model/mask.bin

yn_energy_precond = y
src_index = 1, 3, 1600

sx_min = 5200
