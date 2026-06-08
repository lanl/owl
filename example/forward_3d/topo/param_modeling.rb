
nx = 200
ny = 35
nz = 40

dx = 10
dy = 10
dz = 10

dt = 0.5e-3
tmax = 1.5
data_dt = 1.0e-3

ns = 1

file_geometry = ./geometry/geometry.txt

snaps = 0.0, 0.1, 2.0

verbose = y

yn_free_surface = y
free_surface_dz_refine = 4
measure_source_depth_from_surface = y
measure_receiver_depth_from_surface = y
file_topo = ftopo.txt

f0_factor = 1.2

which_medium = elastic-tti
anisotropy_type = iso
model_name = vp, vs, rho
file_vp = model/vp.bin
file_vs = model/vs.bin
file_rho = model/rho.bin
yn_free_surface = y

dir_synthetic = data_fsg
dir_snapshot = snapshot_fsg

ngroup = 1
rankx = 7
ranky = 3
rankz = 2
