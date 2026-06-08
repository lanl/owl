
nx = 818
ny = 640
nz = 121

dx = 50
dy = 50
dz = 50

dt = 1.0e-3
tmax = 30

ns = 1
file_geometry = ./geometry/geometry.txt

verbose = y

which_medium = elastic-tti
anisotropy_type = iso
model_name = vp, vs, rho
file_vp = model/vp.bin
file_vs = model/vs.bin
file_rho = model/rho.bin

dir_synthetic = data
dir_snapshot = snapshot

yn_free_surface = y
free_surface_dz_refine = 3
measure_source_depth_from_surface = y
measure_receiver_depth_from_surface = y
file_topo = model/topo.txt

snaps = 0.0, 0.5, 30.0

ngroup = 1
rankx = 8
ranky = 4
rankz = 2

rankx = 8
ranky = 8
rankz = 4
