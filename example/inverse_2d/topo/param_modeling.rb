
nx = 301
nz = 101

dx = 10
dz = 10

dt = 1.25e-4
data_dt = 1.0e-3
tmax = 2.5

ns = 60
file_geometry = ./geometry/geometry.txt

dir_synthetic = data

which_medium = elastic-tti
anisotropy_type = iso
model_name = vp, vs, rho
file_vp = model/vp.bin
file_vs = model/vs.bin
file_rho = model/rho.bin

yn_free_surface = y
measure_source_depth_from_surface = y
measure_receiver_depth_from_surface = y
free_surface_dz_refine = 2
file_topo = ./model/ftopo.txt

verbose = y
