

nx = 300
nz = 300

dx = 10
dz = 10

dt = 1e-3
data_dt = 1.0e-3
tmax = 3.0

ns = 1
file_geometry = ./geometry/geometry.txt

snaps = 0.0, 0.2, 3.0

verbose = y

f0_factor = 1.2

which_medium = elastic-tti
anisotropy_type = cij
model_name = c11, c13, c33, c55, rho
file_c11 = model/c11.bin
file_c13 = model/c13.bin
file_c33 = model/c33.bin
file_c55 = model/c55.bin
file_rho = model/rho.bin

dir_synthetic = data_fsg
dir_snapshot = snapshot_fsg
