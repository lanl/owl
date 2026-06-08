
# Dimension parameters
nx = 70
ny = 60
nz = 50

dx = 10
dy = 10
dz = 10

dt = 0.5e-3
tmax = 2

# Medium parameters
which_medium = elastic-tti
anisotropy_type = cij
model_name = rho, c11, c12, c13, c22, c23, c33, c44, c55, c66

file_c11 = ./model/c11.bin
file_c12 = ./model/c12.bin
file_c13 = ./model/c13.bin
file_c22 = ./model/c22.bin
file_c23 = ./model/c23.bin
file_c33 = ./model/c33.bin
file_c44 = ./model/c44.bin
file_c55 = ./model/c55.bin
file_c66 = ./model/c66.bin
file_rho = ./model/rho.bin

# Geometry
ns = 1
file_geometry = ./geometry/geometry.txt

# Snapshots
snaps = 0.0, 0.1, 2.0

verbose = y

yn_free_surface = y
# For simulation without much surface wave, we can use
# a coarse/original mesh in the near-surface region.
free_surface_dz_refine = 1

dir_synthetic = data
dir_snapshot = snapshot

# MPI domain decomposition
ngroup = 1
rankx = 7
ranky = 3
rankz = 2
