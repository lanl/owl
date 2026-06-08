
# Dimension parameters
nx = 201
ny = 201
nz = 81

dx = 50
dy = 50
dz = 50

dt = 4.0e-3
tmax = 10

# Medium parameters
which_medium = elastic-iso
model_name = vp, vs, rho

file_vp = ./model/vp.bin
file_vs = ./model/vs.bin
file_rho = ./model/rho.bin

# Geometry
ns = 1
file_geometry = ./geometry/geometry.txt

# Snapshots
snaps = 0.0, 0.5, 10

verbose = y

yn_free_surface = y
# For simulation without much surface wave, we can use
# a coarse/original mesh in the near-surface region.
free_surface_dz_refine = 2

dir_synthetic = data
dir_snapshot = snapshot

# MPI domain decomposition
ngroup = 1
rankx = 4
ranky = 4
rankz = 3
