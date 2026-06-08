

nx = 300
nz = 50
dx = 10
dz = 10

dt = 0.4e-3
data_dt = 1.0e-3
tmax = 2

ns = 1
file_geometry = ./geometry/geometry.txt

snaps = 0.0, 0.1, 1.5

verbose = y

f0_factor = 1.2

which_medium = elastic-tti
anisotropy_type = thomsen
model_name = vp, vs, epsilon, delta, theta, rho
file_vp = model/vp.bin
file_vs = model/vs.bin
file_epsilon = model/eps.bin
file_delta = model/del.bin
file_theta = model/the.bin
file_rho = model/rho.bin

