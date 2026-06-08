
nx = 300
nz = 200

dx = 10
dz = 10

dt = 0.4e-3
data_dt = 1.0e-3
tmax = 2

ns = 1
file_geometry = ./geometry/geometry.txt

snaps = 0.0, 0.1, 1.5

which_medium = elastic-tti
model_name = vp, vs, rho, epsilon, delta, theta
file_vp = model/vp.bin
file_vs = model/vs.bin
file_rho = model/rho.bin
file_epsilon = model/eps.bin
file_delta = model/del.bin
file_theta = model/the.bin

verbose = y

f0_factor = 1.2
