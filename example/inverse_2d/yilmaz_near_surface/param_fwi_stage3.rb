
nx = 291
nz = 61
dx = 0.333333
dz = 0.333333

dt = 0.25e-4
tmax = 0.5

which_medium = elastic-iso

ns = 49
file_geometry = geometry_fwi/geometry.txt

model_update = vs
model_aux = vp
file_vp = model/vp_fatt.bin

dir_record = data_processed

process_record = top_mute, freq_filt, rms_balance
process_synthetic = top_mute, freq_filt, rms_balance
dp_top_mute_vel = 450
dp_top_mute_width = 0.0
dp_top_mute_taper = 0.005
dp_top_mute_shift = -0.005

data_name = z

yn_free_surface = y
free_surface_dz_refine = 5

yn_energy_precond = y

min_vs = 150
max_vs = 4000
step_max_vs = 1:25, 20:5

verbose = y

# stage 3
process_grad = smooth, rms_balance_z
dp_freq_filt_freqs = 5.0, 10.0, 50.0, 55.0
dp_freq_filt_coefs = 0.0, 1.0, 1.0, 0.0
grad_smooth_x = 1:5.0, 50:2.5
grad_smooth_z = 0.5
misfit_type = corr
file_vs = ./test_s2/iteration_30/model/updated_vs.bin
dir_working = test_s3


