
#####################################################################################
# Create model and geometry
x_runf90 create_test.f90


#####################################################################################
# Modeling
export OMP_NUM_THREADS=32
mpirun -np 80 owl_modeling2 ./param_modeling.rb


#####################################################################################
# FWI

export OMP_NUM_THREADS=32

# Waveform
p=param_waveform.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = waveform ' >> $p
echo 'dir_working = test_waveform ' >> $p
mpirun -np 80 owl_fwi2 $p


# Adaptive
p=param_adaptive.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = adaptive ' >> $p
echo 'dir_working = test_adaptive ' >> $p
echo 'tlag_max = 1.0' >> $p
echo 'deconv_eps = 0.2' >> $p
mpirun -np 80 owl_fwi2 $p


# Adaptive-spacetime
p=param_adaptive_st.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = adaptive-spacetime ' >> $p
echo 'dir_working = test_adaptive_spacetime ' >> $p
echo 'tlag_max = 1.0' >> $p
echo 'deconv_eps = 0.1' >> $p
echo 'adaptive_half_window = 3 ' >> $p
mpirun -np 80 owl_fwi2 $p


# Adaptive-local
p=param_adaptive_local.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = local-adaptive ' >> $p
echo 'dir_working = test_adaptive_local ' >> $p
echo 'adj_nt = 2001' >> $p
echo 'tlag_max = 1.0' >> $p
echo 'deconv_eps = 0.1' >> $p
echo 'lawi_sigma = 0.5' >> $p
mpirun -np 80 owl_fwi2 $p


# Phase
p=param_phase.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = phase ' >> $p
echo 'dir_working = test_phase ' >> $p
mpirun -np 80 owl_fwi2 $p


# DTW amplitude
p=param_dtw_amp.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = dtw ' >> $p
echo 'dir_working = test_dtw ' >> $p
echo 'adj_nt = 501' >> $p
echo 'dtw_form = amp' >> $p
echo 'tlag_max = 0.5' >> $p
echo 'dtw_smooth_median = 1' >> $p
echo 'dtw_smooth_gaussian = 3' >> $p
echo 'dtw_rinst = 0.5 ' >> $p
echo 'dtw_rcuml = 5.0 ' >> $p
echo 'dtw_loss = l0.5 ' >> $p
echo 'jumpout_factor = 1.05 ' >> $p
mpirun -np 80 owl_fwi2 $p


# DTW phase
p=param_dtw_phase.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = dtw ' >> $p
echo 'dir_working = test_dtw_phase ' >> $p
echo 'adj_nt = 501' >> $p
echo 'dtw_form = phase' >> $p
echo 'tlag_max = 0.5' >> $p
echo 'dtw_smooth_median = 1' >> $p
echo 'dtw_smooth_gaussian = 3' >> $p
echo 'dtw_loss = l0.5 ' >> $p
echo 'dtw_rinst = 0.5 ' >> $p
echo 'dtw_rcuml = 5.0 ' >> $p
echo 'jumpout_factor = 1.05' >> $p
mpirun -np 80 owl_fwi2 $p


#####################################################################################
# Plot results
ruby plot.rb
