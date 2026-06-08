
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


# DTW amplitude (vector form)
p=param_dtw_vector.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = dtw-vector ' >> $p
echo 'dir_working = test_dtw_vector ' >> $p
echo 'adj_nt = 1000 ' >> $p
echo 'dtw_form = amplitude' >> $p
echo 'tlag_max = 1.0' >> $p
echo 'dtw_smooth_median = 1' >> $p
echo 'dtw_smooth_gaussian = 3' >> $p
echo 'dtw_rinst = 2.0 ' >> $p
echo 'dtw_rcuml = 10.0 ' >> $p
echo 'dtw_loss = l0.5 ' >> $p
mpirun -np 80 owl_fwi2 $p


#####################################################################################
# Plot results
ruby plot.rb
