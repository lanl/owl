
################################################################################
# Make model and geometry

x_runf90 create_test.f90


################################################################################
# Modeling

export OMP_NUM_THREADS=4
mpirun -np 1760 owl_modeling3 ./param_modeling.rb


################################################################################
# FWI

# L2
p=param_waveform.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = waveform ' >> $p
echo 'dir_working = test_waveform ' >> $p
echo 'jumpout_factor = 1~55:1.0, 56:1.05' >>$p
mpirun -np 1760 owl_fwi3 $p


# Adaptive
p=param_adaptive.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = adaptive ' >> $p
echo 'dir_working = test_adaptive ' >> $p
echo 'tlag_max = 0.5' >> $p
echo 'deconv_eps = 0.2' >> $p
echo 'jumpout_factor = 1~20:1, 21:1.05' >> $p
mpirun -np 1760 owl_fwi3 $p


# DTW amplitude
p=param_dtw.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = dtw ' >> $p
echo 'dir_working = test_dtw ' >> $p
echo 'adj_nt = 500 ' >> $p
echo 'dtw_form = amp' >> $p
echo 'tlag_max = 0.5' >> $p
echo 'dtw_smooth_median = 1' >> $p
echo 'dtw_smooth_gaussian = 3' >> $p
echo 'dtw_rinst = 2 ' >> $p
echo 'dtw_rcuml = 5 ' >> $p
echo 'dtw_loss = l0.5 ' >> $p
echo 'jumpout_factor = 1.05 ' >> $p
mpirun -np 1760 owl_fwi3 $p


################################################################################
# Plot results

ruby plot.rb
