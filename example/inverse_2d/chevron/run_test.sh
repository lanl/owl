
######################################################################################
# Extract Vp and common-shot gathers from original SEG-Y files

python extract_data.py --vp-segy=SEG14.Vpsmoothstarting.segy --vp-raw=vp_raw.bin --data-segy=SEG14.Pisoelastic.segy --csg-dir=./csg --shot-key=FieldRecord


######################################################################################
# Create model and geometry for OWL

make clean
make
./create_test


######################################################################################
# FWI

export OMP_NUM_THREADS=32

# Stage a
p=param.rb
cp -rp param_fwi.rb $p
echo 'file_vp = model/vp.bin ' >> $p
echo 'process_record = freq_filt ' >> $p
echo 'process_synthetic = freq_filt ' >> $p
echo 'grad_smooth_x = 1:500, 50:250' >> $p
echo 'grad_smooth_z = 1:250, 50:125' >> $p
echo 'misfit_type = dtw ' >> $p
echo 'adj_nt = 1001 ' >> $p
echo 'tlag_max = 1 ' >> $p
echo 'dtw_smooth_median = 2' >> $p
echo 'dtw_smooth_gaussian = 3' >> $p
echo 'dtw_rinst = 15.0 ' >> $p
echo 'dtw_rcuml = 15.0 ' >> $p
echo 'dtw_loss = l0.5 ' >> $p
echo 'dtw_form = phase' >> $p
echo 'dp_freq_filt_freqs = 1.0, 2.0, 5.0, 6.0' >> $p
echo 'dp_freq_filt_coefs = 0.0, 1.0, 1.0, 0.0' >> $p
echo 'dir_working = test_a ' >> $p
mpirun -np 200 owl_fwi2 $p


# Stage b
p=param.rb
cp -rp param_fwi.rb $p
echo 'file_vp = test_a/iteration_30/model/updated_vp.bin ' >> $p
echo 'process_record = freq_filt ' >> $p
echo 'process_synthetic = freq_filt ' >> $p
echo 'process_grad = mask, smooth, mask' >> $p
echo 'grad_smooth_x = 1:300, 50:150' >> $p
echo 'grad_smooth_z = 1:75, 50:25' >> $p
echo 'grad_mask = model/mask.bin' >> $p
echo 'misfit_type = corr'  >> $p
echo 'dp_freq_filt_freqs = 1.0, 2.0, 7.0, 8.0' >> $p
echo 'dp_freq_filt_coefs = 0.0, 1.0, 1.0, 0.0' >> $p
echo 'dir_working = test_b ' >> $p
mpirun -np 200 owl_fwi2 $p


# Stage c
p=param.rb
cp -rp param_fwi.rb $p
echo 'file_vp = test_b/iteration_50/model/updated_vp.bin ' >> $p
echo 'process_record = freq_filt ' >> $p
echo 'process_synthetic = freq_filt ' >> $p
echo 'process_grad = mask, smooth, mask' >> $p
echo 'grad_smooth_x = 1:150, 50:75' >> $p
echo 'grad_smooth_z = 1:50, 50:25' >> $p
echo 'grad_mask = model/mask.bin' >> $p
echo 'misfit_type = corr'  >> $p
echo 'dp_freq_filt_freqs = 1.0, 2.0, 9.0, 10.0' >> $p
echo 'dp_freq_filt_coefs = 0.0, 1.0, 1.0, 0.0' >> $p
echo 'dir_working = test_c ' >> $p
mpirun -np 200 owl_fwi2 $p


# Stage d
p=param.rb
cp -rp param_fwi.rb $p
echo 'file_vp = test_c/iteration_100/model/updated_vp.bin ' >> $p
echo 'process_record = freq_filt ' >> $p
echo 'process_synthetic = freq_filt ' >> $p
echo 'process_grad = mask, smooth, mask' >> $p
echo 'grad_smooth_x = 1:150, 50:75' >> $p
echo 'grad_smooth_z = 1:50, 50:25' >> $p
echo 'grad_mask = model/mask.bin' >> $p
echo 'misfit_type = corr'  >> $p
echo 'dp_freq_filt_freqs = 1.0, 2.0, 13.0, 15.0' >> $p
echo 'dp_freq_filt_coefs = 0.0, 1.0, 1.0, 0.0' >> $p
echo 'dir_working = test_d ' >> $p
mpirun -np 200 owl_fwi2 $p


######################################################################################
# Plot results

ruby plot.rb

ln -s ../../../misc/python ./
python plot_misfit.py
