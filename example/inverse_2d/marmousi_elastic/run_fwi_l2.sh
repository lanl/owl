
p=param_waveform.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = waveform ' >> $p
echo 'dir_working = test_waveform ' >> $p
$HOME/intel/mpi/bin/mpirun -np 80 $HOME/bin/owl_fwi2 $p
