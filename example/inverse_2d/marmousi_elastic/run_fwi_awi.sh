
p=param_adaptive.rb
cp -rp param_fwi.rb $p
echo 'misfit_type = adaptive ' >> $p
echo 'dir_working = test_adaptive ' >> $p
echo 'tlag_max = 1.0' >> $p
echo 'deconv_eps = 0.2' >> $p
$HOME/intel/mpi/bin/mpirun -np 80 $HOME/bin/owl_fwi2 $p
