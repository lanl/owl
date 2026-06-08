
make clean
make

./create_test

# fatt
mpirun -np 49 x_fatt2 ./param_fatt.rb

./build_vs

# fwi
mpirun -np 49 owl_fwi2 ./param_fwi_stage1.rb
mpirun -np 49 owl_fwi2 ./param_fwi_stage2.rb
mpirun -np 49 owl_fwi2 ./param_fwi_stage3.rb

# plot
ruby plot.rb
