
system "x_runf90 create_test.f90"

system "mpirun -np 42 owl_modeling3 param_modeling.rb "
