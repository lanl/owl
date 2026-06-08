
##############################################################################3
# Create test and model

x_runf90 ./create_test.f90

##############################################################################3
# Forward modeling

export OMP_NUM_THREADS=6
mpirun -np 8 owl_modeling3 param_modeling.rb

##############################################################################3
# FWI

export OMP_NUM_THREADS=6
mpirun -np 8 owl_fwi3 param_fwi.rb
