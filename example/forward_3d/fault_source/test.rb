
# Create fault source
system "python create_source.py"

# Create test
system "x_runf90 create_test.f90"

# Run test
system "mpirun -np 48 owl_modeling3 param_modeling.rb "
