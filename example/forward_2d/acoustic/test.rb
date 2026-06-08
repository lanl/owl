
# Generate model and geometry
system "x_runf90 create_test.f90"


# Modeling
system "mpirun -np 1 owl_modeling2 param_modeling.rb"


# Plot results

system "mkdir -p ./plot"

system "x_showmatrix -in=model/vp.bin -n1=100 -d1=0.01 -d2=0.01 -label1='Z (km)' -label2='X (km)' -legend=y -unit='P-wave Velocity (m/s)' -tick1d=0.5 -tick2d=0.5 -mtick1=4 -mtick2=4 -out=plot/vp.png &"
system "x_showmatrix -in=model/rho.bin -n1=100 -d1=0.01 -d2=0.01 -label1='Z (km)' -label2='X (km)' -legend=y -unit='Density (kg/m$^3$)' -tick1d=0.5 -tick2d=0.5 -mtick1=4 -mtick2=4 -out=plot/rho.png &"

for i in 1..16

    system "x_showmatrix -in=snapshot/shot_1_forward_wavefield_p_#{i}.bin -n1=100 -d1=0.01 -d2=0.01 -label1='Z (km)' -label2='X (km)' -color=binary -clip=10 -tick1d=0.5 -tick2d=0.5 -mtick1=4 -mtick2=4 -out=plot/snapshot_#{i}.png &"

    puts i

end
