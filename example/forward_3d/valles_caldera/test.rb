
system "mkdir -p ./model"

# Prepare topography
system "python build_topo.py "

# Prepare model
system "x_runf90 build_model.f90 "

# Run simulation
system "mpirun -np 256 owl_modeling3 param_modeling.rb "

abort
# Download results from HPC to local

# Plot result
system "mkdir -p plot"

opts = "-n1=121 -n2=640 -d1=-50 -d2=50 -d3=50 -label1='Elevation (m)' -label2='UTM NAD83 13N Northing (m)' \
	-label3='UTM NAD83 13N Easting (m)' -tick1d=1000 -tick2d=8000 -tick3d=5000 -mtick2=4 -mtick3=4 \
	-legend=y -mask=model/mask.bin -size1=2.2 -size2=5 -size3=6 -tick2format='%6d' -tick3format='%6d' \
	-o1=3.5180000E+03 -o3=342255.56048806757 -o2=3956341.4307338847 -tick1d=-1000 -mtick1=1 \
	-slice1=2500 "

system "x_showslice -in=model/vp_masked.bin -color=rainbowcmyk -lmtick=9 -ld=1000 -unit='P-Wave Velocity (m/s)' \
	-cmin=3000 -cmax=6000 #{opts} -tr=model/vp3d.png -out=plot/vp.png &"

system "x_showslice -in=snapshot/shot_1_forward_wavefield_z_25.bin -color=bwr -lmtick=9 -unit='Amplitude' \
	-clip=0.75e-8 #{opts} -slice2=3982000 -slice3=370000 -tr=snapshot/wave3d.png -out=plot/wave.png &"

system "x_showwiggle -in=./data/data.bin -n1=9 -transpose=y -along=2 -size1=5 -size2=7 -clip=3e-8 -d2=1.0e-3 -label1='Receiver/Component'  -ticks1=0:'Cerro Toledo - X',1:'Cerro Toledo - Y',2:'Cerro Toledo - Z',3:'Valles Caldera - X',4:'Valles Caldera - Y',5:'Valles Caldera - Z',6:'Los Alamos - X',7:'Los Alamos - Y',8:'Los Alamos - Z' -label2='Time (s)' -tracecolor=k,k,k,b,b,b,r,r,r -out=plot/data.pdf &"

system "mkdir -p snapshot "
for i in 1..61

	system "x_showslice -in=snapshot/shot_1_forward_wavefield_z_#{i}.bin -color=bwr -lmtick=9 -unit='Amplitude' \
		-clip=0.75e-8 #{opts} -slice2=3982000 -slice3=370000 -out=snapshot/wave_#{i}.png -title='Time = #{(i - 1)*0.5} s' -titlex=5.5 -titley=5.15 "

	puts i

end
