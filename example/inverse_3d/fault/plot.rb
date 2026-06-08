

opts = "-n1=81 -n2=201 -label1='Depth (km)' -size1=2 -size2=3.2 -size2=3.5 -label2='Horizontal Position Y (km)' -label3='Horizontal Position X (km)' -d1=0.02 -d2=0.02 -d3=0.02 -tick1d=0.5 -mtick1=4 -tick2d=0.5 -mtick2=4 -tick3d=0.5 -mtick3=4 -lloc=right -lmtick=9 -legend=y -color=rainbowcmyk -slice1=0.9 -slice3=2.6 -slice2=2.5 "

pcolor = " -cmin=2200 -cmax=4000 -unit='P-wave Velocity (m/s)' "

system "mkdir -p ./plot"

system "x_showslice -in=model/vp.bin #{opts} -out=./plot/vp_gt.pdf #{pcolor}&"
system "x_showslice -in=model/vp_init.bin #{opts} -out=./plot/vp_init.pdf #{pcolor}&"

for m in ['waveform', 'adaptive', 'dtw']

    system "x_showslice -in=test_#{m}/iteration_140/model/updated_vp.bin #{opts} -out=./plot/vp_#{m}.pdf #{pcolor} &"

end
