

# Elevation
system "x_showgraph -in=ftopo.txt -ftype=ascii -ptype=2 -size1=6 -size2=2 -label1='Horizontal Distance (m)' -label2='Elevation (m)' -x1beg=0 -x1end=3000 -x2beg=-10 -x2end=310 -tick2beg=0 -linewidth=2 -tick1d=500 -mtick1=4 -tick2d=50 -mtick2=4 -out=plot/elevation.pdf -grid2=y  &"


# Inversion results

opts = "-n1=101 -label1='Depth (km)' -size1=2 -size2=6 -label2='Horizontal Position (km)' -d1=0.01 -d2=0.01 -tick1d=0.25 -mtick1=4 -tick2d=0.5 -mtick2=4 -mask=model/mask.bin -lloc=bottom -lmtick=9 -legend=y -color=rainbowcmyk "

pcolor = " -cmin=2800 -cmax=4500 -unit='P-wave Velocity (m/s)' "
scolor = " -cmin=1750 -cmax=2400 -unit='S-wave Velocity (m/s)' "
rhocolor = " -cmin=2150 -cmax=2410 -unit='Density (kg/m$^3$)' "

system "x_showmatrix -in=model/vp.bin #{opts} -out=./plot/vp_gt.pdf #{pcolor}&"
system "x_showmatrix -in=model/vp_init.bin #{opts} -out=./plot/vp_init.pdf #{pcolor}&"

system "x_showmatrix -in=model/vs.bin #{opts} -out=./plot/vs_gt.pdf #{scolor}&"
system "x_showmatrix -in=model/vs_init.bin #{opts} -out=./plot/vs_init.pdf #{scolor}&"

system "x_showmatrix -in=model/rho.bin #{opts} -out=./plot/rho.pdf #{rhocolor}&"

# waveform
system "x_showmatrix -in=test_waveform/iteration_50/model/updated_vp.bin #{opts} -out=./plot/vp_l2.pdf #{pcolor}&"
system "x_showmatrix -in=test_waveform/iteration_50/model/updated_vs.bin #{opts} -out=./plot/vs_l2.pdf #{scolor}&"

# adaptive
system "x_showmatrix -in=test_adaptive/iteration_50/model/updated_vp.bin #{opts} -out=./plot/vp_awi.pdf #{pcolor}&"
system "x_showmatrix -in=test_adaptive/iteration_50/model/updated_vs.bin #{opts} -out=./plot/vs_awi.pdf #{scolor}&"

# dtw
system "x_showmatrix -in=test_dtw_vector/iteration_50/model/updated_vp.bin #{opts} -out=./plot/vp_dtw.pdf #{pcolor}&"
system "x_showmatrix -in=test_dtw_vector/iteration_50/model/updated_vs.bin #{opts} -out=./plot/vs_dtw.pdf #{scolor}&"


# Adjoint source of awi for illustration

dataopts = " -n1=2501 -size1=5 -size2=5 -lloc=bottom -lmtick=9 -legend=y -n1=2501 -d1=1.0e-3 -label1='Time (s)' -label2='Trace Number' -color=binary -tick2d=50 -mtick2=4 -tick1d=0.5 -mtick1=4 -lwidth=4.5 -unit='Adjoint Source Amplitude' "


system "x_suhdrstrip <data/shot_4_seismogram_z.su >data.bin "
system "x_showmatrix -in=data.bin #{dataopts} -out=./plot/data.pdf -unit='Amplitude ($\\times 10^6$)' -cscale=1.0e6 -clip=1 & "

system "x_suhdrstrip <test_adaptive/iteration_50/adjoint_source/shot_4_seismogram_z.su >ad.bin "
system "x_showmatrix -in=ad.bin -out=./plot/adjsrc_awi.pdf #{dataopts} -cscale=0.1 -unit='Adjoint Source Amplitude ($\\times 10^{1}$)' -clip=1 & "

system "x_suhdrstrip <test_dtw_vector/iteration_50/adjoint_source/shot_4_seismogram_z.su >ad2.bin "
system "x_showmatrix -in=ad2.bin -out=./plot/adjsrc_dtw.pdf #{dataopts} -unit='Adjoint Source Amplitude ($\\times 10^7$)' -cscale=1.0e7 -clip=1 & "

system "x_suhdrstrip <test_waveform/iteration_50/adjoint_source/shot_4_seismogram_z.su >ad3.bin "
system "x_showmatrix -in=ad3.bin -out=./plot/adjsrc_l2.pdf #{dataopts} -unit='Adjoint Source Amplitude ($\\times 10^7$)' -cscale=1.0e7 -clip=1 & "


# Waveform comparison

wiggleopts = " -wigglecolor=b,r -every=6 -size1=5 -size2=5 -label1='Time (s)' -label2='Trace Number' -clip=0.5e-6 -plotlabel='Observed':'Synthetic' "

system "x_showwiggle -in=data/shot_4_seismogram_z.su,test_adaptive/iteration_0/synthetic/shot_4_seismogram_z.su #{wiggleopts} -out=./plot/compare_init.pdf &"

system "x_showwiggle -in=data/shot_4_seismogram_z.su,test_waveform/iteration_50/synthetic/shot_4_seismogram_z.su #{wiggleopts} -out=./plot/compare_l2.pdf &"

system "x_showwiggle -in=data/shot_4_seismogram_z.su,test_dtw_vector/iteration_50/synthetic/shot_4_seismogram_z.su #{wiggleopts} -out=./plot/compare_dtw.pdf &"

system "x_showwiggle -in=data/shot_4_seismogram_z.su,test_adaptive/iteration_50/synthetic/shot_4_seismogram_z.su #{wiggleopts} -out=./plot/compare_awi.pdf &"




