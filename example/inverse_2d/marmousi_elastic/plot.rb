

opts = "-n1=111 -label1='Depth (km)' -size1=2.22 -size2=4.01 -label2='Horizontal Position (km)' -d1=0.02 -d2=0.02 -tick1d=1 -mtick1=9 -tick2d=2 -mtick2=9 -lloc=bottom -lmtick=9 -legend=y -color=rainbowcmyk "

pcolor = " -cmin=1000 -cmax=4600 -unit='P-wave Velocity (m/s)' "
scolor = " -cmin=680 -cmax=2700 -unit='S-wave Velocity (m/s)' "

system "mkdir -p ./plot"

system "x_showmatrix -in=model/vp.bin #{opts} -out=./plot/vp_gt.pdf #{pcolor}&"
system "x_showmatrix -in=model/vp_init.bin #{opts} -out=./plot/vp_init.pdf #{pcolor}&"


system "x_showmatrix -in=model/vs.bin #{opts} -out=./plot/vs_gt.pdf #{scolor}&"
system "x_showmatrix -in=model/vs_init.bin #{opts} -out=./plot/vs_init.pdf #{scolor}&"


wopts = " -wigglecolor=b,r -every=10 -size1=4 -size2=4 -wigglecolor=b,r -wigglewidth=0.75,0.75 -label1='Time (s)' -label2='Trace Number' -tick1d=2 -mtick1=9 -clip=8e-8 -tick2d=100 -mtick2=9 -plotlabel='Observed':'Synthetic'  -x1end=6 -tick1d=1  "

l = 0
for m in ['waveform', 'adaptive', 'dtw_vector']

    system "x_showmatrix -in=test_#{m}/iteration_100/model/updated_vp.bin #{opts} -out=./plot/vp_#{m}.pdf #{pcolor} &"
    system "x_showmatrix -in=test_#{m}/iteration_100/model/updated_vs.bin #{opts} -out=./plot/vs_#{m}.pdf #{scolor} &"

    if m == 'waveform'
        system "x_showwiggle #{wopts} -in=data/shot_18_seismogram_z.su,test_#{m}/iteration_0/synthetic/shot_18_seismogram_z.su  -out=plot/compare_init.pdf &"
    end

    system "x_showwiggle #{wopts} -in=data/shot_18_seismogram_z.su,test_#{m}/iteration_100/synthetic/shot_18_seismogram_z.su  -out=plot/compare_#{m}.pdf &"

    if l == 0
        system "head -n 101 test_#{m}/data_misfit.txt > misfit.txt "
    else
        system "head -n 101 test_#{m}/data_misfit.txt >> misfit.txt "
    end

    l = l + 1

end


system "x_showgraph -ptype=2 -select=1,3 -ftype=ascii -x1beg=0 -x1end=100 -tick1d=20 -mtick1=9 -size1=4 -size2=3 -label1='Iteration Number' -label2='Relative Data Misfit' -x2beg=0 -x2end=1 -tick2d=0.2 -mtick2=1 -in=misfit.txt -n1=101,101,101 -plotlabel='L2':'AWI':'GWI-2' -linewidth=2,2,2,2,2,2 -linecolor=k,b,r -out=plot/misfit.pdf &"


#system "x_showgraph -ptype=2 -select=1,3 -ftype=ascii -x1beg=0 -x1end=100 -tick1d=20 -mtick1=9 -size1=4 -size2=3 -label1='Iteration Number' -label2='Relative Data Misfit' -x2beg=0 -x2end=1 -tick2d=0.2 -mtick2=1 -in=test_rwi/data_misfit.txt -n1=101 -plotlabel='RWI' -linewidth=2,2,2,2,2,2 -linecolor=b -out=plot/misfit_rwi.pdf &"



abort


#
#system "x_showgraph -in=ftopo.txt -ftype=ascii -ptype=2 -size1=6 -size2=2 -label1='Horizontal Distance (m)' -label2='Elevation (m)' -x1beg=0 -x1end=3000 -x2beg=-10 -x2end=310 -tick2beg=0 -linewidth=2 -tick1d=500 -mtick1=4 -tick2d=50 -mtick2=4 -out=plot/elevation.pdf -grid2=y  &"
#
#abort




#
#
## adjoint source of awi for illustration
#
#dataopts = " -n1=2501 -size1=5 -size2=5 -lloc=bottom -lmtick=9 -legend=y -n1=2501 -d1=1.0e-3 -label1='Time (s)' -label2='Trace Number' -color=binary -tick2d=50 -mtick2=4 -tick1d=0.5 -mtick1=4 -lwidth=4.5 -unit='Adjoint Source Amplitude' "
#
#
#system "x_suhdrstrip <data/shot_4_seismogram_z.su >data.bin "
#system "x_showmatrix -in=data.bin #{dataopts} -out=./plot/data.pdf -unit='Amplitude ($\\times 10^6$)' -cscale=1.0e6 -clip=1 & "
#
#system "x_suhdrstrip <test_adaptive/iteration_50/adjoint_source/shot_4_seismogram_z.su >ad.bin "
#system "x_showmatrix -in=ad.bin -out=./plot/adjsrc_awi.pdf #{dataopts} -cscale=0.1 -unit='Adjoint Source Amplitude ($\\times 10^{1}$)' -clip=1 & "
#
#system "x_suhdrstrip <test_dtw_vector/iteration_50/adjoint_source/shot_4_seismogram_z.su >ad2.bin "
#system "x_showmatrix -in=ad2.bin -out=./plot/adjsrc_dtw.pdf #{dataopts} -unit='Adjoint Source Amplitude ($\\times 10^7$)' -cscale=1.0e7 -clip=1 & "
#
#system "x_suhdrstrip <test_waveform/iteration_50/adjoint_source/shot_4_seismogram_z.su >ad3.bin "
#system "x_showmatrix -in=ad3.bin -out=./plot/adjsrc_l2.pdf #{dataopts} -unit='Adjoint Source Amplitude ($\\times 10^7$)' -cscale=1.0e7 -clip=1 & "
#
#
#wiggleopts = " -wigglecolor=b,r -every=6 -size1=5 -size2=5 -label1='Time (s)' -label2='Trace Number' -clip=0.5e-6 "
#
#system "x_showwiggle -in=data/shot_4_seismogram_z.su,test_adaptive/iteration_0/synthetic/shot_4_seismogram_z.su #{wiggleopts} -out=./plot/compare_init.pdf &"
#
#system "x_showwiggle -in=data/shot_4_seismogram_z.su,test_waveform/iteration_50/synthetic/shot_4_seismogram_z.su #{wiggleopts} -out=./plot/compare_l2.pdf &"
#
#system "x_showwiggle -in=data/shot_4_seismogram_z.su,test_dtw_vector/iteration_50/synthetic/shot_4_seismogram_z.su #{wiggleopts} -out=./plot/compare_dtw.pdf &"
#
#system "x_showwiggle -in=data/shot_4_seismogram_z.su,test_adaptive/iteration_50/synthetic/shot_4_seismogram_z.su #{wiggleopts} -out=./plot/compare_awi.pdf &"
#
#
#
#
