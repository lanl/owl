

opts = "-n1=111 -label1='Depth (km)' -size1=2.22 -size2=4.01 -label2='Horizontal Position (km)' -d1=0.02 -d2=0.02 -tick1d=1 -mtick1=9 -tick2d=2 -mtick2=9 -lloc=bottom -lmtick=9 -legend=y -color=rainbowcmyk "

pcolor = " -cmin=1000 -cmax=4600 -unit='P-wave Velocity (m/s)' "

system "mkdir -p ./plot"

system "x_showmatrix -in=model/vp.bin #{opts} -out=./plot/vp_gt.pdf #{pcolor}&"
system "x_showmatrix -in=model/vp_init.bin #{opts} -out=./plot/vp_init.pdf #{pcolor}&"


l = 0
for m in ['waveform', 'adaptive', 'adaptive_local', 'adaptive_spacetime', 'dtw_amp', 'dtw_phase']

    system "x_showmatrix -in=test_#{m}/iteration_100/model/updated_vp.bin #{opts} -out=./plot/vp_#{m}.pdf #{pcolor} &"

    if m == 'waveform'

        system "x_showwiggle -in=data/shot_18_seismogram_p.su,test_#{m}/iteration_0/synthetic/shot_18_seismogram_p.su -wigglecolor=b,r -every=10 -size1=4 -size2=4 -wigglecolor=b,r -wigglewidth=0.75,0.75 -label1='Time (s)' -label2='Trace Number' -tick1d=2 -mtick1=9 -clip=30 -tick2d=100 -mtick2=9 -out=plot/compare_init.pdf -plotlabel='Observed':'Synthetic'  -x1end=5 -tick1d=1  &"

    end

    system "x_showwiggle -in=data/shot_18_seismogram_p.su,test_#{m}/iteration_100/synthetic/shot_18_seismogram_p.su -wigglecolor=b,r -every=10 -size1=4 -size2=4 -wigglecolor=b,r -wigglewidth=0.75,0.75 -label1='Time (s)' -label2='Trace Number' -tick1d=2 -mtick1=9 -clip=30 -tick2d=100 -mtick2=9 -out=plot/compare_#{m}.pdf -plotlabel='Observed':'Synthetic'  -x1end=5 -tick1d=1  &"

    if l == 0
        system "head -n 101 test_#{m}/data_misfit.txt > misfit.txt "
    else
        system "head -n 101 test_#{m}/data_misfit.txt >> misfit.txt "
    end

    l = l + 1

end


system "x_showgraph -ptype=2 -select=1,3 -ftype=ascii -x1beg=0 -x1end=100 -tick1d=20 -mtick1=9 -size1=4 -size2=3 -label1='Iteration Number' -label2='Relative Data Misfit' -x2beg=0 -x2end=1 -tick2d=0.2 -mtick2=1 -in=misfit.txt -n1=101,101,101,101,101,101 -plotlabel='L2':'AWI':'LAWI':'ST-AWI':'GWI-1':'GWI-2' -linewidth=2,2,2,2,2,2 -linecolor=k,b,green,cyan,yellow,r -out=plot/misfit.pdf &"
