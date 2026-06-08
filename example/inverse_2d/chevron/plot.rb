
system "mkdir -p ./plot"


system "x_showgraph -in=geometry/chevron_wavelet.txt -ftype=ascii -ptype=2 -x2beg=-35 -x2end=35 -tick2beg=-40 -x1beg=0 -x1end=0.6 -tick1d=0.1 -mtick1=4 -tick2d=20 -mtick2=4 -label1='Time (s)' -label2='Amplitude' -out=plot/wavelet.pdf -size1=4 -size2=3 &"


shot = 1000
iters = [30, 50, 70, 100]
labels = ['a', 'b', 'c', 'd']
clips = [0.2e-6, 0.5e-6, 1e-6, 1.5e-6]


mopts = " -mask=model/mask.bin -n1=480 -tick2beg=5 -x2beg=#{500*0.0125} -x2end=#{3360*0.0125} -size1=2.5 -size2=5 -color=gist_ncar -legend=y -lloc=bottom -unit='P-wave Velocity (m/s)' -cmin=1500 -cmax=4400 -lmtick=4 -ld=500 -d1=0.0125 -d2=0.0125 -label1='Depth (km)' -tick1d=1 -mtick1=9 -label2='Horizontal Position (km)' -tick2d=5 -mtick2=9 "

system "x_showmatrix -in=model/vp.bin #{mopts} -out=./plot/vp_init.pdf &"
system "x_suhdrstrip <field_data/csg/shot_1000_seismogram_p.su >data.bin; \
    x_showmatrix -in=data.bin -n1=2001 -d1=0.004 -clip=8 -label2='Trace Number' -size1=5 -size2=4 -label1='Time (s)' -tick1d=2 -mtick1=4 -tick2d=50 -mtick2=4 -color=binary -legend=y -lloc=right -cscale=1.0e7 -unit='Amplitude ($\\times 10^{-7}$)' -out=plot/data.pdf &"

for l in 1..4

    iter = iters[l - 1]
    dir = 'test_' + labels[l - 1]
    clip = clips[l - 1]

    dopts = "-x1end=6 -every=10 -wigglecolor=b,r -clip=#{clip} -label2='Trace Number' -scaling=1,1 -size1=4 -size2=6 -plotlabel='Observed':'Synthetic' -plotlabelloc=upper_right -wigglewidth=0.75,0.75 -label1='Time (s)' -tick1d=2 -mtick1=4 -tick2d=50 -mtick2=4 "

    system "x_showwiggle -in=#{dir}/record_processed/shot_#{shot}_seismogram_p.su,#{dir}/iteration_0/synthetic_processed/shot_#{shot}_seismogram_p.su #{dopts} -out=./plot/data_init_stage_#{labels[l - 1]}.pdf &"

    system "x_showwiggle -in=#{dir}/record_processed/shot_#{shot}_seismogram_p.su,#{dir}/iteration_#{iter}/synthetic_processed/shot_#{shot}_seismogram_p.su #{dopts} -out=./plot/data_stage_#{labels[l - 1]}.pdf &"

    system "x_showmatrix -in=#{dir}/iteration_#{iter}/model/updated_vp.bin #{mopts} -out=./plot/vp_stage_#{labels[l - 1]}.pdf &"

end
