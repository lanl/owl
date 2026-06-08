
system "mkdir -p ./plot"

# data
system "cat ./data_processed/shot_1.bin ./data_processed/shot_10.bin ./data_processed/shot_20.bin ./data_processed/shot_30.bin ./data_processed/shot_40.bin ./data_processed/shot_49.bin >data.bin "

system "x_showmatrix -in=data.bin -n1=501 -d1=0.001 -label1='Time (s)' -label2='Common-Shot Gather - Trace' -color=binary -size1=3 -size2=9 -clip=3000 -tick1d=0.1 -mtick1=4 -tick2d=40 -mtick2=4 -out=plot/data_raw.pdf & "

# fatt
mopts = "-n1=21 -size1=1.75 -size2=5.5 -tick1d=5 -tick2d=10 -mtick1=4 -mtick2=9 -cmin=350 -cmax=5500 -d1=1 -d2=1 -color=gist_ncar -label1='Depht (m)' -label2='Horizontal Position (m)' -interp=gaussian -legend=y -lloc=bottom -unit='P-wave Velocity (m/s)' -ld=1000 -lmtick=9 "

system "x_showmatrix -in=./fatt/iteration_0/model/vp.bin #{mopts} -out=./plot/vp_init.pdf & "
system "x_showmatrix -in=./fatt/iteration_100/model/updated_vp.bin #{mopts} -out=./plot/vp_fatt.pdf & "

system "x_showgraph -in=geometry_fwi/wavelet.txt -ftype=ascii -ptype=2 -size1=4 -size2=2.5 -label1='Time (s)' -label2='Amplitude' -out=plot/wavelet.pdf -tick1beg=0 -cmin=-40 -cmax=40 -tick1d=0.05 -mtick1=4 -mtick2=1 &"

system "x_showgraph -in=fatt/data_misfit.txt -ftype=ascii -ptype=2 -size1=4 -size2=2.5 -label1='Iteration' -label2='Normalized Data Misfit' -out=plot/fatt_misfit.pdf -tick1beg=0 -cmin=0 -cmax=1 -tick1d=20 -x1end=100 -mtick1=4 -mtick2=1 -tick2d=0.2 -tick2beg=0 -select=1,3 &"


# fwi
iters = [26, 30, 40]

dopts = " -every=2 -wigglecolor=b,r -label1='Time (s)' -wigglewidth=0.75,0.75 -tick1d=0.1 -mtick1=1 -label2='Trace Number' -tick2d=4 -clip=0.05 -scaling=0.25,1 -size1=2.5 -size2=5 "

mopts = "-n1=61 -size1=1.75 -size2=5.5 -tick1d=5 -tick2d=10 -mtick1=4 -mtick2=9 -cmin=150 -cmax=4000 -d1=0.333333 -d2=0.333334 -color=gist_ncar -label1='Depht (m)' -label2='Horizontal Position (m)' -legend=y -lloc=bottom -unit='S-wave Velocity (m/s)' -ld=1000 -lmtick=9 "

system "x_showmatrix -in=model/vs_fatt.bin #{mopts} -out=./plot/vs_init.pdf & "

for i in [1, 2, 3]

    dir = 'test_s' + i.to_s
    iter = iters[i - 1]

    system "x_showmatrix -in=#{dir}/iteration_#{iter}/model/updated_vs.bin #{mopts} -out=./plot/vs_s#{i}.pdf &"

    for shot in [3, 45]

        if shot == 3
            ll = 'lower_left'
        else
            ll = 'lower_right'
        end

        system "x_showwiggle -in=#{dir}/record_processed/shot_#{shot}_seismogram_z.su,#{dir}/iteration_0/synthetic_processed/shot_#{shot}_seismogram_z.su #{dopts} -plotlabel='Observed':'Synthetic' -plotlabelloc=#{ll} -out=./plot/data_init_s#{i}_shot#{shot}.pdf &"

        system "x_showwiggle -in=#{dir}/record_processed/shot_#{shot}_seismogram_z.su,#{dir}/iteration_#{iter}/synthetic_processed/shot_#{shot}_seismogram_z.su #{dopts} -plotlabel='Observed':'Synthetic' -plotlabelloc=#{ll} -out=./plot/data_s#{i}_shot#{shot}.pdf &"

    end

end


