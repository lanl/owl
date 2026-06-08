
system "make"

f0 = [12.5]
sz = [0, 1000]
refine = [4.0]
c = [0.25e-4, 0.1e-4]

for f in f0
    h = 0
    for z in sz

        # specfem -- sz should be max(topo) + sz
        if z == 0
            system "./create_test f0=#{f} sz=#{z} sx=1000 "
            #            system "mpirun -np 6 ./exec f0=#{f} sz=#{z} dir=./data_specfem_f0=#{f}_sz=#{z} sx=1000 src_on_surface=y rec_on_surface=y "
        else
            system "./create_test f0=#{f} sz=#{z} sx=1000 rz=500 "
            #            system "mpirun -np 6 ./exec f0=#{f} sz=#{z} dir=./data_specfem_f0=#{f}_sz=#{z} sx=1000 src_on_surface=n rec_on_surface=n rz=500 free_surface=n "
        end

        for r in refine

            # owl
            for m in ['fsg']

                system "cp -rp param_modeling.rb tmp.rb "
                system "echo 'which_medium = elastic-tti ' >> tmp.rb "

                system "echo 'dir_synthetic = ./data_#{m}_refine=#{r}_f0=#{f}_sz=#{z} ' >> tmp.rb "
                system "echo 'dir_snapshot = ./snapshot_#{m}_refine=#{r}_f0=#{f}_sz=#{z} ' >> tmp.rb "
                system "echo 'free_surface_dz_refine = #{r} ' >> tmp.rb "
                if z == 0
                    system "echo 'yn_free_surface = y' >> tmp.rb "
                    system "echo 'measure_source_depth_from_surface = y' >> tmp.rb "
                    system "echo 'measure_receiver_depth_from_surface = y' >> tmp.rb "
                else
                    system "echo 'yn_free_surface = n' >> tmp.rb "
                end
                system "mpirun -np 1 owl_modeling2 tmp.rb "

#                system "./compare_seis \
#                    dir1=./data_specfem_f0=#{f}_sz=#{z} dir2=./data_#{m}_refine=#{r}_f0=#{f}_sz=#{z} \
#                    label1='specfem_f0=#{f}_sz=#{z}' label2='#{m}_refine=#{r}_f0=#{f}_sz=#{z}' "
#
#                clip = [c[h], c[h]]
#
#                l = 0
#                for data in ['x', 'z']
#
#                    system "x_showwiggle -size1=4 -size2=6 -along=2 -clip=#{clip[l]*4} \
#                        -in=#{data}_specfem_f0=#{f}_sz=#{z}.bin,#{data}_#{m}_refine=#{r}_f0=#{f}_sz=#{z}.bin \
#                        -n1=300 -transpose=y -every=20 -wigglecolor=b,r -plotlabel='SPECFEM':'OWL' -plotlabelloc=upper_right \
#                        -wigglestyle=solid,dashed -d1=10 -tick2d=0.5 -mtick2=4 -wigglewidth=1,1 \
#                        -label1='Horizontal Position (m)' -label2='Time (s)' -d2=1.0e-3 \
#                        -out=c_#{m}_#{data}_refine=#{r}_f0=#{f}_sz=#{z}.pdf & "
#
#                    system "x_showwiggle -size1=4 -size2=6 -along=2 -clip=#{clip[l]/2.0} \
#                        -in=#{data}_specfem_f0=#{f}_sz=#{z}.bin,#{data}_#{m}_refine=#{r}_f0=#{f}_sz=#{z}.bin \
#                        -n1=300 -transpose=y -every=20 -wigglecolor=b,r -plotlabel='SPECFEM':'OWL' -plotlabelloc=upper_right \
#                        -wigglestyle=solid,dashed -d1=10 -tick2d=0.5 -mtick2=4 -wigglewidth=1,1 \
#                        -label1='Horizontal Position (m)' -label2='Time (s)' -d2=1.0e-3 \
#                        -out=c_#{m}_#{data}_refine=#{r}_f0=#{f}_sz=#{z}_zoom.pdf & "
#
#                    if z == 0
#
#                        opts = "-ftype=ascii -ptype=3 \
#                            -x1beg=-150 -x1end=3140 -size1=#{3*3.3/2.15} -size2=3 -clip=#{clip[l]/3.0} -color=binary -reverse2=y \
#                            -label1='Horizontal Position (m)' -label2='Depth (m)' -tick1d=500 -x2beg=0 -x2end=2150 \
#                            -tick2beg=-500 -tick2d=500 -mtick1=4 -mtick2=4 \
#                            -tick1beg=-500 -tick1d=500 -plotorder=2 -markersizemin=5 -markersizemax=5 \
#                            -arrow=0,0,1990,0:1990,0,1990,2990:1990,2990,0,2990 \
#                            -arrowfacecolor=w,w,w -arrowlinestyle=solid,solid,solid -arrowstyle=-,-,- \
#                            -arrowwidth=1.5,1.5,1.5 -arroworder=10,10,10 \
#                            -textloc=1500,2080 -textsize=12 -text='CFS-MPML Region' -textcolor=w \
#                            -ticktop=y -label1loc=top -tickbottom=n -label1pad=10 "
#
#                        system "x_showgraph -in=./snapshot_#{m}_refine=#{r}_f0=#{f}_sz=#{z}/shot_1_forward_wavefield_5.txt \
#                            #{opts} -out=wave_#{m}_#{data}_refine=#{r}_f0=#{f}_sz=#{z}_step1.png & "
#                        system "x_showgraph -in=./snapshot_#{m}_refine=#{r}_f0=#{f}_sz=#{z}/shot_1_forward_wavefield_8.txt \
#                            #{opts} -out=wave_#{m}_#{data}_refine=#{r}_f0=#{f}_sz=#{z}_step2.png & "
#
#                    else
#
#                        opts = "-n1=230 -size1=3 -size2=#{3*3.3/2.3} -o1=-150 -o2=-150 -d1=10 -d2=10 -label1='Depth (m)' \
#                            -label2='Horizontal Position (m)' -color=binary \
#                            -tick1d=500 -tick2d=500 -mtick1=4 -mtick2=4 -tick1beg=-500 -tick2beg=-500 -clip=#{clip[l]/2.0} \
#                            -arrow=0,0,1990,0:1990,0,1990,2990:1990,2990,0,2990:0,2990,0,0 \
#                            -arrowfacecolor=w,w,w,w -arrowlinestyle=solid,solid,solid,solid -arrowstyle=-,-,-,- \
#                            -arrowwidth=1.5,1.5,1.5,1.5 -arroworder=10,10,10,10 \
#                            -textloc=2080,1500 -textsize=12 -text='CFS-MPML Region' -textcolor=w "
#
#                        system "x_showmatrix -in=./snapshot_#{m}_refine=#{r}_f0=#{f}_sz=#{z}/shot_1_forward_wavefield_z_5.bin \
#                            #{opts} -out=wave_#{m}_#{data}_refine=#{r}_f0=#{f}_sz=#{z}_step1.png & "
#                        system "x_showmatrix -in=./snapshot_#{m}_refine=#{r}_f0=#{f}_sz=#{z}/shot_1_forward_wavefield_z_8.bin \
#                            #{opts} -clip=#{clip[l]/4.0} -out=wave_#{m}_#{data}_refine=#{r}_f0=#{f}_sz=#{z}_step2.png & "
#
#                    end
#
#                    l = l + 1
#
#                end

            end

        end

        h = h + 1

    end
end
