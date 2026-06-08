

system "make"

f0 = [12.5]
sz = [0]
refine = [5.0]
c = [0.2e-4]

for f in f0
    h = 0
    for z in sz

        # parameter
        system "./create_test f0=#{f} sz=#{z} sx=300.0 rz=0  "

        ## specfem -- sz should be max(topo) + sz
        #system "mpirun -np 6 ./exec f0=#{f} sz=#{z} dir=./data_specfem_f0=#{f}_sz=#{z} sx=300.0 src_on_surface=y rec_on_surface=y "

        for r in refine

            # owl
            for m in ['fsg']

                system "cp -rp param_modeling.rb tmp.rb "
                system "echo 'which_medium = elastic-tti ' >> tmp.rb "

                system "echo 'free_surface_dz_refine = #{r} ' >> tmp.rb "
                system "echo 'dir_synthetic = ./data_#{m}_refine=#{r}_f0=#{f}_sz=#{z} ' >> tmp.rb "
                system "echo 'dir_snapshot = ./snapshot_#{m}_refine=#{r}_f0=#{f}_sz=#{z} ' >> tmp.rb "
                system "echo 'free_surface_dz_refine = #{r} ' >> tmp.rb "
                if z == 0
                    system "echo 'measure_source_depth_from_surface = y' >> tmp.rb "
                else
                    system "echo 'measure_source_depth_from_surface = n' >> tmp.rb "
                end
                system "echo 'measure_receiver_depth_from_surface = y' >> tmp.rb "
                system "mpirun -np 1 owl_modeling2 tmp.rb "

                system "./compare_seis \
                    dir1=./data_specfem_f0=#{f}_sz=#{z} dir2=./data_#{m}_refine=#{r}_f0=#{f}_sz=#{z} \
                    label1='specfem_f0=#{f}_sz=#{z}' label2='#{m}_refine=#{r}_f0=#{f}_sz=#{z}' "

                clip = [c[h], c[h]]

                l = 0
                for data in ['x', 'z']

                    system "x_showwiggle -size1=4 -size2=6 -along=2 -clip=#{clip[l]} \
                        -in=#{data}_specfem_f0=#{f}_sz=#{z}.bin,#{data}_#{m}_refine=#{r}_f0=#{f}_sz=#{z}.bin \
                        -n1=300 -transpose=y -every=20 -wigglecolor=b,r -plotlabel='SPECFEM':'OWL' -plotlabelloc=lower_right \
                        -wigglestyle=solid,dashed -d1=10 -tick2d=0.5 -mtick2=4 -wigglewidth=1,1 \
                        -label1='Horizontal Position (m)' -label2='Time (s)' -d2=1.0e-3 \
                        -out=c_#{m}_#{data}_refine=#{r}_f0=#{f}_sz=#{z}.pdf & "

                    system "x_showwiggle -size1=4 -size2=6 -along=2 -clip=#{clip[l]/7.0} \
                        -in=#{data}_specfem_f0=#{f}_sz=#{z}.bin,#{data}_#{m}_refine=#{r}_f0=#{f}_sz=#{z}.bin \
                        -n1=300 -transpose=y -every=20 -wigglecolor=b,r -plotlabel='SPECFEM':'OWL' -plotlabelloc=lower_right \
                        -wigglestyle=solid,dashed -d1=10 -tick2d=0.5 -mtick2=4 -wigglewidth=1,1 \
                        -label1='Horizontal Position (m)' -label2='Time (s)' -d2=1.0e-3 \
                        -out=c_#{m}_#{data}_refine=#{r}_f0=#{f}_sz=#{z}_zoom.pdf & "

                    opts = "-ftype=ascii -ptype=3 \
                        -x1beg=-150 -x1end=3140 -size1=6 -size2=2 -clip=#{clip[l]/4.0} -color=binary -reverse2=y \
                        -label1='Horizontal Position (m)' -label2='Depth (m)' -tick1d=500 -x2beg=-120 -x2end=550 \
                        -tick2beg=-200 -tick2d=100 -mtick1=4 -mtick2=1 \
                        -tick1beg=-500 -tick1d=500 -plotorder=2 -markersizemin=10 -markersizemax=10 \
                        -arrow=-200,0,400,0:400,0,400,2990:400,2990,-200,2990 \
                        -arrowfacecolor=w,w,w -arrowlinestyle=solid,solid,solid -arrowstyle=-,-,- \
                        -arrowwidth=1.5,1.5,1.5 -arroworder=10,10,10 \
                        -textloc=2700,480 -textsize=12 -text='CFS-MPML Region' -textcolor=w -tickbottom=n -ticktop=y -label1loc=top -label1pad=10 "

                    if data == 'x'
                        opts = opts + " -select=1,2,3 "
                    else
                        opts = opts + " -select=1,2,4 "
                    end

                    system "x_showgraph -in=./snapshot_#{m}_refine=#{r}_f0=#{f}_sz=#{z}/shot_1_forward_wavefield_4.txt \
                        #{opts} -clip=#{clip[l]/4.0} -out=wave_#{m}_#{data}_refine=#{r}_f0=#{f}_sz=#{z}_step1.png & "
                    system "x_showgraph -select=1,2,4 -in=./snapshot_#{m}_refine=#{r}_f0=#{f}_sz=#{z}/shot_1_forward_wavefield_10.txt \
                        #{opts} -clip=#{clip[l]/15.0} -out=wave_#{m}_#{data}_refine=#{r}_f0=#{f}_sz=#{z}_step2.png & "

                    l = l + 1

                end

            end

        end

        h = h + 1

    end
end
