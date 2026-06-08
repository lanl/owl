
system "x_runf90 create_test.f90"


for v in ['high', 'low']
    for w in ['high', 'low']

        system "cp -rp param_modeling.rb tmp.rb"
        system "echo 'file_vp = model/vp_#{v}.bin' >> tmp.rb "
        system "echo 'file_vs = model/vs_#{w}.bin' >> tmp.rb "
        system "echo 'dir_synthetic = data_#{v}_#{w}' >> tmp.rb "
        system "echo 'dir_snapshot = snapshot_#{v}_#{w}' >> tmp.rb "

        # Data generation
        system "mpirun -np 1 owl_modeling2 tmp.rb"

        # FWI
        for m in ['waveform', 'envelope', 'phase', 'adaptive', 'local-adaptive', 'dtw_1', 'dtw_2']

            system "cat param_fwi.rb ./param_fwi/param_#{m}.rb > tmp.rb"
            system "echo 'dir_record = data_#{v}_#{w}' >> tmp.rb "
            system "echo 'dir_working = test_#{v}_#{w}_#{m}' >> tmp.rb "

            system "mpirun -np 1 owl_fwi2 tmp.rb"

        end

    end
end

# Plotting results
for i in ['low', 'high']
    for j in ['low', 'high']

        system "x_suheadstrip <data_#{i}_#{j}/shot_1_seismogram_x.su >g#{i}#{j}.bin; \
                x_suheadstrip <test_#{i}_#{j}_waveform/iteration_0/synthetic/shot_1_seismogram_x.su >l#{i}#{j}.bin; \
                x_showwiggle -in=g#{i}#{j}.bin,l#{i}#{j}.bin -n1=3 -transp=y -plotlabel='Ground Truth':'Initial Model' \
                -wigglecolor=b,r -d2=1.0e-3 -clip=0.15 -tick2d=0.5 -mtick2=4 -label2='Time (s)' \
                -label1='Trace' -tick1d=1 -o1=1 -out=data_#{i}_#{j}_x.pdf -plotlabelloc=upper_left -along=2 -size1=2 -size2=5 "

        system "x_suheadstrip <data_#{i}_#{j}/shot_1_seismogram_z.su >g#{i}#{j}.bin; \
                x_suheadstrip <test_#{i}_#{j}_waveform/iteration_0/synthetic/shot_1_seismogram_z.su >l#{i}#{j}.bin; \
                x_showwiggle -in=g#{i}#{j}.bin,l#{i}#{j}.bin -n1=3 -transp=y -plotlabel='Ground Truth':'Initial Model' \
                -wigglecolor=b,r -d2=1.0e-3 -clip=0.15 -tick2d=0.5 -mtick2=4 -label2='Time (s)' \
                -label1='Trace' -tick1d=1 -o1=1 -out=data_#{i}_#{j}.pdf -plotlabelloc=upper_left -along=2 -size1=2 -size2=5 "

    end
end

# Plot gradients
system "ln -s ../../../misc/python ./"
system "python plot_result.py"
