
system "x_runf90 create_test.f90"

system "ln -s ../homo/data_low ./"
system "ln -s ../homo/data_high ./"

# FWI
for v in ['low', 'high']

    for m in ['waveform', 'envelope', 'phase', 'adaptive', 'local-adaptive', 'dtw_1', 'dtw_2']

        system "cat param_fwi.rb ./param_fwi/param_#{m}.rb > tmp.rb"
        system "echo 'dir_record = data_#{v}' >> tmp.rb "
        system "echo 'dir_working = test_#{v}_#{m}' >> tmp.rb "

        system "mpirun -np 1 owl_fwi2 tmp.rb"

    end

end

# Plotting data
for i in ['low', 'high']

    system "x_suheadstrip <data_#{i}/shot_1_seismogram_p.su >g#{i}.bin; \
            x_suheadstrip <test_#{i}_waveform/iteration_0/synthetic/shot_1_seismogram_p.su >l#{i}.bin; \
            x_showgraph -in=g#{i}.bin,l#{i}.bin -n1=1501,1501 -plotlabel='Ground Truth':'Initial Model' \
            -linecolor=b,r -d1=1.0e-3 -tick2d=25 -x2beg=-75 -x2end=75 -tick1d=0.5 -mtick1=4 -label1='Time (s)' \
            -label2='Amplitude' -out=data_#{i}.pdf -plotlabelloc=lower_right -size1=5 -size2=2 &"

end


# Plot gradients
system "ln -s ../../../misc/python ./"
system "python plot_result.py"
