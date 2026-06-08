
system "make"
system "./create_test"

#make clean
#make
#mpirun -np 40 ./exec
#
#

system "export OMP_NUM_THREADS=1"
system "mpirun -np 42 owl_modeling3 param_modeling.rb"

system "./compare_seis dir1=data_specfem dir2=data_fsg label1=specfem label2=fsg"

system "x_showcontour -in=topo.bin -n1=35 -d1=10 -d2=10 -size1=3 -size2=5 -color=binary -shading_scale=0.25 -tick1d=50 -mtick1=4 -tick2d=500 -mtick2=4 -clabelsize=0 -contourlevel=10 -mcontour=1 -contourwidth=1.5 -mcontourstyle=dashed -legend=y -unit='Elevation (m)' -out=topo.pdf -contourfill=y -color=terrain -reverse1=y -label1='Horizontal Y (m)' -label2='Horizontal X (m)' -curve=geometry/src.txt,geometry/rec_subset.txt -curveselect=2,1 -curvestyle=scatter*,scatterv -curvefacecolor=r,k -curveedgecolor=none,none -curvesize=50,6 & "


for c in ['x', 'y', 'z']

    system "x_showwiggle -size1=6 -size2=4 -along=2 -n1=200 -transpose=y -clip=1.5e-7 -in=$c''_specfem.bin,$c''_fsg.bin \
    -wigglestyle=solid,dashed -plotlabelloc=lower_right -label1='Trace Number' -label2='Time (s)' \
    -d2=0.001 -mtick1=9 -mtick2=4 -tick2d=0.5 \
    -every=5 -wigglecolor=b,r -wigglewidth=0.75,0.75 -plotlabel='SPECFEM':'OWL' -out=waveform_$c.pdf"

end
