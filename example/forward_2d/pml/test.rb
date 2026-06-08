

#system "make"
#
#
system "x_runf90 ./create_test.f90 f0=5 sz=1000 sx=1000 rz=500 "
#system "mpirun -np 6 ./exec f0=5 sz=1000 dir=./data_specfem sx=1000 src_on_surface=n rec_on_surface=n rz=500 free_surface=n "

system "mpirun -np 1 owl_modeling2 param_modeling.rb "

abort

#
#system "x_runf90 interp_specfem_solution.f90 "
#

system "x_showgraph -in=energy.txt -ftype=ascii -ptype=2 -norm2=log -size1=5 -size2=3 \
	-label1='Wave Propagation Time (s)' -label2='Kinetic Energy $||\\mathbf{v}||_2^2$' \
	-tick1beg=0 -x2beg=1.0e-20 -x2end=1.0e-5 -tick2d=1000 -tick1d=10 -mtick1=4 \
	-linewidth=1.5 -out=energy.pdf &"

abort


opts = "-n1=330 -size1=3.5 -size2=3.5 -o1=-150 -o2=-150 -d1=10 -d2=10 -label1='Depth (m)' \
    -label2='Horizontal Position (m)' -color=binary \
    -tick1d=500 -tick2d=500 -mtick1=4 -mtick2=4 -tick1beg=-500 -tick2beg=-500 -clip=#{1.0e-5} \
    -arrow=0,0,2990,0:2990,0,2990,2990:2990,2990,0,2990:0,2990,0,0 \
    -arrowfacecolor=w,w,w,w -arrowlinestyle=solid,solid,solid,solid -arrowstyle=-,-,-,- \
    -arrowwidth=1.5,1.5,1.5,1.5 -arroworder=10,10,10,10 \
    -textloc=3080,1500 -textsize=12 -text='CFS-MPML Region' -textcolor=w "

system "x_showmatrix -in=./snapshot_fsg/shot_1_forward_wavefield_z_4.bin \
    #{opts} -out=wave_owl_step1.png & "

system "x_showmatrix -in=./snapshot_fsg/shot_1_forward_wavefield_z_11.bin \
    #{opts} -out=wave_owl_step2.png & "

system "x_showmatrix -in=./wave_specfem_1.bin \
    #{opts} -text='CFS-PML Region' -out=wave_specfem_step1.png & "

system "x_showmatrix -in=./wave_specfem_2.bin \
    #{opts} -text='CFS-PML Region' -textloc=1500,3080 -textrot=90 -out=wave_specfem_step2.png & "
