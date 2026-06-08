
opts = "-ftype=ascii -ptype=3 \
    -x1beg=0 -x1end=2990 -size1=6 -size2=2 -color=rainbowcmyk -reverse2=y \
    -label1='Horizontal Position (m)' -label2='Depth (m)' -tick1d=500 -x2beg=-120 -x2end=400 \
    -tick2beg=-200 -tick2d=100 -mtick1=4 -mtick2=1 -lmtick=4 -lwidth=5.5 \
    -tick1beg=-500 -tick1d=500 -plotorder=2 -markersizemin=20 -markersizemax=20 \
    -legend=y -lloc=bottom -tickbottom=n -ticktop=y -label1loc=top -label1pad=10 "

system "x_showgraph -in=model.txt -select=1,2,3 -unit='P-wave Velocity (m/s)' -ld=300 -ltickbeg=2500 #{opts} -out=model_vp.png & "
system "x_showgraph -in=model.txt -select=1,2,4 -unit='S-wave Velocity (m/s)' -ld=250 -ltickbeg=1750 #{opts} -out=model_vs.png & "
system "x_showgraph -in=model.txt -select=1,2,5 -unit='Density (kg/m3)' #{opts} -out=model_rho.png & "
system "x_showgraph -in=model.txt -select=1,2,6 -unit='Thomsen Parameter $\\varepsilon$' #{opts} -out=model_eps.png & "
system "x_showgraph -in=model.txt -select=1,2,7 -unit='Thomsen Parameter $\\delta$' #{opts} -out=model_del.png & "
system "x_showgraph -in=model.txt -select=1,2,8 -unit='TI Tilt Angle $\\theta$ (rad)' #{opts} -out=model_the.png & "
