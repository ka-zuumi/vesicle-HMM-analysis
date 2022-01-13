set term pngcairo size 1800,1200 background rgb 'gray'
set output "polarPlot1.png"
set multiplot layout 1,2
set polar
set zeroaxis
set trange [0:2*pi]
set rrange [0:pi]
set size square
set tics front
unset border
unset key
set label 1 'South Pole' font ',18' at screen 0.22,0.1
set label 2 'North Pole' font ',18' at screen 0.72,0.1
set label 3 'Red: POPC Blue: POPE Green: POPS Yellow: CHOL' font ',18' at screen 0.35,0.02

set palette defined ( 0 'red', 1 'blue', 2 'green', 3 'yellow')
cmx = 0.0
cmy = 0.0
cmz = 0.0
R(x,y,z) = atan2(z-cmz,sqrt(((x-cmx)**2)+((y-cmy)**2)))
theta(x,y,z) = pi + atan2(y-cmy,x-cmx)
unset colorbox
set rtics ('0' 0, 'π/4' pi/4, 'π/2' pi/2, '3π/4' 3*pi/4, 'π' pi)
set rtics offset 1,0
unset xtics
unset ytics
piby2=0.5*pi
piby4=0.5*piby2
f(x)=piby4

plot f(t) lw 1 lc rgb "black", 2*f(t) lw 2 lc rgb "black", 3*f(t) lw 1 lc rgb "black", 4*f(t) lw 2 lc rgb "black", "example_frame-hiddenstates.xyz" u (theta($2,$3,$4)):(($5)>0?R($2,$3,$4)+piby2:1/0):1 w p lc palette ps 1 pt 7

set rtics ('0' pi, 'π/4' 3*pi/4, 'π/2' pi/2, '3π/4' pi/4, 'π' 0)
set rtics offset 1,0

plot f(t) lw 1 lc rgb "black", 2*f(t) lw 2 lc rgb "black", 3*f(t) lw 1 lc rgb "black", 4*f(t) lw 2 lc rgb "black", "example_frame-hiddenstates.xyz" u (theta($2,$3,$4)):(($5)>0?piby2-R($2,$3,$4):1/0):1 w p lc palette ps 1 pt 7

unset multiplot
