set term pngcairo size 1800,1200 background rgb 'gray'
set output "polarPlot2.png"
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
set label 4 'Hidden State 1: White' tc rgb 'white' font ',18' at screen 0.35,0.02 offset character 0,4
set label 5 'Hidden State 2: Black' tc rgb 'black' font ',18' at screen 0.35,0.02 offset character 0,2
#set label 6 'XYZ from: /home/kazuumi/rsun_lts/cakang/anton/anton_traj/gromacs-popc-chl/only_heads/' font ',12' at screen 0.02,0.98
#set label 7 'Most probable sequences from: firstfourthousand-average/4lipid-system-mostprobable.out' font ',12' at screen 0.02,0.98 offset character 0,-1.5
set palette defined ( -2 'black', -1 'white', 0 'red', 1 'blue', 2 'green', 3 'yellow')
cmx = replacecmx
cmy = replacecmy
cmz = replacecmz
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

plot f(t) lw 1 lc rgb "black", 2*f(t) lw 2 lc rgb "black", 3*f(t) lw 1 lc rgb "black", 4*f(t) lw 2 lc rgb "black", "tmp" u (theta($2,$3,$4)):(($5)>0?R($2,$3,$4)+piby2:1/0):(-($5)) w p lc palette ps 2.0 pt 7,\
     f(t) lw 1 lc rgb "black", 2*f(t) lw 2 lc rgb "black", 3*f(t) lw 1 lc rgb "black", 4*f(t) lw 2 lc rgb "black", "tmp" u (theta($2,$3,$4)):(($5)>0?R($2,$3,$4)+piby2:1/0):1 w p lc palette ps 1 pt 7

set rtics ('0' pi, 'π/4' 3*pi/4, 'π/2' pi/2, '3π/4' pi/4, 'π' 0)
set rtics offset 1,0
plot f(t) lw 1 lc rgb "black", 2*f(t) lw 2 lc rgb "black", 3*f(t) lw 1 lc rgb "black", 4*f(t) lw 2 lc rgb "black", "tmp" u (theta($2,$3,$4)):(($5)>0?piby2-R($2,$3,$4):1/0):(-($5)) w p lc palette ps 2.0 pt 7,\
     f(t) lw 1 lc rgb "black", 2*f(t) lw 2 lc rgb "black", 3*f(t) lw 1 lc rgb "black", 4*f(t) lw 2 lc rgb "black", "tmp" u (theta($2,$3,$4)):(($5)>0?piby2-R($2,$3,$4):1/0):1 w p lc palette ps 1 pt 7

unset multiplot
