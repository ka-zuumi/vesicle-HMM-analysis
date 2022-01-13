#!/bin/bash

#
# Example Usage:
# ./makeDeltaTetrahedronGraph.sh <(grep -v '^#' 4lipid-system-timeaverage.dat) <(grep -v '^#' hmm-model1-distribution.dat) hmm-model1-hiddenstates.png "HMM 2-State Comparison, Outer Membrane"
#

tmpfile1="tmp1"
tmpfile2="tmp2"

############################################################################

if [ $# -gt "4" ] || [ $# -lt 3 ]; then
  echo ""
  echo "Error: Requires exactly three or four arguments:"
  echo "       1) 5+ column data file with quaternary values"
  echo "       2) 4+2 column data file with quaternary values"
  echo "       3) Name of output image file (must be .PNG)"
  echo "May have an optional third argument:"
  echo "       4) Text to put on image (to identify the graph)"
  echo "STOPPING"
  exit
fi

if [ $# -eq "4" ]; then
  gnuplotline=$4
else
  gnuplotline=""
fi

imagename=$3
awk '{diff=sqrt((($5)-($6))**2); print $0, diff}' $2 > $tmpfile1

############################################################################

# Get the maximum probability to have a sense of how to scale
# the point sizes
minP=$(awk 'BEGIN {min=1} {diff=($5)-($6); if (diff < min) min = diff} END {print 100*min}' $tmpfile1)
maxP=$(awk 'BEGIN {max=-1} {diff=($5)-($6); if (diff > max) max = diff} END {print 100*max}' $tmpfile1)
maxABSP=$(awk 'BEGIN {max=-1} {if ($7 > max) max = $7} END {print max}' $tmpfile1)

# Now get the distribution-averaged measurements
paste -d' ' <(awk '{$7=""; print $0}' $tmpfile1) <(awk '{$1=""; $2=""; $3=""; $4=""; print $0}' $1) | awk '{for (i=8;i<=NF;i++) {sum[1,i]+=($5)*($i); sum[2,i]+=($6)*($i); sum[3,i]+=($7)*($i)}} END {for (i=8;i<=NF;i++) {print sum[1,i], sum[2,i], sum[3,i]}}' > $tmpfile2

# Let's visualize some properties about the
# local lipid composition

gnuplot <<- EOF
set terminal pngcairo size 2400,1600
set output '$imagename'
#set size ratio 0.866
#set bmargin 3
#set lmargin 3
#set rmargin 3
#set tmargin 3

set multiplot

set yrange [0:0.966]
set xrange [0:1]
set noborder
set noxtics
set noytics

set label 5 "$gnuplotline " font ",24" at screen 0.02,0.98

set label 'Measurement'    font ",24" left at screen 0.52,0.15 offset character 0,0
set label 'Hidden State 1' font ",24" left at screen 0.68,0.15 offset character 0,0 tc rgb "blue"
set label 'Hidden State 2' font ",24" left at screen 0.80,0.15 offset character 0,0 tc rgb "red"
set label 'Total'          font ",24" left at screen 0.92,0.15 offset character 0,0

set label 'Avg Surface Area (A^2)'                             font ",24" left at screen 0.52,0.15 offset character 0,-2
set label '$(head -1 $tmpfile2 | tail -1 | awk "{print \$1}")' font ",24" left at screen 0.68,0.15 offset character 0,-2 tc rgb "blue"
set label '$(head -1 $tmpfile2 | tail -1 | awk "{print \$2}")' font ",24" left at screen 0.80,0.15 offset character 0,-2 tc rgb "red"
set label '$(head -1 $tmpfile2 | tail -1 | awk "{print \$3}")' font ",24" left at screen 0.92,0.15 offset character 0,-2

set label 'CHOL Amount (%)'                                        font ",24" left at screen 0.52,0.15 offset character 0,-4
set label '$(head -2 $tmpfile2 | tail -1 | awk "{print 100*\$1}")' font ",24" left at screen 0.68,0.15 offset character 0,-4 tc rgb "blue"
set label '$(head -2 $tmpfile2 | tail -1 | awk "{print 100*\$2}")' font ",24" left at screen 0.80,0.15 offset character 0,-4 tc rgb "red"
set label '$(head -2 $tmpfile2 | tail -1 | awk "{print 100*\$3}")' font ",24" left at screen 0.92,0.15 offset character 0,-4

set label 'POPC Amount (%)'                                        font ",24" left at screen 0.52,0.15 offset character 0,-6
set label '$(head -3 $tmpfile2 | tail -1 | awk "{print 100*\$1}")' font ",24" left at screen 0.68,0.15 offset character 0,-6 tc rgb "blue"
set label '$(head -3 $tmpfile2 | tail -1 | awk "{print 100*\$2}")' font ",24" left at screen 0.80,0.15 offset character 0,-6 tc rgb "red"
set label '$(head -3 $tmpfile2 | tail -1 | awk "{print 100*\$3}")' font ",24" left at screen 0.92,0.15 offset character 0,-6

set label 'POPE Amount (%)'                                        font ",24" left at screen 0.52,0.15 offset character 0,-8
set label '$(head -4 $tmpfile2 | tail -1 | awk "{print 100*\$1}")' font ",24" left at screen 0.68,0.15 offset character 0,-8 tc rgb "blue"
set label '$(head -4 $tmpfile2 | tail -1 | awk "{print 100*\$2}")' font ",24" left at screen 0.80,0.15 offset character 0,-8 tc rgb "red"
set label '$(head -4 $tmpfile2 | tail -1 | awk "{print 100*\$3}")' font ",24" left at screen 0.92,0.15 offset character 0,-8

set label 'POPS Amount (%)'                                        font ",24" left at screen 0.52,0.15 offset character 0,-10
set label '$(head -5 $tmpfile2 | tail -1 | awk "{print 100*\$1}")' font ",24" left at screen 0.68,0.15 offset character 0,-10 tc rgb "blue"
set label '$(head -5 $tmpfile2 | tail -1 | awk "{print 100*\$2}")' font ",24" left at screen 0.80,0.15 offset character 0,-10 tc rgb "red"
set label '$(head -5 $tmpfile2 | tail -1 | awk "{print 100*\$3}")' font ",24" left at screen 0.92,0.15 offset character 0,-10

 
set style line 1 lt 1 lw 3 pt -1 ps 1 lc rgb "black"
set style line 2 lt 5 lw 1 pt -1 ps 1 lc rgb "cyan"
set style line 3 lt 5 dt 2 lw 3 pt -1 ps 1 lc rgb "gray"

#######################################################################################################################

unset key
set palette defined ($minP "red", 0 "#ffbbbb", 0 "#bbffff", $maxP "blue")
set cblabel "Δ Probability (%)" font ",24" offset character 2.5,0
set cbrange [$minP:$maxP]
set cbtics font ",20"
set colorbox user origin screen 0.01,0.1 size screen 0.015,0.60

set bmargin at screen 0.05
set lmargin at screen 0.05
set rmargin at screen 0.35
set tmargin at screen 0.48

# x
set arrow 1 from 0,0 to 1, 0.0 nohead linestyle 1
set arrow 2 from 1.0/6.0, 0.0 to 1.0-0.5*5.0/6.0, 0.5*sqrt(3)*5.0/6.0 nohead linestyle 2
set arrow 3 from 2.0/6.0, 0.0 to 1.0-0.5*4.0/6.0, 0.5*sqrt(3)*4.0/6.0 nohead linestyle 2
set arrow 4 from 3.0/6.0, 0.0 to 1.0-0.5*3.0/6.0, 0.5*sqrt(3)*3.0/6.0 nohead linestyle 2
set arrow 5 from 4.0/6.0, 0.0 to 1.0-0.5*2.0/6.0, 0.5*sqrt(3)*2.0/6.0 nohead linestyle 2
set arrow 6 from 5.0/6.0, 0.0 to 1.0-0.5*1.0/6.0, 0.5*sqrt(3)*1.0/6.0 nohead linestyle 2

# z
set arrow 7 from 1, 0 to 0.50, 0.866 nohead linestyle 1
set arrow 8 from 1.0/6.0, 0.0 to 0.5*1.0/6.0, 0.5*sqrt(3)*1.0/6.0 nohead linestyle 2
set arrow 9 from 2.0/6.0, 0.0 to 0.5*2.0/6.0, 0.5*sqrt(3)*2.0/6.0 nohead linestyle 2
set arrow 10 from 3.0/6.0, 0.0 to 0.5*3.0/6.0, 0.5*sqrt(3)*3.0/6.0 nohead linestyle 2
set arrow 11 from 4.0/6.0, 0.0 to 0.5*4.0/6.0, 0.5*sqrt(3)*4.0/6.0 nohead linestyle 2
set arrow 12 from 5.0/6.0, 0.0 to 0.5*5.0/6.0, 0.5*sqrt(3)*5.0/6.0 nohead linestyle 2

# y
set arrow 13 from 0.50, 0.866 to 0,0 nohead linestyle 1
set arrow 14 from 0.5*1.0/6.0, 0.5*sqrt(3)*1.0/6.0 to 1.0-0.5*1.0/6.0, 0.5*sqrt(3)*1.0/6.0 nohead linestyle 2
set arrow 15 from 0.5*2.0/6.0, 0.5*sqrt(3)*2.0/6.0 to 1.0-0.5*2.0/6.0, 0.5*sqrt(3)*2.0/6.0 nohead linestyle 2
set arrow 16 from 0.5*3.0/6.0, 0.5*sqrt(3)*3.0/6.0 to 1.0-0.5*3.0/6.0, 0.5*sqrt(3)*3.0/6.0 nohead linestyle 2
set arrow 17 from 0.5*4.0/6.0, 0.5*sqrt(3)*4.0/6.0 to 1.0-0.5*4.0/6.0, 0.5*sqrt(3)*4.0/6.0 nohead linestyle 2
set arrow 18 from 0.5*5.0/6.0, 0.5*sqrt(3)*5.0/6.0 to 1.0-0.5*5.0/6.0, 0.5*sqrt(3)*5.0/6.0 nohead linestyle 2

plot "<awk '{if (\$1 == 0) {print (\$2+2*\$4)/(2*(\$2+\$3+\$4)), sqrt(3)*\$2/(2*(\$2+\$3+\$4)), 100*((\$5)-(\$6)), 1+12*(\$7)/$maxABSP}}' $tmpfile1" u 1:2:4:3 w p pt 7 ps variable lt palette

unset label 1
unset label 2
unset label 3

#######################################################################################################################

unset arrow

set bmargin at screen 0.26
set lmargin at screen 0.28
set rmargin at screen 0.50
set tmargin at screen 0.56

# x
set arrow 1 from 0,0 to 1, 0.0 nohead linestyle 1
set arrow 2 from 1.0/5.0, 0.0 to 1.0-0.5*4.0/5.0, 0.5*sqrt(3)*4.0/5.0 nohead linestyle 2
set arrow 3 from 2.0/5.0, 0.0 to 1.0-0.5*3.0/5.0, 0.5*sqrt(3)*3.0/5.0 nohead linestyle 2
set arrow 4 from 3.0/5.0, 0.0 to 1.0-0.5*2.0/5.0, 0.5*sqrt(3)*2.0/5.0 nohead linestyle 2
set arrow 5 from 4.0/5.0, 0.0 to 1.0-0.5*1.0/5.0, 0.5*sqrt(3)*1.0/5.0 nohead linestyle 2

# z
set arrow 7 from 1, 0 to 0.50, 0.866 nohead linestyle 1
set arrow 8 from 1.0/5.0, 0.0 to 0.5*1.0/5.0, 0.5*sqrt(3)*1.0/5.0 nohead linestyle 2
set arrow 9 from 2.0/5.0, 0.0 to 0.5*2.0/5.0, 0.5*sqrt(3)*2.0/5.0 nohead linestyle 2
set arrow 10 from 3.0/5.0, 0.0 to 0.5*3.0/5.0, 0.5*sqrt(3)*3.0/5.0 nohead linestyle 2
set arrow 11 from 4.0/5.0, 0.0 to 0.5*4.0/5.0, 0.5*sqrt(3)*4.0/5.0 nohead linestyle 2

# y
set arrow 13 from 0.50, 0.866 to 0,0 nohead linestyle 1
set arrow 14 from 0.5*1.0/5.0, 0.5*sqrt(3)*1.0/5.0 to 1.0-0.5*1.0/5.0, 0.5*sqrt(3)*1.0/5.0 nohead linestyle 2
set arrow 15 from 0.5*2.0/5.0, 0.5*sqrt(3)*2.0/5.0 to 1.0-0.5*2.0/5.0, 0.5*sqrt(3)*2.0/5.0 nohead linestyle 2
set arrow 16 from 0.5*3.0/5.0, 0.5*sqrt(3)*3.0/5.0 to 1.0-0.5*3.0/5.0, 0.5*sqrt(3)*3.0/5.0 nohead linestyle 2
set arrow 17 from 0.5*4.0/5.0, 0.5*sqrt(3)*4.0/5.0 to 1.0-0.5*4.0/5.0, 0.5*sqrt(3)*4.0/5.0 nohead linestyle 2

#plot "<awk '{if (\$1 == 1) {print (\$2+2*\$4)/(2*(\$2+\$3+\$4)), sqrt(3)*\$2/(2*(\$2+\$3+\$4)), \$5, 1+8*(\$5)/$maxP}}' $tmpfile" u 1:2:4:3 w p pt 7 ps variable lt palette
plot "<awk '{if (\$1 == 1) {print (\$2+2*\$4)/(2*(\$2+\$3+\$4)), sqrt(3)*\$2/(2*(\$2+\$3+\$4)), 100*((\$5)-(\$6)), 1+12*(\$7)/$maxABSP}}' $tmpfile1" u 1:2:4:3 w p pt 7 ps variable lt palette

#######################################################################################################################

unset arrow

set bmargin at screen 0.41
set lmargin at screen 0.45
set rmargin at screen 0.62
set tmargin at screen 0.63

# x
set arrow 1 from 0,0 to 1, 0.0 nohead linestyle 1
set arrow 2 from 1.0/4.0, 0.0 to 1.0-0.5*3.0/4.0, 0.5*sqrt(3)*3.0/4.0 nohead linestyle 2
set arrow 3 from 2.0/4.0, 0.0 to 1.0-0.5*2.0/4.0, 0.5*sqrt(3)*2.0/4.0 nohead linestyle 2
set arrow 4 from 3.0/4.0, 0.0 to 1.0-0.5*1.0/4.0, 0.5*sqrt(3)*1.0/4.0 nohead linestyle 2

# z
set arrow 7 from 1, 0 to 0.50, 0.866 nohead linestyle 1
set arrow 8 from 1.0/4.0, 0.0 to 0.5*1.0/4.0, 0.5*sqrt(3)*1.0/4.0 nohead linestyle 2
set arrow 9 from 2.0/4.0, 0.0 to 0.5*2.0/4.0, 0.5*sqrt(3)*2.0/4.0 nohead linestyle 2
set arrow 10 from 3.0/4.0, 0.0 to 0.5*3.0/4.0, 0.5*sqrt(3)*3.0/4.0 nohead linestyle 2

# y
set arrow 13 from 0.50, 0.866 to 0,0 nohead linestyle 1
set arrow 14 from 0.5*1.0/4.0, 0.5*sqrt(3)*1.0/4.0 to 1.0-0.5*1.0/4.0, 0.5*sqrt(3)*1.0/4.0 nohead linestyle 2
set arrow 15 from 0.5*2.0/4.0, 0.5*sqrt(3)*2.0/4.0 to 1.0-0.5*2.0/4.0, 0.5*sqrt(3)*2.0/4.0 nohead linestyle 2
set arrow 16 from 0.5*3.0/4.0, 0.5*sqrt(3)*3.0/4.0 to 1.0-0.5*3.0/4.0, 0.5*sqrt(3)*3.0/4.0 nohead linestyle 2

#plot "<awk '{if (\$1 == 2) {print (\$2+2*\$4)/(2*(\$2+\$3+\$4)), sqrt(3)*\$2/(2*(\$2+\$3+\$4)), \$5, 1+8*(\$5)/$maxP}}' $tmpfile" u 1:2:4:3 w p pt 7 ps variable lt palette
plot "<awk '{if (\$1 == 2) {print (\$2+2*\$4)/(2*(\$2+\$3+\$4)), sqrt(3)*\$2/(2*(\$2+\$3+\$4)), 100*((\$5)-(\$6)), 1+12*(\$7)/$maxABSP}}' $tmpfile1" u 1:2:4:3 w p pt 7 ps variable lt palette

#######################################################################################################################

unset arrow

set bmargin at screen 0.52
set lmargin at screen 0.58
set rmargin at screen 0.70
set tmargin at screen 0.68

# x
set arrow 1 from 0,0 to 1, 0.0 nohead linestyle 1
set arrow 2 from 1.0/3.0, 0.0 to 1.0-0.5*2.0/3.0, 0.5*sqrt(3)*2.0/3.0 nohead linestyle 2
set arrow 3 from 2.0/3.0, 0.0 to 1.0-0.5*1.0/3.0, 0.5*sqrt(3)*1.0/3.0 nohead linestyle 2

# z
set arrow 7 from 1, 0 to 0.50, 0.866 nohead linestyle 1
set arrow 8 from 1.0/3.0, 0.0 to 0.5*1.0/3.0, 0.5*sqrt(3)*1.0/3.0 nohead linestyle 2
set arrow 9 from 2.0/3.0, 0.0 to 0.5*2.0/3.0, 0.5*sqrt(3)*2.0/3.0 nohead linestyle 2

# y
set arrow 13 from 0.50, 0.866 to 0,0 nohead linestyle 1
set arrow 14 from 0.5*1.0/3.0, 0.5*sqrt(3)*1.0/3.0 to 1.0-0.5*1.0/3.0, 0.5*sqrt(3)*1.0/3.0 nohead linestyle 2
set arrow 15 from 0.5*2.0/3.0, 0.5*sqrt(3)*2.0/3.0 to 1.0-0.5*2.0/3.0, 0.5*sqrt(3)*2.0/3.0 nohead linestyle 2

#plot "<awk '{if (\$1 == 3) {print (\$2+2*\$4)/(2*(\$2+\$3+\$4)), sqrt(3)*\$2/(2*(\$2+\$3+\$4)), \$5, 1+8*(\$5)/$maxP}}' $tmpfile" u 1:2:4:3 w p pt 7 ps variable lt palette
plot "<awk '{if (\$1 == 3) {print (\$2+2*\$4)/(2*(\$2+\$3+\$4)), sqrt(3)*\$2/(2*(\$2+\$3+\$4)), 100*((\$5)-(\$6)), 1+12*(\$7)/$maxABSP}}' $tmpfile1" u 1:2:4:3 w p pt 7 ps variable lt palette

#######################################################################################################################

unset arrow

set bmargin at screen 0.61
set lmargin at screen 0.68
set rmargin at screen 0.75
set tmargin at screen 0.71

# x
set arrow 1 from 0,0 to 1, 0.0 nohead linestyle 1
set arrow 2 from 1.0/2.0, 0.0 to 1.0-0.5*1.0/2.0, 0.5*sqrt(3)*1.0/2.0 nohead linestyle 2

# z
set arrow 7 from 1, 0 to 0.50, 0.866 nohead linestyle 1
set arrow 8 from 1.0/2.0, 0.0 to 0.5*1.0/2.0, 0.5*sqrt(3)*1.0/2.0 nohead linestyle 2

# y
set arrow 13 from 0.50, 0.866 to 0,0 nohead linestyle 1
set arrow 14 from 0.5*1.0/2.0, 0.5*sqrt(3)*1.0/2.0 to 1.0-0.5*1.0/2.0, 0.5*sqrt(3)*1.0/2.0 nohead linestyle 2

#plot "<awk '{if (\$1 == 4) {print (\$2+2*\$4)/(2*(\$2+\$3+\$4)), sqrt(3)*\$2/(2*(\$2+\$3+\$4)), \$5, 1+8*(\$5)/$maxP}}' $tmpfile" u 1:2:4:3 w p pt 7 ps variable lt palette
plot "<awk '{if (\$1 == 4) {print (\$2+2*\$4)/(2*(\$2+\$3+\$4)), sqrt(3)*\$2/(2*(\$2+\$3+\$4)), 100*((\$5)-(\$6)), 1+12*(\$7)/$maxABSP}}' $tmpfile1" u 1:2:4:3 w p pt 7 ps variable lt palette

#######################################################################################################################

unset arrow

set bmargin at screen 0.67
set lmargin at screen 0.75
set rmargin at screen 0.80
set tmargin at screen 0.74

# x
set arrow 1 from 0,0 to 1, 0.0 nohead linestyle 1

# z
set arrow 7 from 1, 0 to 0.50, 0.866 nohead linestyle 1

# y
set arrow 13 from 0.50, 0.866 to 0,0 nohead linestyle 1

set arrow 20 from screen 0.20, 0.435 to screen 0.90,0.80 nohead linestyle 3
set arrow 21 from screen 0.35, 0.05 to screen 0.90,0.80 nohead linestyle 3

set label 1 '(0,0,6,0) Pure "N" (POPE)' font ",24" at screen 0.10, 0.03 center offset character 0,-1
set label 2'(0,0,0,6) Pure "C" (POPS)' font ",24" at screen 0.30, 0.03 center offset character 2,-1
set label 3 '(0,6,0,0) Pure "P" (POPC)' font ",24" at screen 0.20, 0.45 center offset character 0,1.5
set label 4 '(6,0,0,0) Pure "O" (CHOL)' font ",24" at screen 0.85, 0.80 center offset character 0,1.5

#plot "<awk '{if (\$1 == 5) {print (\$2+2*\$4)/(2*(\$2+\$3+\$4)), sqrt(3)*\$2/(2*(\$2+\$3+\$4)), \$5, 1+8*(\$5)/$maxP}}' $tmpfile" u 1:2:4:3 w p pt 7 ps variable lt palette
plot "<awk '{if (\$1 == 5) {print (\$2+2*\$4)/(2*(\$2+\$3+\$4)), sqrt(3)*\$2/(2*(\$2+\$3+\$4)), 100*((\$5)-(\$6)), 1+12*(\$7)/$maxABSP}}' $tmpfile1" u 1:2:4:3 w p pt 7 ps variable lt palette

#######################################################################################################################


unset multiplot
EOF

exit
