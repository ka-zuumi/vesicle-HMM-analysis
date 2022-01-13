#!/bin/bash

#
# Example usage:
# ./compareHMModels.sh 4lipid-system-5and1-analysis2.out HMMcomparison.png "Simulated Vesicle - 10000 Frames - 4 Lipid - 2 Hidden States"
#

tmpfile1=tmp1
tmpfile2=tmp2


############################################################################
#
# Warning: Right now, this only supports models with
#          exactly two hidden states
#

if [ $# -gt "3" ] || [ $# -lt 2 ]; then
  echo ""
  echo "Error: Requires exactly two or three arguments:"
  echo "       1) HMM output file"
  echo "       2) Name of output image file (must be .PNG)"
  echo "May have an optional third argument:"
  echo "       3) Text to put on image (to identify the graph)"
  echo "STOPPING"
  exit
fi

if [ $# -eq "3" ]; then
  gnuplotline=$3
else
  gnuplotline=""
fi

imagename=$2

# Process the HMM output file to get each model as a line on the tmpfile

# To look at ALL 'best model' iterations, do:
#awk '/best model/ {score=$NF; getline; getline; getline; line=$0; getline; getline; line=line$0; getline; line=line$0; getline; getline; line=line$0; if (getline <= 0) {exit}; line=line$0; print score, line}' $1 > $tmpfile1

# To just get one 'best model' per initial guess, do:
awk 'BEGIN {iteration=0} /iteration/ {iteration=$NF} /best model/ {score=$NF; getline; getline; getline; line=$0; getline; getline; line=line$0; getline; line=line$0; getline; getline; line=line$0; if (getline <= 0) {exit}; line=line$0; print iteration, score, line}' $1 | tac | awk 'BEGIN {iteration=99999} {if ($1 < iteration) {iteration=$1; $1=""; print $0}}' > $tmpfile1

# Sort this by score
sort -g -k1 $tmpfile1 > $tmpfile2

# Then store the scores as a separate variable
awk 'BEGIN {scores=""} {scores=scores" "$1; $1=""; print $0} END {print scores}' $tmpfile2 > $tmpfile1
scores=$(tail -1 $tmpfile1)
sed -i '$ d' $tmpfile1

# The number of parameters to compare between models is
# exactly how many numbers are spit out each line while
# each model takes up a line
Nparameters=$(head -1 $tmpfile1 | awk '{print NF}')
Nmodels=$(cat $tmpfile1 | wc -l)

# 2 parameters for initial probabilities
# 4 parameters for transition probabilities
# and two sets of observable probabilities
let "Nobservables = ($Nparameters - 6) / 2"

############################################################################

# Spit out the highest scoring model for later
line=$(tail -1 $tmpfile2 | awk '{$1=""; print $0}')
pi_i=$(echo "$line" | awk "{print \$1, \$2}")
aij=$(echo "$line" | awk "{print \$3, \$4, \$5, \$6}")
bik=$(echo "$line" | awk "{\$1=\"\"; \$2=\"\"; \$3=\"\"; \$4=\"\"; \$5=\"\"; \$6=\"\"; print \$0}")

echo "$pi_i" | xargs -n1
echo "$aij" | xargs -n2
echo "$bik" | xargs -n$Nobservables

############################################################################

gfortran makeDistanceMatrix.f90 -o makeDistanceMatrix.out
./makeDistanceMatrix.out $Nparameters $Nmodels < $tmpfile1 > $tmpfile2

maxP=$(awk 'BEGIN {maxP=0.0} {for (i=1;i<=NF;i++) {if ($i > maxP) {maxP=$i}}} END {print maxP}' $tmpfile2)

echo "$scores" | xargs -n1 > $tmpfile1
maxScore=$(tail -1 $tmpfile1)
minScore=$(head -1 $tmpfile1)

############################################################################

# Let's visualize the difference matrix


module load vis/gnuplot/5.2.6-foss-2018b
gnuplot <<- EOF
set terminal pngcairo size 1400,1400
set output '$imagename'

set multiplot
set noborder
set noxtics
set noytics
unset key
set label "$gnuplotline " at screen 0.02,0.98
set label "Number of Models: $Nmodels" at screen 0.02,0.02

set bmargin at screen 0.05
set lmargin at screen 0.10
set rmargin at screen 0.85
set tmargin at screen 0.80

set cbrange [0:$maxP]
set cbtics font ",12"
set cblabel "RMSD Between HMM Probabilities" font ",16" offset 2,0
#splot '$tmpfile2' matrix w image
plot '$tmpfile2' matrix w image

set border 3 lw 2

set bmargin at screen 0.85
set lmargin at screen 0.10
set rmargin at screen 0.85
set tmargin at screen 0.95

set xrange [1-0.5:$Nmodels+0.5]
unset xlabel
set yrange [$minScore:$maxScore]
set ylabel "Score (logP)" font ",16" offset -1,0

plot "$tmpfile1" u ((\$0)+1):1 w l lw 2 lc rgb "red"
unset multiplot
EOF
