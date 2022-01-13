#!/bin/bash

#
# Example use of this script:
# ./getHMManalysis.sh
#

tmpfile1=tmp1
tmpfile2=tmp2

# Input file with the observable states of each
# lipid for each frame line-by-line...
# The order of the lipids must be kept consistent
# over each line (each frame)... if a lipid is not
# on the leaflet of interest, record it as a "0",
hmminputfile=hmm-only-sidegroup-head-chl-39995.txt

##############################################################################

# Minimum length of a sequence of observables ...
# this means that lipids that only spend a short
# length of time on the leaflet of interest (e.g.
# cholesterol that switches leaflets) will not
# be optimized over in the HMM
n_measure_min=800

# Number of hidden states
Nhidden=2

# Number of observables states
Nobservable=84

##############################################################################

gfortran HMMsolverMultipleSequences_KF.f95 -o HMMsolverMultipleSequences_KF.o

cat $hmminputfile > $tmpfile1

Nlipids=$(head -1 $tmpfile1 | xargs -n1 | wc -l)
Nlines=$(cat $tmpfile1 | wc -l)

# Get the sequences of states for each lipid;
# If the lipid gets 'cut out', then there may
# be multiple sequences
awk "BEGIN {for (i=1;i<=$Nlipids;i++) {lengths[i]=0;sequences[i]=\"\"}} {sum=0; for (i=1;i<=NF;i++) {if (\$i > 0) {sum+=1;lengths[i]+=1;sequences[i]=sequences[i]\" \"\$i; if (NR == $Nlines) {print i\":\", sequences[i]}} else {if (lengths[i] > 0) {print i\":\", sequences[i]}; lengths[i]=0}}}" $tmpfile1 > $tmpfile2

# Get rid of any sequences that are too short
awk "{if (NF > $n_measure_min) {\$1=\"\"; print \$0}}" $tmpfile2 > $tmpfile1
mv $tmpfile1 $tmpfile2
n_sequence=$(cat $tmpfile2 | wc -l)

date
echo "   Nlipids: $Nlipids"
echo "    Nlines: $Nlines"
echo "Nsequences: $n_sequence"

# Process these and get the HMM model
./HMMsolverMultipleSequences_KF.o $Nhidden $Nobservable $n_sequence < <(awk '{print NF; for (i=1;i<=NF;i++) {print $i}}' $tmpfile2)

date

##############################################################################

