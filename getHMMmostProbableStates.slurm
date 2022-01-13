#!/bin/bash

#
# Example use of this script:
# ./getHMMmostProbableStaes.sh > getHMMmostProbableStaes.out
#

tmpfile1=tmp1
tmpfile2=tmp2

# Supply the original PDB file to figure out
# the ordering of lipid types
originalpdb=example_frame.pdb

# Supply the original HMM input file
hmminputfile=hmm-observables.txt

# Supply the output of the HMM analysis
hmmanalysisfile=getHMManalysis.out

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

# Do some pre-processing

cat $hmminputfile > $tmpfile1

Nlipids=$(grep '^ATOM ' $originalpdb | wc -l)
Nlines=$(cat $tmpfile1 | wc -l)

date
echo "Nlipids: $Nlipids"
echo "Nlines: $Nlines"
echo ""

# Get the sequences of states for each lipid;
# If the lipid gets 'cut out', then there may
# be multiple sequences
awk "BEGIN {for (i=1;i<=$Nlipids;i++) {lengths[i]=0;sequences[i]=\"\"}} {sum=0; for (i=1;i<=NF;i++) {if (\$i > 0) {sum+=1;lengths[i]+=1;sequences[i]=sequences[i]\" \"\$i; if (NR == $Nlines) {print 1+NR-lengths[i], i\":\", sequences[i]}} else {if (lengths[i] > 0) {print NR-lengths[i], i\":\", sequences[i]; sequences[i]=\"\"}; lengths[i]=0}}}" $tmpfile1 > $tmpfile2

mv $tmpfile2 $tmpfile1

# Get rid of any sequences that are too short
awk "{if (NF > $n_measure_min) {\$1=\"\"; \$2=\"\"; print \$0}}" $tmpfile1 > $tmpfile2
n_sequence=$(cat $tmpfile2 | wc -l)

date
echo "Nsequences: $n_sequence"

##############################################################################

gfortran ViterbiMostProbableSequence_KF.f95 -o ViterbiMostProbableSequence_KF.o

# Get the most probable hidden states for
# each sequence

cat $tmpfile1 | while read aline; do
  identifier=$(echo "$aline" | awk '{print $1, $2}')
  line=$(echo "$aline" | awk "{\$1=\"\"; \$2=\"\"; print \$0}")
  n_sequence=$(echo "$line" | awk '{print NF}')
  (echo "$identifier"; 
  ./ViterbiMostProbableSequence_KF.o $Nhidden $Nobservable $n_sequence < <(awk '/best model/ {score=$NF; getline; getline; getline; line=$0; getline; getline; line=line$0; getline; line=line$0; getline; getline; line=line$0; if (getline <= 0) {exit}; line=line$0; print score, line}' $hmmanalysisfile | sort -g -k1 | tail -1 | awk '{printf "%s\n", $2; printf "%s\n", $3; printf "%s %s\n", $4, $5; printf "%s %s\n", $6, $7; n=(NF-7)/2; for (i=8;i<=7+n;i++) {printf "%s ", $i}; printf "\n"; for (i=8+n;i<=NF;i++) {printf "%s ", $i}; printf "\n"}'; echo "$line" | xargs -n1) | sed 1,1d
  ) | xargs
done

date

##############################################################################

