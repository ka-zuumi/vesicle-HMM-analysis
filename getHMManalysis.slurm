#!/bin/bash
#SBATCH --job-name=vesicleAnalysis-4lipid-system-5and1-analysis2
#SBATCH --partition=rsun
#SBATCH --account=rsun
##SBATCH --exclusive
##SBATCH --constraint=ib_qdr
##SBATCH -w, --nodelist=node-0007
#SBATCH --exclude=gpu-[0005-0007]
#SBATCH --nodes=1
#SBATCH --tasks-per-node=9
#SBATCH --cpus-per-task=1
#SBATCH --error="4lipid-system-5and1-analysis2.err"
#SBATCH --time=3-00:00:00 ## time format is DD-HH:MM:SS
#SBATCH --output="4lipid-system-5and1-analysis2.out"
#SBATCH --mem=45000

#
# Example use of this script:
# ./getHMManalysis.sh hmm-only-sidegroup-head-chl.txt
#

tmpfile1=allsequences2
tmpfile2=tmpsequences2

n_measure_min=800

#From: /anton2fs/raw/cakang/gromacs_popc_chl/step_test/only_heads/frame00001
originalpdb=/home/kazuumi/rsun_lts/kazuumi/vesiclestuff/4lipid-system-5and1-2hiddenstates-analysis/frame00001-only-sidegroup-head-chl.pdb

##############################################################################

# bash: supply the arguments from the command line

#if [ $# -ne "2" ]; then
#  echo "ERROR. Two arguments required:"
#  echo "       (1) Lipid HMM states file"
#  echo "       (2) Lipid SA file"
#  exit 1
#fi
#
## The first argument is the input
## file with the HMMs of each lipid
## on each line for each frame
#hmminputfile=$1

###############################################

# sbatch: just supply the arguments

hmminputfile=hmm-only-sidegroup-head-chl-39995.txt


##############################################################################

module load tools/Vim
gfortran HMMsolverMultipleSequences_KF.f95 -o a.out

cat $hmminputfile > $tmpfile1


Nlipids=$(head -1 $tmpfile1 | xargs -n1 | wc -l)
Nlines=$(cat $tmpfile1 | wc -l)

# Number of each lipid type (outer/total):
#   394/562 CHL1 - type 4 - letter "O"
#   113/397 POPC - type 1 - letter "P"
#   392/434 POPE - type 2 - letter "N"
#   294/294 POPS - type 3 - letter "C"

# Get the sequences of states for each lipid;
# If the lipid gets 'cut out', then there may
# be multiple sequences
awk "BEGIN {for (i=1;i<=$Nlipids;i++) {lengths[i]=0;sequences[i]=\"\"}} {sum=0; for (i=1;i<=NF;i++) {if (\$i > 0) {sum+=1;lengths[i]+=1;sequences[i]=sequences[i]\" \"\$i; if (NR == $Nlines) {print i\":\", sequences[i]}} else {if (lengths[i] > 0) {print i\":\", sequences[i]}; lengths[i]=0}}}" $tmpfile1 > $tmpfile2

# Get rid of any sequences that are too short
#awk "{if (NF > $n_measure_min) {\$1=\"\"; print \$0}}" $tmpfile2 | tail -1000 > $tmpfile1 # Only look at the last 1000 sequences
awk "{if (NF > $n_measure_min) {\$1=\"\"; print \$0}}" $tmpfile2 > $tmpfile1
mv $tmpfile1 $tmpfile2
n_sequence=$(cat $tmpfile2 | wc -l)

date
echo "   Nlipids: $Nlipids"
echo "    Nlines: $Nlines"
echo "Nsequences: $n_sequence"

# Process these and get the HMM model
./a.out 2 84 $n_sequence < <(awk '{print NF; for (i=1;i<=NF;i++) {print $i}}' $tmpfile2)

date

##############################################################################

