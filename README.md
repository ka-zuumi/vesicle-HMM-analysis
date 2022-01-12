# vesicle-HMM-analysis

## Purpose

After gathering molecular dynamics (MD) trajectories of a biophysical simulation of a vesicle, further analysis can be done to distinguish hidden states, such as lipid ordered and disordered states. Consider a simplified vesicle of only the lipid head groups; a sample frame "example_frame.pdb" is given which can be visualized with a pre-made VMD tcl script:

```
vmd example_frame.pdb -e <(echo "source lookAtVesicle.tcl")
```

![Alt text](vesicle-analysis-visualaid1-bare1.png?raw=true "Vesicle Description")

After simulating many consecutive frames of this mixture of lipids, it may be discernable to the human eye that small or medium sized domains now and then transiently form. To have a more concrete definition for these domains, a hidden markov state model (HMM) may be made to differentiate lipids of different domains (hidden states) through observation of key properties (obervable states) over time. The coordinates in the PDB files are sufficient to record a few key observables states, namely the local lipid composition (LCC).

## Local Lipid Composition

Each lipid has a number of "neighboring" lipids which it borders. This local lipid composition (LCC) varies in number (roughly six neighbors for hexagonally-packed lipids) and composition (i.e. cholesterol-rich, POPC-rich). To simplify the diversity of LCCs, consider only the five closest lipids and the lipid itself; this composition of six lipid types can be uniquely assigned an integer value distinct from every other composition of six lipid types. Altogether, all LCCs can be described by one of 84 compositions, as is summarized in the figure below.

![Alt text](vesicle-analysis-visualaid1.1.png?raw=true "Vesicle LCC Breakdown")

Each composition has some probability of occurring, as can be approximated by a Bernoulli multivariate distribution over the four lipid types for this vesicle. For this vesicle's outer membrane's specific composition, these probabilities as displayed in the figure above can be calculated and visualized as follows:

```
gfortran getBernoulliProbability.f90 -o getBernoulliProbability.o;
./makeTetrahedronGraph.sh <(seq 0 6 | while read nO; do seq 0 6 | while read nP; do seq 0 6 | while read nN; do let "nC = 6 - ($nO + $nP + $nN)"; if [ "$nC" -ge "0" ]; then P=$(./getBernoulliProbability.out $nP $nN $nC $nO | awk '{print 100 * $1}'); echo "$nO $nP $nN $nC $P"; fi; done; done; done) bernoulli-triangle.png "Bernoulli Probability Distribution"
```

To use this observable for the HMM, coordinates from the PDB files must first be processed to record the LCC of each lipid as frames listed line-by-line as follows:

```
python BLAHBLAH example_frame.pdb > hmm-observables.txt
```
