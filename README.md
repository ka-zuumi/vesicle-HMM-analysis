# vesicle-HMM-analysis

## Purpose

After gathering molecular dynamics (MD) trajectories of a biophysical simulation of a vesicle, further analysis can be done to distinguish hidden states, such as lipid ordered and disordered states. Consider a simplified vesicle of only the lipid head groups; a sample frame "example_frame.pdb" is given which can be visualized with a pre-made VMD tcl script:

```
vmd example_frame.pdb -e <(echo "source lookAtVesicle.tcl")
```

![Alt text](vesicle-analysis-visualaid1-bare1.png?raw=true "Vesicle Description")

After simulating many consecutive frames of this mixture of lipids, it may be discernable to the human eye that small or medium sized domains now and then transiently form. To have a more concrete definition for these domains, a hidden markov state model (HMM) may be made to differentiate lipids of different domains (hidden states) through observation of key properties (obervable states) over time. The coordinates in the PDB files are sufficient to record a few key observables states, namely the local lipid composition (LCC).

## Local Lipid Composition

Each lipid has a number of "neighboring" lipids which it borders. This local lipid composition (LCC) varies in number (e.g. six neighbors for hexagonally-packed lipids) and composition (e.g. cholesterol-rich, POPC-rich). To simplify the diversity of LCCs, consider only the five closest lipids and the lipid itself; this composition of six lipid types can be uniquely assigned an integer value distinct from every other composition of six lipid types. Altogether, all LCCs can be described by one of 84 compositions, as is summarized in the figure below.

![Alt text](vesicle-analysis-visualaid1.1.png?raw=true "Vesicle LCC Breakdown")

Each composition has some probability of occurring and are approximated by a Bernoulli multivariate distribution over the four lipid types for this vesicle in the figure above. For this vesicle's outer membrane's specific composition, these probabilities can be calculated and visualized as follows:

```
gfortran getBernoulliProbability.f90 -o getBernoulliProbability.o;
./makeTetrahedronGraph.sh <(seq 0 6 | while read nO; do seq 0 6 | while read nP; do seq 0 6 | while read nN; do let "nC = 6 - ($nO + $nP + $nN)"; if [ "$nC" -ge "0" ]; then P=$(./getBernoulliProbability.out $nP $nN $nC $nO | awk '{print 100 * $1}'); echo "$nO $nP $nN $nC $P"; fi; done; done; done) bernoulli-triangle.png "Bernoulli Probability Distribution"
```

To use this observable for the HMM, coordinates from the PDB files must first be processed to record the LCC of each lipid as frames listed line-by-line as follows:

```
python BLAHBLAH example_frame.pdb > hmm-observables.txt
```

## Hidden Markov State Model

The full list of observables states of each lipid for all frames (here, not given due to the large size), can then be fed into a HMM. While any number of hidden states may be specified, by specifying two, less variables will be needed to optimize over while still differentiating between the two most distinct hidden states. This can be done with the "getHMManalysis.slurm" script, after specifying the inputfile as:

```
hmminputfile=hmm-observables.txt
```

As this may take some time (an hour to test a single set of initial conditions), execution of this script is normally set up to be done in the background. When a number of different converged models are found with this (each from separate initial conditions), the output of the script "getHMManalysis.out" may be processed to compare the different models:

```
./compareHMModels.sh getHMManalysis.out HMMcomparison.png "4 Lipids - 2 Hidden States"
```

In this example, the output found twenty models before termination. These models are probability distributions and can be compared to one another by their root mean square deviation (RMSD). The score, as displayed in the figure below, is the logarithm of the probability that a model produced the sequence of observable states that were initially given. In this way, the highest-scoring model is the most-probable model. The low RMSD between the most-probable models is one indicator of model consensus agreement (and possible convergence to the global maximally-probable model).

![Alt text](HMMcomparison.png?raw=true "HMM comparison")
