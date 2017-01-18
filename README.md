# Machete
Machete: Likelihood decay indices using PAUP

## Overview

Bootstrap proportion (BP) support remains a commonly used metric of the reliability of genome-scale phylogenetic analyses because sampling error decreases as the length of sequences increase resulting in a trend where BP support approaches 100%. However, not all conflicting phylogenetic signal is due to sampling error; processes such as incomplete lineage sorting and horizontal gene transfer can result in valid alternative genetic histories.

Despite this, with long-enough alignments, 100% BP can be achieved even if 49% of the data supports an alternative topology. The heterogeneous nature of the underlying support for branches with 100% BP requires a novel approach and a change in our notion of "support".

To address this, Machete implements a likelihood decay support value. Based on the premise of Bremmer support, it is the difference in likelihoods of the optimal trees that do or do not include a given split. Likelihood decay represents a novel way to assess support which discriminates between different internal branches and is insensitive to alignment length. 

Machete takes as input a nexus formatted aligned DNA or Amino Acid sequences and uses PAUP to calculate the likelihood decay supports for each internal branch of the resulting tree.

Machete controls and interacts with Paup using a pipe, and not using a predefined script. This allows dataset-specific optimisations to be carried out (as a user would).

## Algorithm

Machete carries out the following steps:

1. Identifies the phylogeny best supported by the inputted alignment by:
  * Carrying out a heuristic search of tree space (unless a tree is provided - see below).
  * Optimising model parameters using best tree.
  * Repeat previous steps until likelihoods stops improving.

2. Define all constraints based on each of the internal branches of the best tree.

3. For each constraint, carry out a search for the the best tree that does not contain that branch.

5. Calculate the likelihood decay support for each internal branch and output the best phylogeny with the support values labelled on each branch


## Installation

Download the machete.c file above, or if you have git installed, type the command 

`git clone https://github.com/ChrisCreevey/machete.git`

Run the following command to create a version of machete for your system (linux or macosx).

```
cc machete.c -o machete -lm
```

You will also need a copy of paup. A binary for your system can be downloaded at [https://people.sc.fsu.edu/~dswofford/paup_test/](https://people.sc.fsu.edu/~dswofford/paup_test/).

This will need to be renamed "paup" and may need to be made executable:

```
mv paupVERSION_SYSTEM paup
chmod a+x paup
```

It would be best to move both paup and machete to somewhere on your path (like ~/bin) to make both available everywhere in the system.

## Using machete

To use machete, it is necessary to pass as input a nexus formatted alignment file using the option '-f':

```
machete -f NEXUSFILE
```

An example file `Primate.nex` has been provided. To run the reverse constraint analysis for this file use the command:

```
machete -f Primate.nex
```
You can provide a pre-calculated phylogeny to machete by adding a "trees" block to the end of the nexus file.
If this is provided, machete will not try to build the optimum tree, but use the provided tree instead.
If you wish to over-ride this functionality, use the -b option (see below).

Other options are:

 -c which tells Machete to print all the commands sent to Paup to the standard error. This can then be redirected to a seperate file using the following syntax:
 
 ```
 machete -f Primate.nex -c 2> paupcommands.txt
 ```
  -t preserves temporary files
  
  -b forces machete to build an optimum tree when a tree has been provided in the nexus file
 
  -h prints a description of the usage of machete.
 
 	-s <constraint number> specifies the constraint at which to start the analysis
  
	 -e <constraint number> specific the constraint at which to end the analysis
  
	 -l list constraints (and do not carry out reverse constraints analysis)
 
## Outputs

Machete produces three standard output files:

```
<NEXUSFILE>.labelledtree.tre
<NEXUSFILE>.sitelike.txt
<NEXUSFILE>.constraint.tre
```

*Where* `<NEXUSFILE>` is the name of the alignment input file.

`<NEXUSFILE>.labelledtree.tre` will contain the optimimum phylogeny either calculated by PAUP or provided by the user with internal branch labels in the following format: `2/15.19/41.0945/25.9042/0.613363` where:

>2 = the internal branch ID (and constraint number)
  
>15.19 = The difference in the -lnL between the unconstrained tree and constrained tree for this branch.
  
>41.0945 = The absolute value of the sum of likelihood differences supporting the unconstrained tree at this branch.
  
>25.9042 = The absolute value of the sum of the likelihood differences supporting the constrained tree at this branch.
  
>0.613363 = Proportion of the absolute likelihood differences supporting the unconstrained tree {calculated in this case as: (41.0945)/(41.0945+25.9042)}

These supports can be viewed on the tree using a phylogeny viewer such as [figtree](http://tree.bio.ed.ac.uk/software/figtree/).

`<NEXUSFILE>.sitelike.txt` contains all the sitelikelihoods calcuated for the unconstrained tree and for each of the constrained trees (labelled by internal branch ID).

`<NEXUSFILE>.constraint.tre` contains all the best constrained trees for each internal branch in phylip format (labelled by internal branch ID).






