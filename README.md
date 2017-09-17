# Machete
Machete: Automated Maximum likelihood phlyogeny construction with PAUP* and calculation of likelihood decay indices.

## Sections

 [Overview](README.md#overview)

 [Rationale](README.md#rationale)

 [Quick Start](README.md#quick-start)
 
 [Algorithm](README.md#algorithm)

 [Installation](README.md#installation)

 [Using machete and options](README.md#using-machete-and-options)

 [Outputs](README.md#outputs)



## Overview
Machete takes as input a nexus formatted aligned DNA or Amino Acid sequences and uses PAUP to automatically calculate maximum likelihood trees (and/or carry out boostrap analyses) while optimising the models. It has been desinged to allow calculation of the likelihood decay supports for each internal branch of the resulting tree.

Machete controls and interacts with Paup using a pipe, and not using a predefined script. This allows dataset-specific optimisations to be carried out (as a user would).



### Rationale

Bootstrap proportion (BP) support remains a commonly used metric of the reliability of genome-scale phylogenetic analyses because sampling error decreases as the length of sequences increase resulting in a trend where BP support approaches 100%. However, not all conflicting phylogenetic signal is due to sampling error; processes such as incomplete lineage sorting and horizontal gene transfer can result in valid alternative genetic histories.

Despite this, with long-enough alignments, 100% BP can be achieved even if 49% of the data supports an alternative topology. The heterogeneous nature of the underlying support for branches with 100% BP requires a novel approach and a change in our notion of "support".

To address this, Machete implements a likelihood decay support value. Based on the premise of Bremmer support, it is the difference in likelihoods of the optimal trees that do or do not include a given split. Likelihood decay represents a novel way to assess support which discriminates between different internal branches and is insensitive to alignment length. 



## Quick Start

### To generate a maximum likeihood tree (DNA or Amino Acid) in PAUP* (while optimising for best models etc) 
```
machete -f alignment.nexus -n
```

### To generate a maximum likelihood tree AND carry out a bootstrap analysis in PAUP* (while optimising for best models etc)
```
machete -f alignment.nexus -n -r 100
```

### To calculate likelihood decay for best tree (given only an alignment) - will also calculate ML tree
```
machete -f alignment.nexus 
```

### To calculate likelihood decay for Given tree - will optimise models to tree provided.
```
machete -f alignment_with_treeblock.nexus
```
See the options description below for more details.




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

You will also need a copy of paup. A binary for your system can be downloaded at [http://paup.phylosolutions.com/](http://paup.phylosolutions.com/).

This will need to be renamed "paup" and may need to be made executable:

```
mv paupVERSION_SYSTEM paup
chmod a+x paup
```

It would be best to move both paup and machete to somewhere on your path (like ~/bin) to make both available everywhere in the system.




## Using machete and options

Calling Machete without any options will result in the output of the options avaiable:

```
Machete: Likelihood reverse constraint analysis using PAUP

 Usage: "machete -f <nexus file> -[cthln] [-s INTEGER] [-e INTEGER] [-r INTEGER]"

	Where: <nexus file> is a nexus formatted alignment file of DNA sequences
	-c commands sent to Paup to be also printed to standard error
	-t preserves temporary files
	-h prints this message
	-b force build optimum tree (when a tree has been provided in the nexus file)
	-s <constraint number> specifies the constraint to start at
	-e <constraint number> specific the constraint to end at
	-l list constraints (and do not carry out reverse constraints analysis)
	-n tells machete NOT to carry out the reverse constraint analysis (Just build the best tree)
	-r tells machete how many boostrap replicates to carry out (by default = 0)
 
 ```
 
In more detail:

### -f 

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


### -c

-c which tells Machete to print all the commands sent to Paup to the standard error. This can then be redirected to a seperate file using the following syntax:
 
 ```
 machete -f Primate.nex -c 2> paupcommands.txt
 ```
 
### -t
 
  -t preserves temporary files
 
### -b
 
  -b forces machete to build an optimum tree when a tree has been provided in the nexus file
 
### -h
 
  -h prints a description of the usage of machete.
 
### -s
 
  -s [constraint number] specifies the constraint at which to start the analysis
 
### -e
 
  -e [constraint number] specifies the constraint at which to end the analysis
 
### -l
 
  -l lists constraints (and do not carry out reverse constraints analysis)
  
### -n

  -n tells machete to NOT carry out the reverse constraints analysis. This is useful if you want to use machete to just build ML trees or carry out ML boostrap analyses in PAUP*
  
### -r
 
  -r [INTEGER] turns on the bootstrapping algorithm carrying out the specified number of repetitions (default = 0)
  NOTE: Only the ML tree is used to define the constraints for the likelihood decay analyses, the boostraps are outputted into a sperate file for information only.
  
 
 
 
 
## Outputs

Machete produces three standard output files:

```
[NEXUSFILE].bestMLtree.tre
[NEXUSFILE].likelihood.decays.tre
[NEXUSFILE].sitelike.txt
[NEXUSFILE].constraint.tre
```

*Where* `[NEXUSFILE]` is the name of the alignment input file.

`[NEXUSFILE].bestMLtree.tre` will contain the best ML tree calculated (or the one provided) with branchlengths optimised based onthe the model optimisation carried out.

`[NEXUSFILE].likelihood.decays.tre` will again contain the phylogeny either calculated by PAUP or provided by the user with internal branch labels containg the results of the likelihood decay analysis in the following format: `2/15.19/41.0945/25.9042/0.613363` where:

>2 = the internal branch ID (and constraint number)
  
>15.19 = The difference in the -lnL between the unconstrained tree and constrained tree for this branch.
  
>41.0945 = The absolute value of the sum of likelihood differences supporting the unconstrained tree at this branch.
  
>25.9042 = The absolute value of the sum of the likelihood differences supporting the constrained tree at this branch.
  
>0.613363 = Proportion of the absolute likelihood differences supporting the unconstrained tree {calculated in this case as: (41.0945)/(41.0945+25.9042)}

These supports can be viewed on the tree using a phylogeny viewer such as [figtree](http://tree.bio.ed.ac.uk/software/figtree/).

`[NEXUSFILE].sitelike.txt` contains all the sitelikelihoods calcuated for the unconstrained tree and for each of the constrained trees (labelled by internal branch ID).

`[NEXUSFILE].constraint.tre` contains all the best constrained trees for each internal branch in phylip format (labelled by internal branch ID).

Two other output files are created if the "-b" bootstrapping option is selected:

```
[NEXUSFILE].contree.tre
[NEXUSFILE].boottrees.tre
```

`[NEXUSFILE].contree.tre` will contain the results of a majority-rule consensus (with minor compatible minor components) of the boostrapped trees in nexus format.
`[NEXUSFILE].boottrees.tre` will contain all the bootstrapped trees in nexus format






