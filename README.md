# machete
Machete: Likelihood reverse constraint analysis using PAUP

##Overview
Machete implements a likelihood reverse constraint method for assessing the support in an alignment for each internal branch of a tree.

Machete takes as input a nexus formatted aligned DNA or Amino Acid sequences and uses PAUP to calculate the reverse constraint supports.

Machete controls and interacts with Paup using a pipe, and not using a predefined script. This allows dataset-specific optimisations to be carried out (as a user would).

## Algorithm

Machete carries out the following steps:

1. Identifies the phylogeny best supported by the inputted alignment by:
  * Carrying out a heuristic search of tree space.
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

To use machete, it is necessary to pas as input a nexus formatted alignment file using the option '-f':

```
machete -f NEXUSFILE
```

An example file `Primate.nex` has been provided. To run the reverse constraint analysis for this file use the command:

```
machete -f Primate.nex
```

Other options are:

 -c which tells Machete to print all the commands sent to Paup to the standard error. This can then be redirected to a sperate file using the following syntax:
 
 ```
 machete -f Primate.nex -c 2> paupcommands.txt
 ```
 
  -h prints a description of the usage of machete.
 
 


