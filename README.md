# machete
Machete: Likelihood reverse constraint analysis using PAUP

##Overview
Machete implements a likelihood reverse constraint method for assessing the support in an alignment for each internal branch of a tree.

Machete takes as input a nexus formatted aligned DNA or Amino Acid sequences and uses PAUP to calculate the reverse constraint supports.

Machete controls and interacts with Paup using a pipe, and not using a predefined script. This allows dataset-specific optimisations to be carried out (as a user would).

##Algorithm

Machete carrys out the following steps:

1. Identifies the phylogeny best supported by the inputted alignment by:
  * Carrying out a heuristic search of tree space.
  * Optimising model parameters using best tree.
  * Repeat previous steps until likelihoods stops improving.

2. Define all constraints based on each of the internal branches of the best tree.

3. For each constraint, carry out a search for the the best tree that does not contain that branch.

5. Calculate the likelihood decay support for each internal branch and output the best phylogeny with the support values labelled on each branch


# Installation

Download the machete.c file above, or if you have git installed, type the command `git clone 

```
cc machete.c -o machete -lm
```




