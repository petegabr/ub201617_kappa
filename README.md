# ub201617_kappa

PROJECT A1 (research)

Our task was to define hidden and observed variables of a Hidden Markov Model, that models binding sites of RBP, called TDP-43. We were given a training 
and a test dataset and coordinates, where TDP-43 is reported to bind.

We first used some statistics to better understand our data. Our findings were the following:
1. Not all binding sites have a corresponding DNA sequence. About 800 binding sites are not matched, which means they cannot be used in the learning process
2. After calculating the nucleotide content within the known binding sites, it turns out, that binding sites are rich with T and G nucleotides. This is also confirmed
by the articles, listed as reading material in the Project Presentations PDF

Based on the latter, we decided on observed variables. We decided to try and train our model with 3 different sets of observed variables, namely nucleotide sequences
of length 1 (single nucleotide), 2 and 3 (k = 1, 2, 3). 

As hidden variables, we decided to go with binary alphabet, where state N represents a non-binding site and state B represents a binding site.

Each of our records had a strand attribute, which is what we had to take into consideration when training our model. That is why we decided to make a reverse complement
of all the records on negative strand and recalculate binding site start and end coordinates accordingly. 

To train our model we used Known State Trainer, which is already implemented in biopython. For training purposes, we generated a state path for each DNA sequence using 
known binding sites coordinates. 

After training our model on training dataset, we applied it to test dataset, where we used our own implementation of Viterbi Decoding. We calculated some statistics about 
the quality of our model. We calculated precision, recall and F1 score. Below are the results for our model:

k	all		precision	recall		F1 score
---------------------------------------------
1	no		0.0384		0.1837		0.0636
1	yes		0.0192		0.0964		0.0320
2	no		0.2740		0.0831		0.1275
2	yes		0.0796		0.0086		0.0155
3	no		0.3906		0.1116		0.1736
3	yes		0.0508		0.0017		0.0033

When evaluating our model, it turned out that there were a lot of test samples, where our model discovered zero binding sites. Therefore, we calculated our performance
statistics on all results (all = yes) and on the subset, where there were more than zero binding sites discovered (all = no). We can see, that by using larger k, our model 
finds less binding sites, but with higher precision. 

TODO: katerih nukleotidov je v odkritih binding-site-ih največ? Kater je motif, na katerega se veže TDP-43? Future work and improvements?
