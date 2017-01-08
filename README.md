# ub201617_kappa

## PROJECT A1 (research)

Our task was to define hidden and observed variables of a Hidden Markov Model, that models binding sites of RBP, called TDP-43. We were given a training 
and a test dataset and coordinates, where TDP-43 is reported to bind.

### Data

Our data consisted of genes and binding sites of RNA-binding proteins on given genes. Genes were then separated on training set and testing set.

We first used some statistics to better understand our data. Our findings were the following:
1. Not all binding sites have a corresponding DNA sequence. About 800 binding sites are not matched, which means they cannot be used in the learning process
2. After calculating the nucleotide content within the known binding sites, it turns out, that binding sites are rich with T and G nucleotides. This is also confirmed
by the articles, listed as reading material in the Project Presentations PDF

Based on the latter, we decided on observed variables. We decided to try and train our model with 3 different sets of observed variables, namely nucleotide sequences
of length 1 (single nucleotide), 2, 3 and 4 (k = 1, 2, 3, 4). 

Each of our records had a strand attribute, which is what we had to take into consideration when training our model. That is why we decided to make a reverse complement
of all the records on negative strand and recalculate binding site start and end coordinates accordingly.

---

### Hidden Markov Model

As hidden variables, we decided to go with binary alphabet, where state N represents a non-binding site and state B represents a binding site.

---

### Training

To train our model we used Known State Trainer, which is already implemented in biopython. For training purposes, we generated a state path for each DNA sequence using
known binding sites coordinates. With Known State Trainer in biopython we could not analyze k-mers (k > 1), so we implemented our own Known State Trainer, which returned
almost the same results for 1-mers (difference of 0.5%). Our implementation was also faster and we were able to use arbitrary k-mers.

---

### Testing and Evaluating

After training our model on training dataset, we applied it to test dataset, where we used our own implementation of Viterbi Decoding. We calculated some statistics about 
the quality of our model. We calculated precision, recall and F1 score. For TP, TN, FP, FN we took the following:

* True positive (TP): model recognised binding site as binding site (B -> B)
* True negative (TN): model recognised non-binding site as non-binding site (N -> N)
* False positive (FP): model recognised non-binding site as binding site (N -> B)
* False negative (FN): model recognised binding site as non-binding site (B -> N)

Below are the results for our model:

k | Binds content | TP | FP
--- | --- | ---
1 | A: 17% C: 12%  G: 20%  T: 49%	| 12128 | 622951
2 | A: 14% C: 9%   G: 27%  T: 48%	| 1060 | 11556
3 | A: 15% C: 9%   G: 22%  T: 52%	| 216| 2976
4 | A: 15% C: 9%   G: 23%  T: 50%	| 129 | 2142
5 | A: 15% C: 10%  G: 23%  T: 50%	| 113| 1691
6 | A: 14% C: 11%  G: 25%  T: 48%   | 100 | 1424
9 | A: 15% C: 11%  G: 26%  T: 46%   | 3 | 432
10 | A: 16% C: 9%   G: 23%  T: 50%   | 0 | 171

![alt text](https://github.com/petergabrovsek/ub201617_kappa/blob/master/Figures/k-F1.png "Relation between k and F1")

When evaluating our model, it turned out that there were a lot of test samples, where our model discovered zero binding sites. Therefore, we calculated our performance
statistics on all results (all = yes) and on the subset, where there were more than zero binding sites discovered (all = no). We can see, that by using larger k, our model 
finds less binding sites, but with higher precision. Overall, it seems like k=3 is the best choice for our observed variables length.

Consistently with the articles and our findings on training dataset, the discovered binding sites are rich with T and G nucleotides. We can, therefore, conclude, that TDP-43
prefers binding to TG rich sequences.

### FUTURE WORK AND IMPROVEMENTS
