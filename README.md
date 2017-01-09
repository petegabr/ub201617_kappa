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
--- | --- | --- | ---
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

### Most probable binding motifs

According to our model, these are the most probable binding motifs for different k-mers:

k-mer | 1. most common | 2. most common | 3. most common | 4. most common
--- | --- | --- | --- | ---
1-mer | T 0.42 | G 0.25 | A 0.20 | C 0.13
2-mer | TG 0.15 | TT 0.14 | GT 0.13 | AT 0.08
3-mer | TGT 0.10 | GTG 0.08 | TTT 0.07 | ATG 0.04
4-mer | TGTG 0.07 | GTGT 0.06 | TTTT 0.04 | ATGT 0.02
5-mer | TGTGT 0.05 | GTGTG 0.05 | TTTTT 0.02 | TGTAT 0.01
6-mer | TGTGTG 0.04 | GTGTGT 0.04 |  TTTTTT 0.01 | TGTGTA 0.01

The two most probable motifs for 5-mer "TGTGT" and "GTGTG" are also presented in the article "Tollervey, James R., Curk, Tomaž et al. "Characterizing the RNA targets and
position-dependent splicing regulation by TDP-43." Nature Neuroscience 14.4(2011): 452-458."

### Future Work and Improvements

The models we trained have some limitations.

1. We always trained models using fixed kmers. Perhaps we could obtain better results if we tried mixing the alphabets and different kmer lengths e.g. so we would differentiate the sequences AC from ACG for instance.
2. We used a first order HMM, perhaps using a second or third order Markov model may produce better results.
3. Our model was limited to two states, a Binding and a Non-binding state. This meant that the model had to generalise the non-binding state very broadly, although we know the genome contains several different states (e.g. promoter regions, CpG islands, ...). Our idea was to use a Baum-Welsh model to model the genome to get _n_ best states that would produce the initial state paths, but that turned out to be far to slow. Perhaps using data available on the Internet to get the initial states would produce better results.
4. Different models may produce better results. Perhaps an HMM is too simplistic, and some more complex models could be used e.g. deep neural networks are especially popular at the moment and are producing very good results in most areas, but perhaps even simpler models such as an SVM might produce better results as well.
