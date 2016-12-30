import sys

from collections import deque, defaultdict
from itertools import product

from Bio import SeqIO
from Bio.Alphabet import Alphabet, NucleotideAlphabet
from Bio.HMM import MarkovModel, Trainer
from Bio.Seq import Seq
from collections import Counter

from parse import read_training_data, read_testing_data, read_binding_sites


class BinaryStateAlphabet(Alphabet):
    letters = ['N', 'B']  # N - non binding site, B - binding site


class Kmer1Alphabet(NucleotideAlphabet):
    letters = ['A', 'C', 'G', 'T']


class Kmer2Alphabet(Alphabet):
    nucs = ['A', 'C', 'G', 'T']
    letters = ["".join(p) for p in product(nucs, repeat=2)]


class Kmer3Alphabet(Alphabet):
    nucs = ['A', 'C', 'G', 'T']
    letters = ["".join(p) for p in product(nucs, repeat=3)]


class SeqContent:
    """Calculate some useful statistics of a sequence."""
    def __init__(self, seq):
        self._counter = Counter(seq)
        self.length = sum(self._counter.values())
        self.A = self._counter['A'] / self.length
        self.C = self._counter['C'] / self.length
        self.G = self._counter['G'] / self.length
        self.T = self._counter['T'] / self.length

    def __str__(self):
        return 'A: %d%%\tC: %d%%\tG: %d%%\tT: %d%%\t' % (
            self.A * 100, self.C * 100, self.G * 100, self.T * 100)


def stop_condition(log_likelihood_change, num_iterations):
    """Tell the training model when to stop"""
    if log_likelihood_change < 0.01:
        return 1
    elif num_iterations >= 10:
        return 1
    else:
        return 0

        
def site_contained_within(sequence, bind):
    return sequence.start <= bind.start and sequence.end >= bind.end and sequence.direction == bind.direction

        
def get_binding_sites(binding_sites, data):
    """
    Get a list of binding site sequences from a dataset.

    Parameters
    ----------
    binding_sites : List[TDPBindingSites]
    data : List[SeqRecord]

    Returns
    -------
    List[Seq]

    """
    
    sites = []

    for record in data:
        for site in binding_sites:
            if site_contained_within(record, site):
                # Calculate the offset start and end of the binding site
                offset_start, offset_end = site.start - record.start, site.end - record.start
                sites.append(record.seq[offset_start:offset_end])
                
    return sites


def data_analysis(binding_sites, data):
    # find genes that extend each other (end from one gene == start from another gene)
    if True:
        a, n = [], 0
        for seq in data:
            a.extend([(seq.start, "s", seq.direction), (seq.end, "e", seq.direction)])
        
        a.sort()
        
        for i in range(1, len(a)):
            if a[i][0] == a[i-1][0]:
                #print(a[i], "\t", a[i-1], "\t", a[i][2] == a[i-1][2])
                n += 1
        print("GENES THAT EXTEND EACH OTHER:")
        print("\tnum occurrences:", n)
        # n > 0 iff we put train + test data in
    
    # find how many binds have a corresponding gene
    if True:
        n_found, n_not_found = 0, 0
        for bind in binding_sites:
            for seq in data:
                if site_contained_within(seq, bind):
                    #print("Gene for bind found", "(" + str(bind) + ")")
                    n_found += 1
                    break;
            else:
                #print("Gene for bind NOT FOUND", "(" + str(bind) + ")")
                n_not_found += 1
                
        print("BINDS WITH CORRESPONDING GENES:")
        print("\tFound:", n_found, "Not found:", n_not_found)
    
    
    # extract binding nucleotides
    if True:
        site_seqs = get_binding_sites(binding_sites, data)
        alpha, cum = ["A", "C", "G", "T"], defaultdict(int)
        for ss in site_seqs:
            print(ss)
            ss = SeqContent(ss)
            for a in alpha:
                cum[a] += ss._counter[a]
        
        for a in alpha:
            print(a + ":", cum[a] / sum(cum.values()))
            
    
def make_path(record, bs_data, state_alph):
    """
    Generate state path of emitted record, based on binding sites data and alphabets used

    Parameters
    ----------
    record : Seq ... emitted sequence
    bs_data: List[TDPBindingSites] ... binding sites data
    state_alph: binary state alphabet

    Returns
    -------
    Seq

    """
    
    path = [state_alph.letters[0]] * len(record.seq)
    
    for bind in bs_data:
        if site_contained_within(record, bind):
            s, e = bind.start - record.start, bind.end - record.start
            path[s:e] = [state_alph.letters[1]] * (e - s)

    return Seq("".join(path), state_alph)


def get_training_seq(train_data, bs_data, state_alph):
    """
    Get a list of N training sequences with corresponding state paths

    Parameters
    ----------
    train_data : List[SeqRecord] ... training dataset
    bs_data: List[TDPBindingSites] ... binding sites data
    state_alph: Alphabet ... state alphabet

    Returns
    -------
    List[TrainingSequence]

    """

    training_seqs = []
    for record in train_data:
        path = make_path(record, bs_data, state_alph)
        training_seqs.append(Trainer.TrainingSequence(record.seq, path))
    
    return training_seqs


def train_hmm(train_data, bs_data, state_alph, em_alph, trainer='BW'):
    """
    get trained markov model using BW or Known State Trainer (KST)

    Parameters
    ----------
    train_data: List[Seq] ... training dataset
    bs_data: List[TDPBindingSites] ... binding sites dataset
    state_alph: state alphabet
    em_alph: emission alphabet
    trainer : string (BW ali KTS)

    Returns
    -------
    MarkovModel.HiddenMarkovModel
    """
    
    # example HMM: https://github.com/biopython/biopython/blob/master/Tests/test_HMMCasino.py
    
    mm_builder = MarkovModel.MarkovModelBuilder(state_alph, em_alph)
    mm_builder.allow_all_transitions()
    mm_builder.set_random_probabilities()
    mm_model = mm_builder.get_markov_model()

    #make training sequence with corresponding state paths from first n emissions
    training_seq = get_training_seq(train_data, bs_data, state_alph)
    
    if trainer == 'BW': # ne pride v postev ker imamo znane poti
        bw_trainer = Trainer.BaumWelchTrainer(mm_model)
        trained_bw = bw_trainer.train(training_seq, stop_condition)
        return trained_bw
    elif trainer == 'KST':
        kst_trainer = Trainer.KnownStateTrainer(mm_model)
        trained_kst = kst_trainer.train(training_seq)
        return trained_kst
    else:
        return None


def viterbi_decode(hmm_model, sequences, state_alph=BinaryStateAlphabet):
    """
    return hidden paths for sequences using trained hmm

    Parameters
    ----------
    hmm_model: HiddenMarkovModel
    sequences: List[SeqRecord]

    Returns
    -------
    List[Seq]
    """
    
    decoded_paths = []
    
    for seq in sequences:
        decoded_paths.append(hmm_model.viterbi(seq.seq, state_alph)[0])
    
    return decoded_paths

    
def evaluate_model(hmm_model, test_data, bs_data, state_alph=BinaryStateAlphabet):
    # basic machine learning statistics... not the best evaluation :D
    decoded_paths = viterbi_decode(hmm_model, test_data, state_alph)
    real_paths = []
    
    for seq in test_data:
        real_paths.append(make_path(seq, bs_data, state_alph))
    
    TP, TN, FP, FN = 0, 0, 0, 0
    
    for dp, rp in zip(decoded_paths, real_paths):
        for i in range(len(dp)):
            if dp[i] == "N" and rp[i] == "N": TN += 1
            elif dp[i] == "N" and rp[i] == "B": FN += 1
            elif dp[i] == "B" and rp[i] == "N": FP += 1
            elif dp[i] == "B" and rp[i] == "B": TP += 1
    
    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    #F1 = 2 * precision * recall / (precision + recall)
    
    print("TP:", TP)
    print("TN:", TN)
    print("FP:", FP)
    print("FN:", FN)
    print("precision:", precision)
    print("recall:", recall)
    #print("F1:", F1)
        
if __name__ == '__main__':
    # Read training data
    train = list(read_training_data(Kmer1Alphabet()))
    # Read testing data
    test = list(read_testing_data(Kmer1Alphabet()))
    # Read the binding sites
    binding_sites = list(read_binding_sites())
    
    #data_analysis(binding_sites, train + test)
    
    # train model with these params
    print("training model")
    trained_model = train_hmm(train[:20], binding_sites, BinaryStateAlphabet(), Kmer1Alphabet(), 'KST')
    print("evaluating model")
    evaluate_model(trained_model, train[:100], binding_sites)
    
    #print(trained_model.transition_prob)
    #print(trained_model.emission_prob)
    #paths = viterbi_decode(trained_model, train[:10])
    
    # TODO: evaluacija: primerjava decoded testnih primerov z binding sites podatki
    
    sys.exit(0)