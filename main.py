import sys

from collections import deque
from itertools import permutations

from Bio import SeqIO
from Bio.Alphabet import Alphabet, NucleotideAlphabet
from Bio.HMM import MarkovModel, Trainer
from Bio.Seq import Seq
from collections import Counter

from parse import read_training_data, read_testing_data, read_binding_sites


class BinaryStateAlphabet(Alphabet):
    letters = ['N', 'B']  # N - non binding site, B - binding site


class Kmer1Alphabet(NucleotideAlphabet):
    letters = ['A', 'C', 'T', 'G']


class Kmer2Alphabet(Alphabet):
    nucs = ['A', 'C', 'T', 'G']
    letters = ["".join(p) for p in list(permutations(nucs, 2))]


class Kmer3Alphabet(Alphabet):
    nucs = ['A', 'C', 'T', 'G']
    letters = ["".join(p) for p in list(permutations(nucs, 3))]


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


def get_binding_sites(binding_sites, data):
    """
    Get a list of binding site sequences from a dataset.

    Parameters
    ----------
    binding_sites : List[TDPBindingSites]
    data : List[Seq]

    Returns
    -------
    List[Seq]

    """
    sites = deque()

    def site_contained_within(rec_start, rec_end, binding_start, binding_end):
        return rec_start <= binding_start and rec_end >= binding_end

    for record in data:
        *_, start, end = record.id.split(',')
        start, end = int(start), int(end)
        for site in binding_sites:
            if site_contained_within(start, end, site.start, site.end):
                # Calculate the offset start and end of the binding site
                offset_start, offset_end = site.start - start, site.end - start
                # We take the complement if the direction is negative
                if site.direction == site.DIRECTION_POSITIVE:
                    sites.append(record.seq[offset_start:offset_end + 1])
                else:
                    sites.append(record.seq.complement()[offset_start:offset_end + 1])

    return list(sites)


def get_training_seq(data, strd, n, alphabet, trainer='BW', mm_model=None):
    """
    Get a list of N training sequences from strand (-1 or 1) for BW model training. If trainer used is KST - construct
    state path using viterbi decoding (required)

    Parameters
    ----------
    data : List[Seq]
    strand: int
    n : int
    alphabet: state alphabet
    trainer: string (which trainer is used: BW or KST)
    mm_model: model for viterbi decoding for path construction

    Returns
    -------
    List[TrainingSequence]

    """

    training_seq_multi = []
    for t in data[0:n]:
        chrom, strnd, start, end = \
            t.id.strip().split(',')
        if int(strnd) == strd:
            if trainer == 'BW':
                training_seq_multi.append(Trainer.TrainingSequence(t.seq, Seq('', alphabet)))
            # TODO: else: kst training sequence with known state path
    return training_seq_multi


def train_mm(state_alph, emm_alph, strnd, n, trainer='BW'):
    """
    get trained markov model using BW or Known State Trainer (KST)

    Parameters
    ----------
    state_alph: state alphabet
    emm_alph: emission alphabet
    strnd: int (na katerem strandu ucimo model)
    n: int (stevilo training nizov za ucenje)
    trainer : string (BW ali KTS)

    Returns
    -------
    MarkovModel.HiddenMarkovModel
    """

    mm_builder = MarkovModel.MarkovModelBuilder(
        state_alph, emm_alph)
    mm_builder.allow_all_transitions()
    mm_builder.set_random_probabilities()
    mm_model = mm_builder.get_markov_model()

    training_seq = get_training_seq(train, strnd, n, state_alph, trainer, mm_model)

    if trainer == 'BW':
        bw_trainer = Trainer.BaumWelchTrainer(mm_model)
        trained_bw = bw_trainer.train(training_seq, stop_condition)
        return trained_bw
    else:
        # TODO
        """
        kst_trainer = Trainer.KnownStateTrainer(mm_model)
        trained_kst = kst_trainer.train(training_seq)
        return trained_kst
        """
        return None

if __name__ == '__main__':
    # Read training data
    train = list(read_training_data(Kmer1Alphabet()))
    # Read testing data
    test = list(read_testing_data(Kmer1Alphabet()))

    # Read the binding sites
    binding_sites = list(read_binding_sites())

    """
    # Show what actually occurs at the binding sites
    site_seqs = get_binding_sites(binding_sites, train[:50])
    for ss in site_seqs:
        print(ss)
        print(SeqContent(ss))
    """

    # sys.exit(0)

    # train model with these params
    trained_model = train_mm(BinaryStateAlphabet, Kmer1Alphabet, -1, 2, 'BW')

    """
    # vitrebi decoder multi na testnem primeru
    decoded = trained_model.viterbi(test[0].seq, BinaryStateAlphabet)
    print(decoded)
    """