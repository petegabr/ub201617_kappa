from collections import deque
from itertools import permutations

from Bio import SeqIO
from Bio.Alphabet import Alphabet, NucleotideAlphabet
from Bio.HMM import MarkovModel, Trainer
from Bio.Seq import Seq
from collections import Counter

from parse import read_training_data, read_testing_data, read_binding_sites


class StateAlphabet(Alphabet):
    letters = [str(i) for i in range(3)]


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
    """Get a list of binding site sequences from a dataset.

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

if __name__ == '__main__':
    # Read training data
    train = list(read_training_data(Kmer1Alphabet()))
    # Read testing data
    # test = read_testing_data(Kmer1Alphabet())
    # Read the binding sites
    binding_sites = list(read_binding_sites())
    # Show what actually occurs at the binding sites
    site_seqs = get_binding_sites(binding_sites, train[:50])
    for ss in site_seqs:
        print(ss)
        # print(SeqContent(ss))

    import sys
    sys.exit(0)

    mm_builder = MarkovModel.MarkovModelBuilder(
        StateAlphabet(), Kmer1Alphabet())
    mm_builder.allow_all_transitions()
    mm_builder.set_random_probabilities()

    bw_model = mm_builder.get_markov_model()

    trainer = Trainer.BaumWelchTrainer(bw_model)

    # Convert the first gene to a training sequence for the HMM
    # training_seq = Trainer.TrainingSequence(
    #     train[0].seq, Seq('', StateAlphabet()))
    # trained_mm = trainer.train([training_seq], stop_condition)

    # Convert the first `n` genes to training sequnces for the HMM
    # n = 2
    # training_seq_multi = [
    #     Trainer.TrainingSequence(t.seq, Seq('', StateAlphabet())) for t in
    #     train[0:n]]
    # trainer on multiple training sequences
    # trained_mm_multi = trainer.train(training_seq_multi, stop_condition)
    # vitrebi decoder multi na testnem primeru
    # decoded_multi = trained_mm_multi.viterbi(test[0].seq, StateAlphabet())
    # print(decoded_multi)
