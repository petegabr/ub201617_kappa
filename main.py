from itertools import permutations

from Bio import SeqIO
from Bio.Alphabet import Alphabet
from Bio.HMM import MarkovModel, Trainer
from Bio.Seq import Seq


class StateAlphabet(Alphabet):
    letters = [str(i) for i in range(3)]


class NucleotideAlphabet(Alphabet):
    letters = ['A', 'C', 'T', 'G']


class Kmer2Alphabet(Alphabet):
    nucs = ['A', 'C', 'T', 'G']
    letters = ["".join(p) for p in list(permutations(nucs, 2))]


class Kmer3Alphabet(Alphabet):
    nucs = ['A', 'C', 'T', 'G']
    letters = ["".join(p) for p in list(permutations(nucs, 3))]


def stop_condition(log_likelihood_change, num_iterations):
    """Tell the training model when to stop"""
    if log_likelihood_change < 0.01:
        return 1
    elif num_iterations >= 10:
        return 1
    else:
        return 0


if __name__ == '__main__':
    # Read training data
    with open('data/chr1.genes.train.filtered.fa') as handle:
        train = list(SeqIO.parse(
            handle, format='fasta', alphabet=NucleotideAlphabet()))

    # Read testing data
    with open('data/chr1.genes.test.filtered.fa') as handle:
        test = list(SeqIO.parse(
            handle, format='fasta', alphabet=NucleotideAlphabet()))

    mm_builder = MarkovModel.MarkovModelBuilder(
        StateAlphabet(), NucleotideAlphabet())
    mm_builder.allow_all_transitions()
    mm_builder.set_random_probabilities()

    bw_model = mm_builder.get_markov_model()

    training_seq = Trainer.TrainingSequence(
        train[0].seq, Seq('', StateAlphabet()))

    # list of first n training sequences
    n = 2
    training_seq_multi = [
        Trainer.TrainingSequence(t.seq, Seq('', StateAlphabet())) for t in
        train[0:n]]

    trainer = Trainer.BaumWelchTrainer(bw_model)

    trained_mm = trainer.train([training_seq], stop_condition)
    # trainer on multiple training sequences
    trained_mm_multi = trainer.train(training_seq_multi, stop_condition)

    # vitrebi decoder multi na testnem primeru
    decoded_multi = trained_mm_multi.viterbi(test[0].seq, StateAlphabet())
    # print(decoded_multi)
