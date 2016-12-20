from Bio import SeqIO
from Bio.Alphabet import Alphabet
from Bio.HMM import MarkovModel, Trainer, Utilities
from Bio.Seq import Seq


class StateAlphabet(Alphabet):
    letters = [str(i) for i in range(3)]


class NucleotideAlphabet(Alphabet):
    letters = ['A', 'C', 'T', 'G']


if __name__ == '__main__':
    # Read training data
    with open('data/chr1.genes.train.filtered.fa') as handle:
        train = list(SeqIO.parse(
            handle, format='fasta', alphabet=NucleotideAlphabet()))
    # Read testing data
    # with open('data/chr1.genes.test.filtered.fa') as handle:
    #     test = list(SeqIO.parse(
    #         handle, format='fasta', alphabet=NucleotideAlphabet()))


    mm_builder = MarkovModel.MarkovModelBuilder(
        StateAlphabet(), NucleotideAlphabet())
    mm_builder.allow_all_transitions()
    mm_builder.set_random_probabilities()

    bw_model = mm_builder.get_markov_model()


    def stop_training(log_likelihood_change, num_iterations):
        """Tell the training model when to stop"""
        if log_likelihood_change < 0.01:
            return 1
        elif num_iterations >= 10:
            return 1
        else:
            return 0

    training_seq = Trainer.TrainingSequence(
        train[0].seq, Seq('', StateAlphabet()))
    trainer = Trainer.BaumWelchTrainer(bw_model)
    trained_mm = trainer.train([training_seq], stop_training)
