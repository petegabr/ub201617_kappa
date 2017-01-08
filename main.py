from collections import defaultdict, Counter, Iterable
from itertools import product
from Bio.Alphabet import Alphabet, NucleotideAlphabet, SingleLetterAlphabet
from Bio.HMM import MarkovModel, Trainer
from Bio.Seq import Seq
import sys

from parse import read_training_data, read_testing_data, read_binding_sites
from data.trained_data import *


class SeqContent:
    """Calculate some useful statistics of a sequence."""

    def __init__(self, sequences):
        # Prefer treating everything like a list of sequences
        if not isinstance(sequences, Iterable):
            sequences = (sequences,)
        # Initialize variables
        self._counter = Counter()
        self.length = 0

        for seq in sequences:
            seq_counter = Counter(seq)
            self._counter.update(seq_counter)
            self.length += sum(seq_counter.values())

        self.A = self._counter['A'] / (self.length or 1)
        self.C = self._counter['C'] / (self.length or 1)
        self.G = self._counter['G'] / (self.length or 1)
        self.T = self._counter['T'] / (self.length or 1)

    def __str__(self):
        return 'A: %d%%\tC: %d%%\tG: %d%%\tT: %d%%\t' % (
            self.A * 100, self.C * 100, self.G * 100, self.T * 100)


def extract_binding_sequences(binding_sites, sequences):
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

    for seq in sequences:
        for bind in binding_sites:
            if bind.is_contained_within(seq):
                # Calculate the offset start and end of the binding site
                offset_start = bind.start - seq.start
                offset_end = bind.end - seq.start
                sites.append(seq.seq[offset_start:offset_end])

    return sites


def data_analysis(binding_sites, data):
    """Do some simple data analysis."""
    
    # find genes that extend each other (end from one gene == start from another gene)
    # only happens if data = train + test data
    if False:
        a, n = [], 0
        for seq in data:
            a.extend([(seq.start, "s", seq.direction), (seq.end, "e", seq.direction)])

        a.sort()

        for i in range(1, len(a)):
            if a[i][0] == a[i - 1][0]:
                # print(a[i], "\t", a[i-1], "\t", a[i][2] == a[i-1][2])
                n += 1
        print("GENES THAT EXTEND EACH OTHER:")
        print("\tnum occurrences:", n)
        # n > 0 iff we put train + test data in

    # find how many binds have a corresponding gene
    # all bind have a corresponding gene thanks to Peters filtering
    if False:
        n_found, n_not_found = 0, 0
        for bind in binding_sites:
            for seq in data:
                if bind.is_contained_within(seq):
                    # print("Gene for bind found", "(" + str(bind) + ")")
                    n_found += 1
                    break;
            else:
                # print("Gene for bind NOT FOUND", "(" + str(bind) + ")")
                n_not_found += 1

        print("BINDS WITH CORRESPONDING GENES:")
        print("\tFound:", n_found, "Not found:", n_not_found)

    # extract binding nucleotides content
    if True:
        binding_sequences = extract_binding_sequences(binding_sites, data)
        print("BINDING NUCLEOTIDES CONTENT:")
        print("\t", SeqContent(binding_sequences))


def get_training_seq(sequences, binding_sites):
    """
    Get a list of training sequences with corresponding state paths

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
    for seq in sequences:
        path = ["N"] * len(seq)

        for bind in binding_sites:
            if bind.is_contained_within(seq):
                start, end = bind.start - seq.start, bind.end - seq.start
                path[start:end] = ["B"] * (end - start)
        
        path = Seq("".join(path))
        training_sequence = Trainer.TrainingSequence(seq.seq, path)
        training_seqs.append(training_sequence)

    return training_seqs


def train_model(training_data, binding_sites, mer_len):
    """
    Trains a HMM model from training data and binding sites, using k-mer emissions.

    Parameters
    ----------
    training_data: List[SeqRecord]
    binding_sites: List[TDPBindingSites]
    mer_len: integer

    Returns
    -------
    Dict, Dict
    """

    transition_alphabet = list('BN')
    emission_alphabet = ["".join(p) for p in product('ACGT', repeat=mer_len)]
    transition_probability = {(a, b): 0 for a in transition_alphabet for b in transition_alphabet}
    emission_probability = {(a, b): 0 for a in transition_alphabet for b in emission_alphabet}

    sequences = get_training_seq(training_data, binding_sites)
    for i, sequence in enumerate(sequences):
        print("\r%d/%d" % (i + 1, len(sequences)), end="", flush=True)
        states = str(sequence.states)
        emissions = str(sequence.emissions)
        for i in range(len(states) - mer_len):
            state, next_state = states[i:i + 2]
            transition_probability[(state, next_state)] += 1
            k_mer = emissions[i:i + mer_len]
            emission_probability[(state, k_mer)] += 1
    print()

    # normalize transitions
    n_transitions = list(filter(lambda x: x[0] == "N", transition_probability))
    b_transitions = list(filter(lambda x: x[0] == "B", transition_probability))
    n_count = sum(map(lambda x: transition_probability[x], n_transitions))
    b_count = sum(map(lambda x: transition_probability[x], b_transitions))
    for n, b in zip(n_transitions, b_transitions):
        transition_probability[n] /= n_count
        transition_probability[b] /= b_count

    # normalize emissions
    n_emissions = list(filter(lambda x: x[0] == "N", emission_probability))
    b_emissions = list(filter(lambda x: x[0] == "B", emission_probability))
    n_count = sum(map(lambda x: emission_probability[x], n_emissions))
    b_count = sum(map(lambda x: emission_probability[x], b_emissions))
    for n, b in zip(n_emissions, b_emissions):
        emission_probability[n] /= n_count
        emission_probability[b] /= b_count

    return transition_probability, emission_probability
    

def viterbi_decode(emissions, mer_len, transition_probability, emission_probability):
    """
    Does a Viterbi decode on emissions, using k-mer emissions and previously 
    learned transition and emission probabilities.

    Parameters
    ----------
    emissions: String
    mer_len: integer
    transition_probability: Dict
    emission_probability: Dict

    Returns
    -------
    String
    """

    lin = [{"B": (0, None), "N": (1, None)}]
    for i in range(len(emissions) - mer_len):
        emitted_symbols = emissions[i:i + mer_len]
        l, r = lin[-1]["B"], lin[-1]["N"]
        # B
        eb = emission_probability[("B", emitted_symbols)]
        t_lb = (l[0] * eb * transition_probability[("B", "B")], "B")
        t_rb = (r[0] * eb * transition_probability[("N", "B")], "N")
        a = max(t_lb, t_rb)
        # N
        en = emission_probability[("N", emitted_symbols)]
        t_ln = (l[0] * en * transition_probability[("B", "N")], "B")
        t_rn = (r[0] * en * transition_probability[("N", "N")], "N")
        b = max(t_ln, t_rn)
        if (a[0] < 0.1) and (b[0] < 0.1):  # because of large values
            x, y = a
            a = (x * 10, y)
            x, y = b
            b = (x * 10, y)
        lin.append({"B": a, "N": b})

    lin = lin[::-1]
    decoded_path = []

    current = max(lin[0]["B"], lin[0]["N"])[1]
    for l in lin:
        decoded_path.append(current)
        current = l[current][1]

    return "".join(decoded_path[::-1]) + "".join(["N"] * (mer_len - 1))


def evaluate_model(testing_data, binding_sites, transition_probability, emission_probability, mer_len):
    """Performs model evaluation."""
    
    real_paths = [str(seq.states) for seq in get_training_seq(testing_data, binding_sites)]
    decoded_paths = []
    
    for i, seq in enumerate(testing_data):
        print("\r%d/%d" % (i + 1, len(testing_data)), end="", flush=True)
        decoded_paths.append(viterbi_decode(str(seq.seq), mer_len, transition_probability, emission_probability))
    print()
    
    binds = extract_binding_sequences_from_paths(decoded_paths, [str(seq.seq) for seq in testing_data])
    print("Binds content:")
    print("\t", SeqContent(binds))
    
    TP, TN, FP, FN = 0, 0, 0, 0

    for dp, rp in zip(decoded_paths, real_paths):
        for d, r in zip(dp, rp):
            if d == "N" and r == "N": TN += 1
            elif d == "N" and r == "B": FN += 1
            elif d == "B" and r == "N": FP += 1
            elif d == "B" and r == "B": TP += 1

    precision = TP / ((TP + FP) or 1)
    recall = TP / ((TP + FN) or 1)
    F1 = 2 * precision * recall / ((precision + recall) or 1)
    
    print("TP:", TP)
    print("TN:", TN)
    print("FP:", FP)
    print("FN:", FN)
    print("precision:", precision)
    print("recall:", recall)
    print("F1:", F1)
    print()
    

def extract_binding_sequences_from_paths(paths, sequences):
    """
    Get a list of binding site sequences from a dataset using path as a mask.

    Parameters
    ----------
    paths: List[String]
    sequences: List[String]

    Returns
    -------
    List[String]

    """

    nucleotides = []

    for path, seq in zip(paths, sequences):
        nucl = []
        for p, s in zip(path, seq):
            if p == "B":
                nucl.append(s)
        nucleotides.append("".join(nucl))

    return nucleotides


def get_most_probable_paths(records, num_states):
    """Generate a function that can infer the most probable path of a sequence.

    Seeing as only two states cannot model DNA very well, it can make sense to
    use multiple states. These states are trained using a BW model and a fixed
    number of states. We then use this model on any sequences to infer the
    states for its state path.

    Of course, this doesn't account for any binding site states, so we have to
    include those states on top of these afterwards.

    Parameters
    ----------
    records : SeqRecord
        Training records to train the BW model
    num_states : int
        Number of different states to construct. This is currently limited to
        twice (lowercase and uppercase) the number of letters in the ASCII
        English alphabet.

    Returns
    -------
    List[string]
        Most probably state paths of the input sequences.

    """
    # Create a state alphabet with ASCII letters as states
    class StateAlphabet(SingleLetterAlphabet):
        import string
        letters = list(enumerate(string.ascii_letters[:num_states]))

    alphabet = StateAlphabet()

    # Construct training sequences with unknown state paths for the BW model
    # to infer the best states
    training_seqs = [
        Trainer.TrainingSequence(rec.seq, Seq('', alphabet=alphabet))
        for rec in records]

    # Train the BW model
    mm_builder = MarkovModel.MarkovModelBuilder(
        alphabet, records[0].seq.alphabet)
    mm_builder.allow_all_transitions()
    mm_builder.set_random_probabilities()
    mm_model = mm_builder.get_markov_model()

    def stop_condition(log_likelihood_change, num_iterations):
        """Tell BW when to stop"""
        if log_likelihood_change < 0.01:
            return 1
        elif num_iterations >= 3:
            return 5
        else:
            return 0

    bw_trainer = Trainer.BaumWelchTrainer(mm_model)
    trained_bw = bw_trainer.train(training_seqs, stop_condition)

    # Infer the most probable paths
    def get_most_probable_paths(records):
        return [trained_bw.viterbi(rec.seq, alphabet)[0] for rec in records]
    return get_most_probable_paths


if __name__ == '__main__':
    # Read training data
    train_data = list(read_training_data())
    # Read testing data
    test_data = list(read_testing_data())
    # Read the binding sites
    binding_sites = list(read_binding_sites())

    # STRENGTH_THRESHOLD = 0
    # train = list(filter(lambda x: x.direction == DIRECTION, train))
    # test = list(filter(lambda x: x.direction == DIRECTION, test))
    # binding_sites = list(filter(
    #         lambda x: (x.direction == x.DIRECTION_NEGATIVE) and (int(x.strength) >= STRENGTH_THRESHOLD), binding_sites))
    # binding_sites = list(filter(
    #         lambda x: (x.direction == x.DIRECTION_POSITIVE) and (int(x.strength) >= STRENGTH_THRESHOLD), binding_sites))
    
    # Analyse data
    """
    data_analysis(binding_sites, train_data)
    data_analysis(binding_sites, test_data)
    data_analysis(binding_sites, train_data + test_data)
    """
    
    # TRAINING
    print("training")
    mer_len = 2
    print("mer_len:", mer_len)
    transition_probability, emission_probability = train_model(train_data, binding_sites, mer_len)
    #emission_probability = emission_probability_1
    #print(transition_probability)
    #print(emission_probability)

    # TESTING
    print("testing")
    evaluate_model(test_data, binding_sites, transition_probability, emission_probability, mer_len)

