from collections import defaultdict, Counter, Iterable
from itertools import product
from Bio.Alphabet import Alphabet, NucleotideAlphabet, SingleLetterAlphabet
from Bio.HMM import MarkovModel, Trainer
from Bio.Seq import Seq
import sys

from parse import read_training_data, read_testing_data, read_binding_sites
from data.trained_data import *

class BinaryStateAlphabet(Alphabet):
    letters = ['N', 'B']  # N - non binding site, B - binding site


class Kmer1Alphabet(NucleotideAlphabet):
    letters = ['A', 'C', 'G', 'T']


class Kmer2Alphabet(Alphabet):
    # 16 states
    letters = ["".join(p) for p in product('ACGT', repeat=2)]


class Kmer3Alphabet(Alphabet):
    # 64 states
    letters = ["".join(p) for p in product('ACGT', repeat=3)]


class SeqContent:
    """Calculate some useful statistics of a sequence."""

    def __init__(self, seqs):
        # Prefer treating everything like a list of sequences
        if not isinstance(seqs, Iterable):
            seqs = (seqs,)
        # Initialize variables
        self._counter = Counter()
        self.length = 0

        for seq in seqs:
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


def get_binding_sites(binding_sites, data):
    """
    Get a list of binding site sequences from a dataset.

    In case the direction of the binding site is negative i.e. the binding site
    occurs on the complementary strand, we take the complement of the gene
    sequence. We DO NOT take the reverse complement since this also reverses
    the sequence, which invalidate the start and ending coordinates of the
    binding site.

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
        for binding_site in binding_sites:
            if binding_site.is_contained_within(record):
                # Calculate the offset start and end of the binding site
                offset_start = binding_site.start - record.start
                offset_end = binding_site.end - record.start
                sites.append(record.seq[offset_start:offset_end])

    return sites


def get_gene_sequence_in(records, start, end):
    """Get the gene sequence within a start and end range."""
    return [record.seq[(start - record.start):(end - record.start)]
            for record in records
            if record.start <= start and record.end >= end]


def data_analysis(binding_sites, data):
    # find genes that extend each other (end from one gene == start from another gene)
    if True:
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
    if True:
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
    Seq, Seq (Seq with path, Seq with (new padded) record

    """

    # EMISSION IN PATH LEN MORATA BITI ENAKA !!!
    # vir: http://biopython.org/DIST/docs/api/Bio.HMM.Trainer-pysrc.html

    def most_common(lst):
        # return most common hidden state.
        n = lst.count('N')
        b = lst.count('B')
        if b >= n:
            return 'B'
        elif n > b:
            return 'N'

    def chunks(l, n):
        # split list l in n long chunks
        return [most_common(l[i:i + n]) for i in range(0, len(l), n)]

    emission = record.seq
    em_alph = emission.alphabet

    k = len(em_alph.letters[0])  # k of k-mers
    if len(emission) % k != 0:
        # pad emitted sequence with r nucleotides (A) for len%k == 0
        r = k - (len(emission) % k)
        emission += ''.join(['A'] * r)

    # new padded emission
    record.seq = emission

    # make 1-mer path of all Ns
    path = [state_alph.letters[0]] * len(emission)

    for bind in bs_data:
        if (bind.start >= record.start) and (bind.end <= record.end) and (bind.direction == record.direction):
            s, e = bind.start - record.start, bind.end - record.start
            path[s:e] = [state_alph.letters[1]] * (e - s)

    if k > 1:
        # shorten path by joining k adjacent states - take the most common hidden state of chunk
        path = chunks(path, k)

    return Seq("".join(path), state_alph), record


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
        path, record = make_path(record, bs_data, state_alph)
        training_seqs.append(Trainer.TrainingSequence(record.seq, path))

    return training_seqs


def train_hmm(train_data, bs_data, state_alph, em_alph, trainer='BW'):
    """
    Get trained markov model using BW or Known State Trainer (KST)

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

    See Also
    --------
    https://github.com/biopython/biopython/blob/master/Tests/test_HMMCasino.py
    """

    def stop_condition(log_likelihood_change, num_iterations):
        """Tell the training model when to stop"""
        if log_likelihood_change < 0.01:
            return 1
        elif num_iterations >= 10:
            return 1
        else:
            return 0

    mm_builder = MarkovModel.MarkovModelBuilder(state_alph, em_alph)
    mm_builder.allow_all_transitions()
    mm_builder.set_random_probabilities()
    mm_model = mm_builder.get_markov_model()

    # make training sequence with corresponding state paths from first n emissions
    training_seq = get_training_seq(train_data, bs_data, state_alph)

    #training_seq for 2-mer
    if len(em_alph.letters[0]) == 2:
        for seq in training_seq:
            seq.emissions = list(map(lambda t: "".join(t), zip(seq.emissions, seq.emissions[1:])))
            seq.states = seq.states[:-1]
    #training_seq for 3-mer
    elif len(em_alph.letters[0])==3:
        for seq in training_seq:
            seq.emissions = list(map(lambda t: "".join(t),zip(map(lambda t: "".join(t), zip(seq.emissions, seq.emissions[1:])), seq.emissions[2:])))
            seq.states = seq.states[:-2]

    if trainer == 'BW':  # ne pride v postev ker imamo znane poti
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
    real_paths = [make_path(seq, bs_data, state_alph)[0] for seq in test_data]   #make path returns 2 Seqs: path and new record

    TP, TN, FP, FN = 0, 0, 0, 0

    for dp, rp in zip(decoded_paths, real_paths):
        for i in range(len(dp)):
            if dp[i] == "N" and rp[i] == "N":
                TN += 1
            elif dp[i] == "N" and rp[i] == "B":
                FN += 1
            elif dp[i] == "B" and rp[i] == "N":
                FP += 1
            elif dp[i] == "B" and rp[i] == "B":
                TP += 1

    precision = TP / ((TP + FP) or 1)
    recall = TP / ((TP + FN) or 1)
    F1 = 2 * precision * recall / ((precision + recall) or 1)

    # print("TP:", TP)
    # print("TN:", TN)
    # print("FP:", FP)
    # print("FN:", FN)
    # print("precision:", precision)
    # print("recall:", recall)
    # print("F1:", F1)
    print(precision, recall, F1, sep=",")


def train_model_new(training_data, binding_sites, mer_len):
    transition_alphabet = list('BN')
    emission_alphabet = ["".join(p) for p in product('ACGT', repeat=mer_len)]
    transition_probability = {(a, b): 0 for a in transition_alphabet for b in transition_alphabet}
    emission_probability = {(a, b): 0 for a in transition_alphabet for b in emission_alphabet}

    sequences = get_training_seq(training_data, binding_sites, BinaryStateAlphabet())
    for i, sequence in enumerate(sequences):
        print("\r%d" % i, end="", flush=True)
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


def viterbi_decode_new(emissions, mer_len, transition_probability, emission_probability):
    lin = [{"B": (0.5, None), "N": (0.5, None)}]
    for i in range(len(emissions) - mer_len + 1):
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

    return "".join(decoded_path[::-1])


def evaluate_paths(path, decoded_path):
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    for i in range(min(len(path), len(decoded_path))):
        if (path[i] == "B") and (decoded_path[i] == "B"):
            TP += 1
        elif (path[i] == "N") and (decoded_path[i] == "N"):
            TN += 1
        elif (path[i] == "N") and (decoded_path[i] == "B"):
            FP += 1
        elif (path[i] == "B") and (decoded_path[i] == "N"):
            FN += 1

    precision = TP / ((TP + FP) or 1)
    recall = TP / ((TP + FN) or 1)
    F1 = 2 * precision * recall / ((precision + recall) or 1)

    # return TP, TN, FP, FN

    # print(precision, recall, F1, sep=",")
    # if F1 > 0:
    #     print(TP, TN, FP, FN, sep=",")
    # return [TP, TN, FP, FN]
    if F1 > 0:
        return [TP, TN, FP, FN]
    else:
        return [0, 0, 0, 0]

def decoded_stats(emissions, decoded_path):
    A_B = 0
    C_B = 0
    T_B = 0
    G_B = 0
    A_N= 0
    C_N = 0
    T_N = 0
    G_N = 0
    b_len = 0
    n_len = 0
    min_len = min(len(emissions), len(decoded_path))
    for i in range(min(len(emissions), len(decoded_path))):
        if decoded_path[i] == "B":
            b_len += 1
            if emissions[i] == "A":
                A_B += 1
            elif emissions[i] == "C":
                C_B += 1
            elif emissions[i] == "T":
                T_B += 1
            elif emissions[i] == "G":
                G_B += 1
        else:
            n_len += 1
            if emissions[i] == "A":
                A_N += 1
            elif emissions[i] == "C":
                C_N += 1
            elif emissions[i] == "T":
                T_N += 1
            elif emissions[i] == "G":
                G_N += 1

    if b_len == 0:
        b_len = 1

    if n_len == 0:
        n_len = 1

    ab = (A_B/b_len) * 100
    cb = (C_B/b_len) * 100
    gb = (G_B/b_len) * 100
    tb = (T_B/b_len) * 100
    an = (A_N/n_len) * 100
    cn = (C_N/n_len) * 100
    gn = (G_N/n_len) * 100
    tn = (T_N/n_len) * 100

    return [ab, cb, gb, tb], [an, cn, gn, tn]


def test_model_new(testing_data, binding_sites, transition_probability, emission_probability, mer_len):
    x = [0, 0, 0, 0]
    stats_B = [0, 0, 0, 0]
    stats_N = [0, 0, 0, 0]
    sequences = get_training_seq(testing_data, binding_sites, BinaryStateAlphabet())

    num = 0.0
    for i, sequence in enumerate(sequences):
        num += 1
        states = str(sequence.states)
        emissions = str(sequence.emissions)
        decoded_states = viterbi_decode_new(emissions, mer_len, transition_probability, emission_probability)
        y = evaluate_paths(states, decoded_states)
        stats_Bi, stats_Ni = decoded_stats(emissions, decoded_states)
        for i in range(4):
            x[i] += y[i]
            stats_B[i] += stats_Bi[i]
            stats_N[i] += stats_Ni[i]
        # print("\r%d" % i, end="", flush=True)
        # print(states.count("B"), decoded_states.count("B"))

    for i in range(4):
        stats_B[i] = stats_B[i]/num
        stats_N[i] = stats_N[i]/num

    strB = 'Binding stats: A: %d%%\tC: %d%%\tG: %d%%\tT: %d%%\t' % (
        stats_B[0], stats_B[1], stats_B[2], stats_B[3])
    strN = 'Nonbinding stats: A: %d%%\tC: %d%%\tG: %d%%\tT: %d%%\t' % (
        stats_N[0], stats_N[1], stats_N[2], stats_N[3])
    print(strB)
    print(strN)

    precision = x[0] / (x[0] + x[2])
    recall = x[0] / (x[0] + x[3])
    F1 = 2 * precision * recall / (precision + recall)
    print(precision, recall, F1)



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
    ALPHABET = Kmer3Alphabet()
    STRENGTH_THRESHOLD = 0
    # DIRECTION = 1  # change also filtering binding_sites
    # Read training data
    train = list(read_training_data(Kmer1Alphabet()))
    # Read testing data
    test = list(read_testing_data(Kmer1Alphabet()))
    # Read the binding sites
    binding_sites = list(read_binding_sites())

    # train = list(filter(lambda x: x.direction == DIRECTION, train))
    # test = list(filter(lambda x: x.direction == DIRECTION, test))
    # binding_sites = list(filter(
    #         lambda x: (x.direction == x.DIRECTION_NEGATIVE) and (int(x.strength) >= STRENGTH_THRESHOLD), binding_sites))
    # binding_sites = list(filter(
    #         lambda x: (x.direction == x.DIRECTION_POSITIVE) and (int(x.strength) >= STRENGTH_THRESHOLD), binding_sites))

    # data_analysis(binding_sites, train + test)

    # Until we know how to properly parse the negative strand, we'll only use
    # the positive one since the other one does not contain the correct seqs.
    # binding_sites = [x for x in binding_sites
    #                 if x.direction == x.DIRECTION_POSITIVE]

    # In case we're interested in the nucleotide content within the binding
    # sites
    # print(SeqContent(get_binding_sites(binding_sites, train)))

    # TRAINING
    mer_len = 4
    print("training")
    # transition_probability, emission_probability = train_model_new(train, binding_sites, mer_len)
    emission_probability = emission_probability_4
    print(transition_probability)
    print(emission_probability)


    # TESTING
    print("testing")
    test_model_new(test, binding_sites, transition_probability, emission_probability, mer_len)

    # train model with these params
    # print("training model")
    # trained_model = train_hmm(train[:50], binding_sites, BinaryStateAlphabet(),
    #                           ALPHABET, 'KST')
    # print(trained_model.transition_prob)
    # print(trained_model.emission_prob)

    # print("evaluating model")
    #
    # for i, t in enumerate(test):
    #     print(i, end=",")
    #     evaluate_model(trained_model, [t], binding_sites)
    # paths = viterbi_decode(trained_model, train[:10])
    # print(paths)

