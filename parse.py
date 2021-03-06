from Bio import SeqIO
from Bio import Alphabet

# TRAIN_DATA = 'data/chr1.genes.train.filtered.fa'
# TEST_DATA = 'data/chr1.genes.test.filtered.fa'
# KNOWN_BINDING_SITES_DATA = 'data/iCLIP_TDP-43_tollervey2011_hg19_chr1.bed'
TRAIN_DATA = 'data/new_chr1.genes.train.filtered.fa'
TEST_DATA = 'data/new_chr1.genes.test.filtered.fa'
KNOWN_BINDING_SITES_DATA = 'data/binding_sites.bed'

class TDPBindingSite:
    DIRECTION_NEGATIVE, DIRECTION_POSITIVE = -1, 1

    def __init__(self, chromosome, start, end, method, strength, direction):
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.method = method
        self.strength = strength
        self.direction = direction

    @classmethod
    def from_string(cls, string):
        chrom, start, end, method, strength, direction = \
            string.strip().split('\t')
        site_direction = [cls.DIRECTION_NEGATIVE, cls.DIRECTION_POSITIVE][int(direction == '+')]
        return cls(chrom, start, end, method, strength, site_direction)

    def is_contained_within(self, gene_record):
        """Check if the binding site is contained within a gene sequence."""
        return (gene_record.start <= self.start and
                gene_record.end >= self.end and
                gene_record.direction == self.direction)

    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}'.format(
            self.chromosome, self.start, self.end, self.strength,
            self.direction)


def read_binding_sites():
    with open(KNOWN_BINDING_SITES_DATA) as handle:
        for line in handle:
            yield TDPBindingSite.from_string(line)


def read_data(filename):
    with open(filename) as handle:
        for seq in SeqIO.parse(handle, format='fasta', alphabet=Alphabet.DNAAlphabet()):
            _, direction, start, end = seq.id.split(",")
            direction, start, end = int(direction), int(start), int(end)
            seq.start, seq.end, seq.direction = start, end, direction
            yield seq


def read_training_data():
    yield from read_data(TRAIN_DATA)


def read_testing_data():
    yield from read_data(TEST_DATA)
