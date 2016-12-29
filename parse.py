from Bio import SeqIO

TRAIN_DATA = 'data/chr1.genes.train.filtered.fa'
TEST_DATA = 'data/chr1.genes.test.filtered.fa'
KNOWN_BINDING_SITES_DATA = 'data/iCLIP_TDP-43_tollervey2011_hg19_chr1.bed'


class TDPBindingSites:
    DIRECTION_NEGATIVE, DIRECTION_POSITIVE = range(2)

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
        site_direction = [cls.DIRECTION_NEGATIVE, cls.DIRECTION_POSITIVE][
            int(direction == '+')]
        return cls(chrom, start, end, method, strength, site_direction)

    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}'.format(
            self.chromosome, self.start, self.end, self.strength,
            self.direction)


def read_binding_sites():
    with open(KNOWN_BINDING_SITES_DATA) as handle:
        for line in handle:
            yield TDPBindingSites.from_string(line)


def read_training_data(alphabet=None):
    with open(TRAIN_DATA) as handle:
        yield from SeqIO.parse(handle, format='fasta', alphabet=alphabet)


def read_testing_data(alphabet=None):
    with open(TEST_DATA) as handle:
        yield from SeqIO.parse(handle, format='fasta', alphabet=alphabet)

if __name__ == '__main__':
    print("\n".join(map(str, list(read_binding_sites()))))
