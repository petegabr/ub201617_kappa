from Bio import SeqIO

if __name__ == '__main__':
    # Read training data
    with open('data/chr1.genes.train.filtered.fa') as handle:
        train = list(SeqIO.parse(handle, format='fasta'))
    # Read testing data
    with open('data/chr1.genes.test.filtered.fa') as handle:
        test = list(SeqIO.parse(handle, format='fasta'))
