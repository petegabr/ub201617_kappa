from Bio import SeqIO


def reverse_complement(s):
    r = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(list(map(lambda x: r[x], list(s)))[::-1])


def split_60(s):
    r = ""
    while len(s) > 60:
        r += s[:60] + "\n"
        s = s[60:]
    if len(s) > 0:
        r += s + "\n"
    return r


b_s = open("iCLIP_TDP-43_tollervey2011_hg19_chr1.bed", "r").read()
binding_sites = b_s.split("\n")[:-1]
binding_sites = list(map(lambda x: x.split("\t"), binding_sites))
binding_sites = list(map(lambda x: [x[0], int(x[1]), int(x[2])] + x[3:], binding_sites))

binding_sites_corrected = open("binding_sites.bed", "w")


def save_binding_sites(begin, end, direction):
    d = "+" if direction == 1 else "-"
    binding_sites_contained = filter(lambda x: (x[1] > begin) and (x[2] < end) and (x[-1] == d), binding_sites)

    if direction == 1:
        for line in binding_sites_contained:
            b_start, b_end = line[1:3]
            l = "chr1\t%s\t%s\tiCLIP#TDP-43_tollervey2011_hg19*TDP-43\t%s\t%s\n" % (
                str(b_start), str(b_end), line[-2], line[-1])
            binding_sites_corrected.write(l)
    else:
        for line in binding_sites_contained:
            b_start, b_end = line[1:3]
            print(begin, end, b_start, b_end, end=" ")
            b_len = b_end - b_start
            b_start = begin + end - b_end
            b_end = b_start + b_len
            print(b_start, b_end)
            # chr1	790704	790718	iCLIP#TDP-43_tollervey2011_hg19*TDP-43	6	+
            l = "chr1\t%s\t%s\tiCLIP#TDP-43_tollervey2011_hg19*TDP-43\t%s\t%s\n" % (
                str(b_start), str(b_end), line[-2], line[-1])
            binding_sites_corrected.write(l)
        # shift binding site positions
        pass


def convert_data(filename):
    f = open("new_" + filename, "w")

    with open(filename) as handle:
        for seq in SeqIO.parse(handle, format='fasta'):
            _, direction, start, end = seq.id.split(",")
            direction, start, end = int(direction), int(start), int(end)
            s = str(seq.seq)
            if direction == -1:
                s = reverse_complement(s)
            s = split_60(s)
            f.write("> chr1,")
            f.write(str(direction) + ",")
            f.write(str(start) + ",")
            f.write(str(end) + "\n")
            f.write(s)

            save_binding_sites(start, end, direction)


filenames = ["chr1.genes.test.filtered.fa", "chr1.genes.train.filtered.fa"]
for filename in filenames:
    convert_data(filename)

binding_sites_corrected.close()
