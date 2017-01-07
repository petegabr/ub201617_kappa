from parse import read_data2


def get_binding_sites_between(filename, a, b):
    binding_sites = open(filename, "r").read().split("\n")[:-1]
    binding_sites = map(lambda x: x.split("\t"), binding_sites)
    binding_sites = map(lambda x: [x[0], int(x[1]), int(x[2])] + x[3:], binding_sites)
    binding_sites = list(filter(lambda x: (x[1] > a) and (x[2] < b), binding_sites))
    for site in binding_sites:
        print(site)


# open original and converted train data
original_train_filename = "../data/chr1.genes.train.filtered.fa"
converted_train_filename = "../data/new_chr1.genes.train.filtered.fa"

# open binding sites
original_binding_sites_filename = "../data/iCLIP_TDP-43_tollervey2011_hg19_chr1.bed"
converted_binding_sites_filename = "../data/binding_sites.bed"

# get one original arbitrary gene on negative strand
original_genes = read_data2(original_train_filename)
original_arbitrary_gene, strand, start, end = None, None, None, None
for gene in original_genes:
    _, strand, start, end = gene.id.split(",")
    if strand == "-1":
        original_arbitrary_gene = gene
        print(gene)
        break

# find associated binding_sites
get_binding_sites_between(original_binding_sites_filename, int(start), int(end))

print()
print()

# find associated converted gene
converted_genes = read_data2(converted_train_filename)
for gene in converted_genes:
    _, converted_strand, converted_start, converted_end = gene.id.split(",")
    if (converted_strand == strand) and (converted_start == start) and (converted_end == end):
        print(gene)
        break

get_binding_sites_between(converted_binding_sites_filename, int(converted_start), int(converted_end))
