import sys


def extract_gene_id(seq_name):
    return seq_name.split('_')[0]


def generate_gene_pairs(fasta_file):
    gene_seqs = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                seq_name = line[1:]
                gene_id = extract_gene_id(seq_name)
                if gene_id not in gene_seqs:
                    gene_seqs[gene_id] = [seq_name]
                else:
                    gene_seqs[gene_id].append(seq_name)

    gene_pairs = []
    for gene_id, seq_names in gene_seqs.items():
        if len(seq_names) > 1:
            pairs = [(seq_names[i], seq_names[j]) for i in range(len(seq_names)) for j in range(i + 1, len(seq_names))]
            gene_pairs.extend(pairs)

    return gene_pairs


def write_gene_pairs_to_file(gene_pairs, output_file):
    with open(output_file, 'w') as f:
        for pair in gene_pairs:
            f.write(pair[0] + ' ' + pair[1] + '\n')


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gene_pairs.py fasta_file output_file")
        sys.exit(1)
    fasta_file = sys.argv[1]
    output_file = sys.argv[2]
    gene_pairs = generate_gene_pairs(fasta_file)
    write_gene_pairs_to_file(gene_pairs, output_file)
