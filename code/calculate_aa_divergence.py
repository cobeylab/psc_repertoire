#!/usr/bin/python

"""Main function computes amino acid divergence for all pairs of sequences and also relative to a reference sequence
"""
import sys
from Bio import AlignIO
from Bio.Seq import Seq
import os
import itertools
import re
import numpy
import glob


# Function for removing sites with over 95% gaps (hard-coded threshold)
def remove_mostly_gaps(alignment):
    retained_columns = []
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        fraction_gaps = float(column.count('-')) / len(column)
        if fraction_gaps < 0.95:
            retained_columns.append(alignment[:, i:i + 1])
    alignment = retained_columns[0]

    for i in range(1, len(retained_columns)):
        alignment = alignment + retained_columns[i]

    return alignment


def translate_with_gaps(seq):
    aa_seq = ''
    for codon in range(len(seq)/3):
        nt_seq = seq[3*codon:3*codon + 3]
        if nt_seq.find('-') == -1:
            aa = str(Seq(nt_seq).translate())
        else:
            aa = 'X'
        aa_seq = aa_seq + aa
    return aa_seq


def calc_pairwise_divergence(alignment):
    alignment = list(alignment)
    seq_ids = [seq.id for seq in alignment]
    aa_divergence = {}
    seq_pairs = list(itertools.combinations(seq_ids, 2))
    for pair in seq_pairs:
        seq1 = [str(seq.seq) for seq in alignment if seq.id == pair[0]][0]
        seq2 = [str(seq.seq) for seq in alignment if seq.id == pair[1]][0]

        aa_seqs = [translate_with_gaps(seq) for seq in [seq1, seq2]]

        diff_sites = [i for i in range(len(aa_seqs[0])) if aa_seqs[0][i] != aa_seqs[1][i]]
        sites_with_gaps = [i for i in range(len(aa_seqs[0])) if aa_seqs[0][i] == 'X' or aa_seqs[1][i] == 'X']
        # Count only differences in sites where neither sequence has gaps
        diff_sites = [i for i in diff_sites if i not in sites_with_gaps]

        n_aa_diffs = len(diff_sites)
        n_nongap_sites = len(aa_seqs[0]) - len(sites_with_gaps)

        aa_divergence[pair] = float(n_aa_diffs) / n_nongap_sites
    return aa_divergence


def main(argv):

    alignment_path = argv[1]
    reference_seq_id = argv[2]

    alignment = AlignIO.read(alignment_path, 'fasta')
    #alignment = remove_mostly_gaps(alignment)

    #parent_dir = os.path.dirname(alignment_path)
    alignment_name = re.search(r'[^/]+\.fasta', alignment_path).group().replace('_alignment','').replace('.fasta','')

    aa_divergence = calc_pairwise_divergence(alignment)
    aa_divergence_from_ref = [aa_divergence[pair] for pair in aa_divergence.keys() if reference_seq_id in pair]

    n_sites = str(len(translate_with_gaps(str(alignment[0].seq))))
    # Number of sites in ref sequence (excluding gaps)
    ref_seq = [translate_with_gaps(str(seq.seq).replace('-', '')) for seq in alignment if seq.id == reference_seq_id]
    n_sites_ref = str(len(ref_seq[0]))

    csv_header = 'reference,mean_divergence,median_divergence,min_divergence,max_divergence,sd'
    csv_header = csv_header + ',n_sites_in_alignment,n_sites_in_reference_seq\n'

    with open('../results/aa_divergence/' + alignment_name + '_pairwise_divergence.csv', 'w') as output:
        output.write(csv_header)
        mean_divergence = str(numpy.mean(aa_divergence.values()))
        median_divergence = str(numpy.median(aa_divergence.values()))
        min_divergence = str(min(aa_divergence.values()))
        max_divergence = str(max(aa_divergence.values()))
        sd = str(numpy.std(aa_divergence.values()))

        mean_divergence_ref = str(numpy.mean(aa_divergence_from_ref))
        median_divergence_ref = str(numpy.median(aa_divergence_from_ref))
        min_divergence_ref = str(min(aa_divergence_from_ref))
        max_divergence_ref = str(max(aa_divergence_from_ref))
        sd_ref = str(numpy.std(aa_divergence_from_ref))

        output.write(','.join(['all_pairs', mean_divergence, median_divergence, min_divergence, max_divergence, sd, n_sites, n_sites_ref+'\n']))

        output.write(','.join([reference_seq_id, mean_divergence_ref, median_divergence_ref, min_divergence_ref,max_divergence_ref, sd_ref, n_sites, n_sites_ref + '\n']))


if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)



