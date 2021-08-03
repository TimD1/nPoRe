import os, sys, subprocess
import argparse
import numpy as np

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from vcf import *

from cig import expand_cigar, standardize_cigar
import cfg as cfg



def main():

    print(f"> splitting vcf '{cfg.args.vcf}'")
    vcf1, vcf2 = split_vcf(cfg.args.vcf, 'out/vcf')

    print(f"> indexing '{vcf1}' and '{vcf2}'")
    subprocess.run(['tabix', '-p', 'vcf', vcf1])
    subprocess.run(['tabix', '-p', 'vcf', vcf2])

    print(f"> applying '{vcf1}'")
    seq1, cig1 = apply_vcf(vcf1, cfg.args.ref)
    print(f"> applying '{vcf2}'")
    seq2, cig2 = apply_vcf(vcf2, cfg.args.ref)

    print(f"> reading reference")
    ref = get_fasta(cfg.args.ref, cfg.args.contig)

    print(f"> standardizing cigar1")
    cigar1_data = standardize_cigar(("hap1", "chr19", 0, cig1, ref, seq1))
    print(f"\n> standardizing cigar2")
    cigar2_data = standardize_cigar(("hap2", "chr19", 0, cig2, ref, seq2))

    print('\n'

    # print(f"> merging vcfs: '{vcf1}' and '{vcf2}'")
    # merge_vcfs(vcf1, vcf2)



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("--vcf", type=str, default='/x/gm24385/chr19/guppy_4_0_11/clair3_ra/pileup.vcf.gz')
    parser.add_argument("--ref", type=str, default='/x/gm24385/chr19/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta')
    parser.add_argument("--contig", type=str, default="chr19")
    parser.add_argument("--contig_beg", type=int, default=1)
    parser.add_argument("--contig_end", type=int, default=40_000_000)
    parser.add_argument("--indels_only", default=True, action="store_true")
    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
