import os, sys
import argparse
import numpy as np

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import pysam

from vcf import *


def main(args):

    # vcf = pysam.VariantFile(args.vcf, 'rb')
    # for record in vcf.fetch('chr19', 40_000_000, 58_592_616):
    #     print(record.pos, ':', record.alleles)

    print('> reading reference')
    ref = get_fasta(cfg.args.ref, cfg.args.ctg)

    print('> splitting reference VCF')
    vcf1a, vcf1b = split_vcf(cfg.args.vcf1, cfg.args.out)
    print('> indexing reference VCF')
    subprocess.run(['tabix', '-p', 'vcf', vcf1a])
    subprocess.run(['tabix', '-p', 'vcf', vcf1b])
    print('> reference VCF to hap seqs')
    seq1a, cig1a = apply_vcf(vcf1a, ref)
    seq1b, cig1b = apply_vcf(vcf1b, ref)
    print(f'lengths: {len(seq1a)},{len(seq1b)}')

    print('> splitting standardized VCF')
    vcf2a, vcf2b = split_vcf(cfg.args.vcf2, cfg.args.out)
    print('> indexing standardized VCF')
    subprocess.run(['tabix', '-p', 'vcf', vcf2a])
    subprocess.run(['tabix', '-p', 'vcf', vcf2b])
    print('> standardized VCF to hap seqs')
    seq2a, cig2a = apply_vcf(vcf2a, ref)
    seq2b, cig2b = apply_vcf(vcf2b, ref)
    print(f'lengths: {len(seq2a)},{len(seq2b)}')

    count = 0
    for base1a, base2a in zip(seq1a, seq2a):
        count += 1
        if base1a != base2a:
            print(f"ERROR: difference at position {count}")
            break




def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("--vcf1", type=str, default="/x/gm24385/chr19/ref/ra_hap_gstd/ref/unified.vcf.gz")
    parser.add_argument("--vcf2", type=str, default="/x/gm24385/chr19/guppy_4_0_11/ra_hap_gstd/ref/std.vcf.gz")

    parser.add_argument("--ref", type=str, default="/x/gm24385/chr19/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta")
    parser.add_argument("--ctg", type=str, default="chr19")

    parser.add_argument("--out", type=str, default="out/vcf")

    return parser



if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main(args)
