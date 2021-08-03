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

    print(f"> merging vcfs: '{vcf1}' and '{vcf2}'")
    merge_vcfs(vcf1, vcf2)



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("--vcf", type=str, default='/x/gm24385/chr19/guppy_4_0_11/clair3_ra/pileup.vcf.gz')
    parser.add_argument("--contig", type=str, default="chr19")
    parser.add_argument("--contig_beg", type=int, default=1)
    parser.add_argument("--contig_end", type=int, default=40_000_000)
    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
