import os, sys, subprocess
import pickle
import argparse
import numpy as np

from vcf import *
from cig import *
from bam import *
import cfg as cfg

def main():

    print(f"> splitting SUB and INDEL VCFs into haploid VCFs")
    sub_vcf1, sub_vcf2 = split_vcf(cfg.args.sub_vcf, cfg.args.out_prefix+"_subs_hap")
    indel_vcf1, indel_vcf2 = split_vcf(cfg.args.indel_vcf, cfg.args.out_prefix+"_indels_hap")

    print(f"> indexing vcfs")
    subprocess.run(['tabix', '-f', '-p', 'vcf', sub_vcf1])
    subprocess.run(['tabix', '-f', '-p', 'vcf', sub_vcf2])
    subprocess.run(['tabix', '-f', '-p', 'vcf', indel_vcf1])
    subprocess.run(['tabix', '-f', '-p', 'vcf', indel_vcf2])

    print(f"> merging haploid SUB and INDEL VCFs")
    vcf1 = merge_vcfs(sub_vcf1, indel_vcf1, cfg.args.out_prefix+"_hap1")
    vcf2 = merge_vcfs(sub_vcf2, indel_vcf2, cfg.args.out_prefix+"_hap2")

    print(f"> indexing vcfs")
    subprocess.run([ "tabix", "-p", "vcf", vcf1 ])
    subprocess.run([ "tabix", "-p", "vcf", vcf2 ])

    print(f"> merging haplotypes")
    merge_hap_vcfs(vcf1, vcf2, cfg.args.out_prefix)



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )

    parser.add_argument("sub_vcf", type=str)
    parser.add_argument("indel_vcf", type=str)
    parser.add_argument("ref", type=str)
    parser.add_argument("out_prefix", type=str)

    parser.add_argument("--contig", type=str)
    parser.add_argument("--contig_beg", type=int)
    parser.add_argument("--contig_end", type=int)

    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
