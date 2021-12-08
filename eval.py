import os, sys, subprocess
import pickle
import argparse
import numpy as np

from vcf import *
from cig import *
from bam import *
import cfg as cfg

def main():

    print(f"> splitting CALLS and TRUTH VCFs into haploid VCFs")
    calls_vcf1, calls_vcf2 = split_vcf(cfg.args.calls_vcf, 
            cfg.args.out_prefix+"_calls_hap")
    truth_vcf1, truth_vcf2 = split_vcf(cfg.args.truth_vcf, 
            cfg.args.out_prefix+"_truth_hap")

    print(f"> indexing VCFs")
    subprocess.run(['tabix', '-f', '-p', 'vcf', calls_vcf1])
    subprocess.run(['tabix', '-f', '-p', 'vcf', calls_vcf2])
    subprocess.run(['tabix', '-f', '-p', 'vcf', truth_vcf1])
    subprocess.run(['tabix', '-f', '-p', 'vcf', truth_vcf2])

    print(f"> reading reference")
    ref = get_fasta(cfg.args.ref, cfg.args.contig)

    print(f"> applying CALLS and TRUTH VCFs to reference")
    truth_ref1, truth_cig1 = apply_vcf(truth_vcf1, ref)
    truth_ref2, truth_cig2 = apply_vcf(truth_vcf2, ref)
    calls_seq1, calls_cig1 = apply_vcf(calls_vcf1, ref)
    calls_seq2, calls_cig2 = apply_vcf(calls_vcf2, ref)

    calls1_data = ("hap1", "chr19", 0, 0, calls_cig1, 
            truth_cig1, ref, truth_ref1, calls_seq1, 1)
    calls2_data = ("hap2", "chr19", 0, 0, calls_cig2, 
            truth_cig2, ref, truth_ref2, calls_seq2, 2)
    data = [calls1_data, calls2_data]

    with mp.Pool() as pool:

        print('> changing read basis ref->hap')
        with cfg.counter.get_lock(): cfg.counter.value = 0
        data = pool.map(to_truth_ref, data)

        print(f"\n> realigning hap sequences")
        data = pool.map(simple_realign_read, data)
        print(" ")

    for datum in data:
        read_id, ref_name, start, stop, cigar, ref, seq, hap = datum
        print(f"{read_id}:\t{cigar.count('=')}=\t{cigar.count('X')}X\t"
              f"{cigar.count('I')}I\t{cigar.count('D')}D")



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )

    parser.add_argument("calls_vcf", type=str)
    parser.add_argument("truth_vcf", type=str)
    parser.add_argument("ref", type=str)
    parser.add_argument("out_prefix", type=str)

    parser.add_argument("--contig", type=str)
    parser.add_argument("--contig_beg", type=int)
    parser.add_argument("--contig_end", type=int)
    parser.add_argument("--min_qual", type=int, default=10)

    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
