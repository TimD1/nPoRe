import argparse, os, subprocess
import numpy as np
from collections import defaultdict

import pysam

import cfg
from aln import *
from bam import *


def argparser():

    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )

    # mandatory args
    parser.add_argument("bam")
    parser.add_argument("ref")
    parser.add_argument("out")

    # region of interest
    parser.add_argument("--contig", type=str, default="chr19")
    parser.add_argument("--contig_beg", type=int, default=1)
    parser.add_argument("--contig_end", type=int, default=58592616)

    # algorithm parameters
    parser.add_argument("--max_hp", type=int, default=100)
    parser.add_argument("--chunk_width", type=int, default=10000)

    # path
    parser.add_argument("--stats_dir", default="./stats")

    # boolean options
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--recalc_cms", action="store_true")
    parser.add_argument("--min_snp_qual", type=int, default=0)

    parser.add_argument("--apply_vcf")
    parser.add_argument("--std_cigar", action="store_true")
    parser.add_argument("--indels_only", action="store_true")

    return parser



def main():

    os.makedirs(cfg.args.stats_dir, exist_ok=True)
    subs, hps = get_confusion_matrices()

    if cfg.args.plot:
        print("\n> plotting confusion matrices")
        plot_confusion_matrices(subs, hps)
        print("\n> plotting distributions")
        plot_hp_len_dists(hps)

    print("\n> calculating score matrices")
    cfg.args.sub_scores, cfg.args.hp_scores = calc_score_matrices(subs, hps)

    if cfg.args.plot:
        print("> plotting score matrices")
        plot_hp_score_matrix(cfg.args.hp_scores)

    print('> extracting read data from BAM')
    read_data = get_read_data(cfg.args.bam)

    if cfg.args.apply_vcf:
        print(f"\n> splitting vcf '{cfg.args.apply_vcf}'")
        vcf1, vcf2 = split_vcf(cfg.args.apply_vcf, 
                f"{os.extsep.join(cfg.args.apply_vcf.split(os.extsep)[:-2])}")

        print(f"> indexing '{vcf1}' and\n        '{vcf2}'")
        subprocess.run(['tabix', '-p', 'vcf', vcf1])
        subprocess.run(['tabix', '-p', 'vcf', vcf2])

        print(f"> reading reference: '{cfg.args.ref}'")
        cfg.args.reference = get_fasta(cfg.args.ref, cfg.args.contig)

        print(f"> applying '{vcf1}' to reference")
        cfg.args.hap1, cfg.args.hap1_cig = apply_vcf(vcf1, cfg.args.reference)
        print(f"> applying '{vcf2}' to reference")
        cfg.args.hap2, cfg.args.hap2_cig = apply_vcf(vcf2, cfg.args.reference)

        if cfg.args.std_cigar:
            print(f"> standardizing haplotype cigars")
            _, _, _, _, _, cfg.args.hap1_cig, _, _, _, _ = \
                    standardize_cigar(("1", "chr19", 0, 0, "", cfg.args.hap1_cig, 
                        cfg.args.reference, "", cfg.args.hap1, 2))
            _, _, _, _, _, cfg.args.hap2_cig, _, _, _, _ = \
                    standardize_cigar(("2", "chr19", 0, 0, "", cfg.args.hap2_cig, 
                        cfg.args.reference, "", cfg.args.hap2, 2))
        else:
            cfg.args.hap1_cig = cfg.args.hap1_cig.replace('X','M').replace('=','M')
            cfg.args.hap2_cig = cfg.args.hap2_cig.replace('X','M').replace('=','M')

        print(f"> precomputing haplotype positions")
        cfg.args.ref_poss_hap1, cfg.args.hap1_poss = \
                get_refseq_positions(cfg.args.hap1_cig)
        cfg.args.ref_poss_hap2, cfg.args.hap2_poss = \
                get_refseq_positions(cfg.args.hap2_cig)

    with cfg.read_count.get_lock(): cfg.read_count.value = 0
    with mp.Pool() as pool:
        if cfg.args.apply_vcf: 
            print('> adding haplotype data to reads')
        read_data = list(filter(None, pool.map(add_haplotype_data, read_data)))
        if cfg.args.apply_vcf:
            read_data = pool.map(to_haplotype_ref, read_data)

    print('\n> computing read realignments')
    with cfg.read_count.get_lock(): cfg.read_count.value = 0
    with mp.Pool() as pool:
        print(f"\r    0 reads realigned.", end='', flush=True)
        read_data = pool.map(realign_read, read_data)

    with cfg.read_count.get_lock(): cfg.read_count.value = 0
    with mp.Pool() as pool:
        if cfg.args.std_cigar:
            print('\n> converting to standard INDEL format')
            read_data = pool.map(standardize_cigar, read_data)
        if cfg.args.apply_vcf:
            read_data = pool.map(from_haplotype_ref, read_data)

    print(f"\n> saving results to '{cfg.args.out}'")
    write_results(read_data, cfg.args.out)
    print("\n")



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()

    if cfg.args.indels_only and not cfg.args.std_cigar:
        print("ERROR: cannot set 'indels_only' without 'std_cigar'.")
        exit(1)

    main()
