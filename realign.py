import argparse, os, subprocess
import numpy as np
from collections import defaultdict
from time import perf_counter

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
    parser.add_argument("out_prefix")

    # region of interest
    parser.add_argument("--contig", type=str, default="chr19")
    parser.add_argument("--contig_beg", type=int, default=1)
    parser.add_argument("--contig_end", type=int, default=58592616)
    parser.add_argument("--max_reads", type=int, default=0)
    parser.add_argument("--consensus_itrs", type=int, default=3)

    # algorithm parameters
    parser.add_argument("--max_np", type=int, default=10)
    parser.add_argument("--max_np_len", type=int, default=100)
    parser.add_argument("--chunk_width", type=int, default=10000)

    # path
    parser.add_argument("--stats_dir", default="./stats")

    # boolean options
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--recalc_cms", action="store_true")
    parser.add_argument("--recalc_pileups", action="store_true")
    parser.add_argument("--indel_cigar", action="store_true")


    return parser



def main():

    print("> loading confusion matrices")
    os.makedirs(cfg.args.stats_dir, exist_ok=True)
    subs, nps, inss, dels = get_confusion_matrices()

    print("> calculating score matrices")
    cfg.args.sub_scores, cfg.args.np_scores, cfg.args.ins_scores, cfg.args.del_scores = \
            calc_score_matrices(subs, nps, inss, dels)

    if cfg.args.plot:
        print("> plotting confusion matrices")
        plot_confusion_matrices(subs, nps, inss, dels)

        print("> plotting score matrices")
        plot_np_score_matrices(cfg.args.np_scores)
        exit(0)

    print('> extracting read data from BAM')
    start = perf_counter()
    read_data = get_read_data(cfg.args.bam)
    print(f'\n    runtime: {perf_counter()-start:.2f}s')

    with cfg.counter.get_lock(): cfg.counter.value = 0
    with mp.Pool() as pool:

        print('> computing individual read realignments')
        with cfg.counter.get_lock(): cfg.counter.value = 0
        start = perf_counter()
        print(f"\r    0 reads realigned.", end='', flush=True)
        read_data = pool.map(realign_read, read_data)
        print(f'\n    runtime: {perf_counter()-start:.2f}s')

        with cfg.counter.get_lock(): cfg.counter.value = 0
        print('> converting to standard CIGAR format')
        start = perf_counter()
        read_data = pool.map(standardize_cigar, read_data)
        print(f'\n    runtime: {perf_counter()-start:.2f}s')

    
    cfg.args.itr = 0
    print(f"> saving intermediary results to '{cfg.args.out_prefix}{cfg.args.itr}.bam'")
    write_results(read_data, f'{cfg.args.out_prefix}{cfg.args.itr}.bam')

    while cfg.args.itr < cfg.args.consensus_itrs:
        get_pileup_info()
        cfg.args.itr += 1
        with mp.Pool() as pool:
            print(f'> computing consensus read realignments, iteration {cfg.args.itr}')
            start = perf_counter()
            with cfg.counter.get_lock(): cfg.counter.value = 0
            print(f"\r    0 reads realigned.", end='', flush=True)
            read_data = pool.map(realign_read2, read_data)
            print(f'\n    runtime: {perf_counter()-start:.2f}s')

        print(f"> saving results to '{cfg.args.out_prefix}{cfg.args.itr}.bam'")
        write_results(read_data, f'{cfg.args.out_prefix}{cfg.args.itr}.bam')
        print("\n")



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    try:
        main()
    except KeyboardInterrupt:
        print("\nERROR: Program terminated.")
        exit(1)
