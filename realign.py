import argparse, os, subprocess
import numpy as np
from collections import defaultdict
from time import perf_counter
import multiprocessing as mp
from ctypes import c_wchar_p

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
    parser.add_argument("--contig", type=str)
    parser.add_argument("--contigs", type=str)
    parser.add_argument("--contig_beg", type=int)
    parser.add_argument("--contig_end", type=int)
    parser.add_argument("--max_reads", type=int, default=0)

    # algorithm parameters
    parser.add_argument("--max_np", type=int, default=6)
    parser.add_argument("--max_np_len", type=int, default=100)
    parser.add_argument("--chunk_width", type=int, default=1000000)

    # path
    parser.add_argument("--stats_dir", default="./stats")

    # boolean options
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--recalc_cms", action="store_true")
    parser.add_argument("--recalc_exit", action="store_true")
    parser.add_argument("--indel_cigar", action="store_true")


    return parser



def main():

    print("> selecting BAM regions")
    get_bam_regions()

    # loading confusion matrices
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

    print('> creating output SAM')
    create_header(f'{cfg.args.out_prefix}.sam')

    print('> extracting read data from BAM')
    read_data = get_read_data(cfg.args.bam)

    start = perf_counter()
    with mp.Pool() as pool:
        print('> computing individual read realignments')
        with cfg.counter.get_lock(): cfg.counter.value = 0
        print(f"\r    0 reads processed.", end='', flush=True)
        read_data = pool.imap_unordered(realign_read, read_data, chunksize=100)
        pool.close()
        pool.join()
    print(f'\n    runtime: {perf_counter()-start:.2f}s')



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()

    try:
        main()
    except KeyboardInterrupt:
        print("\nERROR: Program terminated.")
        exit(1)
