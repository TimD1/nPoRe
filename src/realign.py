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
    )

    # mandatory args
    parser.add_argument("--bam", required=True,
            help="Input BAM to be realigned.")
    parser.add_argument("--ref", required=True,
            help="Input reference FASTA.")
    parser.add_argument("--out_prefix", required=True,
            help="Output SAM file prefix.")

    # region of interest
    parser.add_argument("--contig", type=str,
            help="Allows specifying a single contig to realign, can be used "
            "in combination with '--contig_beg' and '--contig_end'. Reads outside "
            "this region are ignored and not included in the output.")
    parser.add_argument("--contig_beg", type=int,
            help='Start of realigned region.')
    parser.add_argument("--contig_end", type=int,
            help='End of realigned region.')
    parser.add_argument("--contigs", type=str,
            help="Allows specifying multiple contigs to realign.")
    parser.add_argument("--max_reads", type=int, default=0,
            help="Allows limiting the number of realigned reads. If zero, all "
            "reads in the selected region are realigned.")
    parser.add_argument("--bed", type=str,
            help="Allows specifying arbitrary realigned regions using a BED file.")

    # algorithm parameters
    parser.add_argument("--max_n", type=int, default=6,
            help="Maximum n-polymer length (period of repeating sequence) "
            "considered during read realignment.")
    parser.add_argument("--max_l", type=int, default=100,
            help="Maximum length (number of times a repeated sequence occurs)"
            " considered during read realignment.")
    parser.add_argument("--chunk_width", type=int, default=100000,
            help="BAM is considered in chunks of size '--chunk_width' "
            "when calculating confusion matrices.")

    # path
    parser.add_argument("--stats_dir", default="./stats",
            help="Directory containing confusion matrices storing measured "
            "probabilities of SUBs, INDELs, and N-polymer CNVs.")

    # boolean options
    parser.add_argument("--plot", action="store_true",
            help="Plot confusion and score matrices, exit without realigning.")
    parser.add_argument("--recalc_cms", action="store_true",
            help="Recalculate confusion matrices using provided BAM.")
    parser.add_argument("--recalc_exit", action="store_true",
            help="Exit after recalculating confusion matrices, without realigning.")


    return parser



def main():

    print("> selecting BAM regions")
    get_bam_regions()

    print("> reading reference")
    if cfg.args.recalc_cms:
        cfg.args.refs = {}
        for ctg, _, _ in cfg.args.regions:
            if ctg not in cfg.args.refs.keys():
                cfg.args.refs[ctg] = get_fasta(cfg.args.ref, ctg)

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
        pool.imap_unordered(realign_read, read_data, chunksize=100)
        pool.close()
        pool.join()
    print(f'    runtime: {perf_counter()-start:.2f}s')



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()

    # if all stats not available, recalculate cms
    if not (os.path.isfile(f'{cfg.args.stats_dir}/subs_cm.npy') and \
            os.path.isfile(f'{cfg.args.stats_dir}/nps_cm.npy') and \
            os.path.isfile(f'{cfg.args.stats_dir}/inss_cm.npy') and \
            os.path.isfile(f'{cfg.args.stats_dir}/dels_cm.npy')):
        cfg.args.recalc_cms = True

    try:
        main()
    except KeyboardInterrupt:
        print("\nERROR: Program terminated.")
        exit(1)
