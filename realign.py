import argparse, os
import numpy as np
from collections import defaultdict

import pysam

import cfg
from vcf import get_positions
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
    parser.add_argument("vcf")
    parser.add_argument("out")

    # region of interest
    parser.add_argument("--contig", type=str, default="")
    parser.add_argument("--contig_beg", type=int, default=0)
    parser.add_argument("--contig_end", type=int, default=55000000)

    # algorithm parameters
    parser.add_argument("--min_qual", type=int, default=0)
    parser.add_argument("--max_hp", type=int, default=100)
    parser.add_argument("--window", type=int, default=25)
    parser.add_argument("--chunk_width", type=int, default=10000)

    # path
    parser.add_argument("--stats_dir", default="./stats")

    # boolean options
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--recalc_cms", action="store_true")
    parser.add_argument("--splice_subs", action="store_true")
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

    print("\n> computing BAM realignments")
    alignments = realign_bam()

    print(f"\n\n> saving results to '{cfg.args.out}'")
    write_results(alignments, cfg.args.out)
    print("\n")



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
