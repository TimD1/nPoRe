import argparse, os
import numpy as np
from collections import defaultdict

import pysam

import cfg
from vcf import get_positions
from bam import realign_bam, BamStats


def argparser():

    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )

    parser.add_argument("bam")
    parser.add_argument("ref")
    parser.add_argument("vcf")
    parser.add_argument("out")

    parser.add_argument("--contig", default="chr19")
    parser.add_argument("--contig_beg", type=int, default=100000)
    parser.add_argument("--contig_end", type=int, default=200000)

    parser.add_argument("--min_qual", type=int, default=0)
    parser.add_argument("--max_hp", type=int, default=10)
    parser.add_argument("--window", type=int, default=25)

    parser.add_argument("--stats_dir", default="./stats")

    parser.add_argument("--force", action="store_true")

    return parser



def main():

    print("> calculating BAM statistics")
    os.makedirs(cfg.args.stats_dir, exist_ok=True)
    stats = BamStats(cfg.args.bam)
    exit(0)

    print("> getting DeepVariant positions")
    positions = get_positions(cfg.args.vcf, cfg.args.min_qual, cfg.args.window)

    print("> computing realignments")
    alignments = realign_bam(positions)

    print("> saving results")
    # write_results(alignments, bam, args.out)



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
