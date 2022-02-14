import os, sys, subprocess
import argparse
import numpy as np

from vcf import *
from cig import *
from bam import *
import cfg as cfg

def main():

    print(f"> selecting vcf regions")
    get_vcf_regions()

    print("> calculating score matrices")
    subs, nps, inss, dels = get_confusion_matrices()
    cfg.args.sub_scores, cfg.args.np_scores, \
            cfg.args.ins_scores, cfg.args.del_scores = \
            calc_score_matrices(subs, nps, inss, dels)

    print(f"> splitting vcf")
    vcf1, vcf2 = split_vcf(cfg.args.vcf, cfg.args.out_prefix+"pre")

    print(f"> converting vcfs and ref to sequences")
    hap1_data = apply_vcf(vcf1, 1)
    hap2_data = apply_vcf(vcf2, 2)

    print(f"> realigning hap sequences")
    with cfg.counter.get_lock(): cfg.counter.value = 0
    with mp.Pool(10) as pool:
        data = pool.map(realign_hap, hap1_data + hap2_data)
    hap1_data = [ x for x in data if x[1] == 1 ]
    hap2_data = [ x for x in data if x[1] == 2 ]

    print('\n> generating standardized vcfs')
    vcf1 = gen_vcf(hap1_data, 1, cfg.args.out_prefix)
    vcf2 = gen_vcf(hap2_data, 2, cfg.args.out_prefix)

    print(f"> merging vcfs")
    out_fn = f"{cfg.args.out_prefix}.vcf.gz"
    merge_vcfs(vcf1, vcf2, out_fn)
    subprocess.run(['tabix', '-f', '-p', 'vcf', out_fn])



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )

    parser.add_argument("vcf", type=str)
    parser.add_argument("ref", type=str)
    parser.add_argument("out_prefix", type=str)

    parser.add_argument("--contig", type=str)
    parser.add_argument("--contig_beg", type=int)
    parser.add_argument("--contig_end", type=int)
    parser.add_argument("--contigs", type=str)

    parser.add_argument("--stats_dir", default="./stats")
    parser.add_argument("--recalc_cms", action="store_true")
    parser.add_argument("--recalc_pileups", action="store_true")

    parser.add_argument("--max_n", type=int, default=6)
    parser.add_argument("--max_l", type=int, default=100)
    parser.add_argument("--chunk_width", type=int, default=10000)

    parser.add_argument("--min_qual", type=int, default=0)
    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
