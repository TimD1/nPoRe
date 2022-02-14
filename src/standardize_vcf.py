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
    )

    parser.add_argument("--vcf", type=str, required=True,
            help="Input VCF to standardize.")
    parser.add_argument("--ref", type=str, required=True,
            help="Input reference FASTA corresponding to VCF.")
    parser.add_argument("--out_prefix", type=str, required=True,
            help="Output VCF prefix.")

    parser.add_argument("--contig", type=str,
            help="Allows specifying a single contig to standardize; it can be used "
            "in combination with '--contig_beg' and '--contig_end'.")
    parser.add_argument("--contig_beg", type=int,
            help='Start of standardized region.')
    parser.add_argument("--contig_end", type=int,
            help='End of standardized region.')
    parser.add_argument("--contigs", type=str,
            help="Allows specifying multiple contigs to standardize.")

    parser.add_argument("--stats_dir", default="./stats",
            help="Directory containing confusion matrices storing measured "
            "probabilities of SUBs, INDELs, and N-polymer CNVs.")

    parser.add_argument("--max_n", type=int, default=6,
            help="Maximum n-polymer length (period of repeating sequence) "
            "considered during reference standardization.")
    parser.add_argument("--max_l", type=int, default=100,
            help="Maximum length (number of times a repeated sequence occurs)"
            " considered during reference standardization.")
    parser.add_argument("--chunk_width", type=int, default=100000,
            help="BAM is considered in chunks of size '--chunk_width' at a time "
            "when calculating confusion matrices.")

    parser.add_argument("--min_qual", type=int, default=0,
            help="Only apply variants in VCF with quality above this threshold.")
    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    cfg.args.recalc_cms = False
    main()
