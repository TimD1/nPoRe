import os, sys
import argparse
import numpy as np

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import pysam



def main(args):

    vcf = pysam.VariantFile(args.vcf, 'rb')
    for record in vcf.fetch('chr19', 40_000_000, 58_592_616):
        print(record.pos, ':', record.alleles)



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("--vcf", type=str, default="/x/gm24385/chr19/ref/HG002_GRCh38_1_22_v4.1_draft_benchmark_std.vcf.gz")
    return parser



if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main(args)
