import argparse

from vcf import *

def main(args):
    filter_overlaps(args.vcf, args.out)



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--vcf", type=str, required=True,
            help="Input VCF from which to filter overlaps.")
    parser.add_argument("--out", type=str, required=True,
            help="Output VCF filename.")
    return parser



if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main(args)
