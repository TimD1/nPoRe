import os, sys, subprocess
import pickle
import argparse
import numpy as np

from vcf import *

def main(args):
    filter_overlaps(args.vcf, args.out)



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("vcf", type=str)
    parser.add_argument("out", type=str)
    return parser



if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main(args)
