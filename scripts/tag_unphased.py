import argparse, sys, random
import pysam
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp

def main(args):

    in_bam = pysam.AlignmentFile(args.in_bam, 'rb')
    out_bam = pysam.AlignmentFile(args.out_bam, 'wb', header=in_bam.header)

    for read in in_bam.fetch():
        if not read.has_tag('HP'):
            read.set_tag('HP', 0, 'i')
        out_bam.write(read)



if __name__ == "__main__":

    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_bam')
    parser.add_argument('--out_bam')
    args = parser.parse_args()

    # call main
    main(args)
