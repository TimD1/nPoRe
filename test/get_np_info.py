import os, sys
import argparse
import numpy as np

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from aln import get_np_info
import cfg

# create test cases
seqs = []
seqs.append("ATATATTTTTTTAAA")
seqs.append("ATATATATATATATATATATTTAA")
seqs.append("ACGATCTCTAGGCAGTTAGCCGAGCAG")
seqs.append("ACCGGCGCAGCAGCAGCAG")
seqs.append("TATATATATGCGCGCGGGGATATA")

def print_results(seq, np_info):
    print("\nseq:    ", end="")
    print("  ".join([base for base in seq]))
    print("N:     ", end="")
    print(" ".join([f"{i:2}" for i in np_info[0]]))
    print("MAX:   ", end="")
    print(" ".join([f"{i:2}" for i in np_info[1]]))
    print("RPT:   ", end="")
    print(" ".join([f"{i:2}" for i in np_info[2]]))
    print("IDX:   ", end="")
    print(" ".join([f"{i:2}" for i in np_info[3]]))


def main():

    for seq in seqs:
        int_seq = np.zeros(len(seq), dtype=np.uint8)
        for i in range(len(seq)): 
            int_seq[i] = cfg.base_dict[seq[i]]
        np_info = get_np_info(int_seq)
        print_results(seq, np_info)

def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("--max_n", type=int, default=10)
    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
