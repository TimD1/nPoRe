import os, sys
import argparse
import numpy as np

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from aln import get_ns, get_np_lens
import cfg

# create test cases
seqs = []
seqs.append("ATATATTTTTTTAAA")
seqs.append("ATATATATATATATATATATTTAA")
seqs.append("ACGATCTCTAGGCAGTTAGCCGAGCAG")
seqs.append("ACCGGCGCAGCAGCAGCAG")
seqs.append("TATATATATGCGCGCGGGGATATA")

def print_results(seq, ns, np_lens):
    print("\nseq:      ", end="")
    print("  ".join([base for base in seq]))
    print("ns:      ", end="")
    print(" ".join([f"{n:2}" for n in ns]))
    print("np_lens: ", end="")
    print(" ".join([f"{l:2}" for l in np_lens]))


def main():

    for seq in seqs:
        int_seq = np.zeros(len(seq), dtype=np.uint8)
        for i in range(len(seq)): 
            int_seq[i] = cfg.base_dict[seq[i]]
        ns = get_ns(int_seq)
        np_lens = get_np_lens(int_seq)
        print_results(seq, ns, np_lens)

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
