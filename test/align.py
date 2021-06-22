import os, sys
import argparse

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from aln import align, calc_score_matrices, dump
from cig import expand_cigar
import cfg as cfg
from bam import get_confusion_matrices



def main():

    # create test cases
    ref_rd_cigs = []
    # ref_rd_cigs.append(("ACCAGGCAT", "ACCAGGCAT", "9="))
    # ref_rd_cigs.append(("ACCAGGCAT", "ACAGGCA", "2=1D5=1D"))
    ref_rd_cigs.append(("ACCAGGCAT", "ACCCAGGAT", "1=1I5=1D2="))

    # load existing score mats
    subs, hps = get_confusion_matrices()
    sub_scores, hp_scores = calc_score_matrices(subs, hps)

    # test alignment algorithm
    print('\n> aligning...')
    for ref, read, cigar in ref_rd_cigs:
        new_cigar = align(ref, read, ref, cigar, sub_scores, hp_scores, verbose=True)
        print(f"Cigar: {expand_cigar(cigar)}")
        dump(ref, read, new_cigar)



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("--recalc_cms", action="store_true")
    parser.add_argument("--max_hp", type=int, default=100)
    parser.add_argument("--stats_dir", type=str, default="../stats")
    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
