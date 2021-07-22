import os, sys
import argparse
import numpy as np

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from aln import align, calc_score_matrices, dump
from cig import expand_cigar, standardize_cigar
import cfg as cfg
from bam import get_confusion_matrices



def main():

    # create test cases
    ref_rd_cigs = []
    # ref_rd_cigs.append(("ACCAGGCAT", "ACCAGGCAT", "9="))
    # ref_rd_cigs.append(("ACCAGGCAT", "ACAGGCA", "2=1D5=1D"))
    # ref_rd_cigs.append(("ACCAGGCAT", "ACCCAGGAT", "1=1I5=1D2="))
    # ref_rd_cigs.append(("AAAACCAGGCA", "AAACCAGGCA", "1D10="))
    # ref_rd_cigs.append(("TAAACCAGGCA", "AAACCAGGCA", "1D10="))
    # ref_rd_cigs.append(("AAAACCAGGCA", "AAAAACCAGGCA", "1I11="))
    # ref_rd_cigs.append(("AAAACCAGGCA", "TAAAACCAGGCA", "1I11="))
    ref_rd_cigs.append(("CCAAAAAATTTTTCC", "CCAAAAATTTTTTCC", "7=1X7="))
    ref_rd_cigs.append(("CACACACATATATATAGG", "CACACACATATATAGG", "14=2D2="))
    ref_rd_cigs.append(("CACACACATATATATAGG", "CACACACATATATATATAGG", "16=2I2="))
    ref_rd_cigs.append(("AACAACAACAACAAAAA", "AACAACAACAAAAA", "10=3D4="))

    # load existing score mats
    subs, hps = get_confusion_matrices()
    sub_scores, hp_scores = calc_score_matrices(subs, hps)

    # test alignment algorithm
    print('\n> aligning...')
    for ref, seq, cigar in ref_rd_cigs:

        int_ref = np.zeros(len(ref), dtype=np.uint8)
        for i in range(len(ref)): 
            int_ref[i] = cfg.base_dict[ref[i]]
        int_seq = np.zeros(len(seq), dtype=np.uint8)
        for i in range(len(seq)): 
            int_seq[i] = cfg.base_dict[seq[i]]

        # new_cigar = align(int_ref, int_seq, cigar, sub_scores, hp_scores, r=4)
        # print(f"Cigar: {expand_cigar(cigar)}")
        # dump(ref, seq, new_cigar)

        read_data = ("read_id", "ref", 0, cigar, ref, seq)
        read_id, ref_name, start, new_cigar, ref, seq = standardize_cigar(read_data)

        print(f"\nCigar: {expand_cigar(cigar)}")
        dump(ref, seq, expand_cigar(new_cigar))



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("--recalc_cms", action="store_true")
    parser.add_argument("--max_hp", type=int, default=100)
    parser.add_argument("--stats_dir", type=str, default="../stats")
    parser.add_argument("--indels_only", type=bool, default=True)
    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
