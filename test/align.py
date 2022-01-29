import os, sys
import argparse
import numpy as np

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
    # ref_rd_cigs.append(("ACCAGGCAT", "ACCCAGGAT", "1=1I5=1D2="))
    # ref_rd_cigs.append(("AAAACCAGGCA", "AAACCAGGCA", "1D10="))
    # ref_rd_cigs.append(("TAAACCAGGCA", "AAACCAGGCA", "1D10="))
    # ref_rd_cigs.append(("AAAACCAGGCA", "AAAAACCAGGCA", "1I11="))
    # ref_rd_cigs.append(("AAAACCAGGCA", "TAAAACCAGGCA", "1I11="))
    # ref_rd_cigs.append(("CCAAAAAATTTTTCC", "CCAAAAATTTTTTCC", "7=1X7="))
    # ref_rd_cigs.append(("CACACACATATATATAGG", "CACACACATATATAGG", "14=2D2="))
    # ref_rd_cigs.append(("CACACACATATATATAGG", "CACACACATATATATATAGG", "16=2I2="))
    # ref_rd_cigs.append(("AACAACAACAACAAAAA", "AACAACAACAAAAA", "10=3D4="))
    # ref_rd_cigs.append(("GCACAGCAGTC", "GCACAGTC", "1=2D2=1D5="))
    # ref_rd_cigs.append(("AAAAAAAA", "AAAAAA", "1=1D3=1D2="))
    ref_rd_cigs.append(("CAAAGAAAGAAAG", "CAAAGAAAGAAG", "9=1D3="))
    ref_rd_cigs.append(("CAAAGAAAGAAAG", "CAAAGAAAAGAAAG", "5=1I8="))
    ref_rd_cigs.append(("CAAAGAAAGAAAG", "CAAAGAAAAG", "5=4D1I4="))
    ref_rd_cigs.append(("CAAAGAAAGAAAG", "CAAGAAAG", "1=5D7="))
    ref_rd_cigs.append(("CGAAAGAAAGAAAG", "CGAAGAAAG", "2=5D7="))
    ref_rd_cigs.append(("CGAAAGAAAGAAAC", "CGAAGAAAC", "2=5D7="))
    # ref_rd_cigs.append(("ATATATATTTTTTAAAGCGCGC", "ATATATATTTTTTAAAGCGCGC", "22="))

    # load existing score mats
    subs, nps, inss, dels = get_confusion_matrices()
    sub_scores, np_scores, ins_scores, del_scores = calc_score_matrices(subs, nps, inss, dels)

    # test alignment algorithm
    print('\n> aligning...')
    for ref, seq, cigar in ref_rd_cigs:
        print("\n\n\n")
        cigar = expand_cigar(cigar)

        int_ref = np.zeros(len(ref), dtype=np.uint8)
        for i in range(len(ref)): 
            int_ref[i] = cfg.base_dict[ref[i]]
        int_seq = np.zeros(len(seq), dtype=np.uint8)
        for i in range(len(seq)): 
            int_seq[i] = cfg.base_dict[seq[i]]

        print("-"*80)
        new_cigar = align(int_ref, int_seq, cigar, sub_scores, np_scores, verbose=True, r=10)
        print(f"Cigar: {cigar}")
        dump(ref, seq, new_cigar)



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("--recalc_cms", action="store_true")
    parser.add_argument("--max_n", type=int, default=6)
    parser.add_argument("--max_l", type=int, default=100)
    parser.add_argument("--stats_dir", type=str, default="../old_4stats")
    parser.add_argument("--indels_only", type=bool, default=True)
    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
