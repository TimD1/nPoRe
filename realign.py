import argparse
import numpy as np
import multiprocessing as mp
from collections import defaultdict
from numba import njit

import pysam




class Cigar():
    '''
    Enum for pysam's cigartuples encoding.
    '''
    M = 0
    I = 1
    D = 2
    N = 3
    S = 4
    H = 5
    P = 6
    E = 7
    X = 8
    B = 9



def get_positions(vcf_filename: str, min_qual: int, window: int):
    '''
    Parse VCF file and extract candidate variant positions.
     - VCF filename must be zipped (supply <filename>.vcf.gz)
     - VCF index must be present (<filename>.vcf.gz.tbi)
    '''

    # get candidate variants from VCF
    vcf = pysam.VariantFile(vcf_filename, 'r')
    all_pos = [record.start for record in vcf.fetch() if record.qual > min_qual]

    # filter nearby variant candidates (within window)
    last_pos = 0
    pos = []
    for p in all_pos:
        if p >= last_pos + 2*window:
            pos.append(p)
            last_pos = p
    return pos



def realign_bam(positions):
    
    # with mp.Pool() as pool:
    #     read_alignments = pool.map(realign_pos, positions)

    read_alignments = realign_pos(positions[0])
    



def realign_pos(pos):

    alignments = defaultdict(list)
    bam = pysam.AlignmentFile(args.bam, 'rb')

    ref_start = pos-args.window
    ref_end = pos+args.window+args.model.kmer_len-1
    print(f'pos {pos}: {ref_start}-{ref_end}')

    for read in bam.fetch(args.contig, pos, pos+1):

        # skip read if it doesn't fully cover the window
        ref = read.get_reference_sequence().upper() \
                [ref_start-read.reference_start : ref_end-read.reference_start]
        if len(ref) < args.window*2:
            continue

        # find read substring overlapping region
        cigar_types = [ c[0] for c in read.cigartuples ]
        cigar_counts = [ c[1] for c in read.cigartuples ]
        read_idx, ref_idx = 0, read.reference_start
        read_start, read_end, first = None, None, True

        while ref_idx < ref_end:

            # first read position overlapping region
            if first and ref_idx >= ref_start:
                first = False
                read_start = read_idx

            read_move, ref_move = None, None
            cigar = cigar_types[0]

            # determine whether to move on read/ref
            if cigar == Cigar.S:    # soft-clipped
                read_move = True
                ref_move = False
            elif cigar == Cigar.H:    # hard-clipped
                read_move = False
                ref_move = False
            elif cigar == Cigar.X:    # substitution
                read_move = True
                ref_move = True
            elif cigar == Cigar.I:    # insertion
                read_move = True
                ref_move = False
            elif cigar == Cigar.D:    # deletion
                read_move = False
                ref_move = True
            elif cigar == Cigar.E:    # match
                read_move = True
                ref_move = True
            elif cigar == Cigar.M:    # match/sub
                read_move = True
                ref_move = True
            else:
                print("ERROR: unexpected CIGAR type for {}.".format(read.query_name))

            # shift reference index by one base or deleted section
            if ref_move:
                if read_move:
                    ref_idx += 1
                else:
                    ref_idx += cigar_counts[0]

            # shift read index
            if read_move:
                cigar_counts[0] -= 1
                read_idx += 1
            else:
                cigar_counts[0] = 0

            # move to next CIGAR section of interest
            if cigar_counts[0] == 0:
                del cigar_counts[0]
                del cigar_types[0]

        # extract read section
        read_end = read_idx
        seq = read.query_sequence[read_start:read_end]
        print(f'\nref:\t{ref[:-args.model.kmer_len+1]}')
        print(f'seq:\t{seq[:-args.model.kmer_len+1]}')

        # convert both sequences to expected current levels
        raw_ref = args.model.convert(ref)
        raw_seq = args.model.convert(seq)
        # print(f'raw_ref: {raw_ref}\traw_seq: {raw_seq}')

        # perform DTW alignment in signal space
        raw_cigar = dtw_align(raw_seq, raw_ref, seq, ref)
        print(f'cigar:\t{raw_cigar}')






        # alignments[read.alignment.query_name].append((ref_start, ref_end, cigar))

    return alignments


@njit()
def dtw_align(raw_seq, raw_ref, seq, ref):

    # define path mat encoding
    EQ = 0
    SUB = 1
    INS = 2
    DEL = 3

    # initialize path matrix
    path_mat = np.zeros((len(raw_seq), len(raw_ref))) # ignore k-mer padding
    for i in range(len(raw_seq)):
        path_mat[i, 0] = INS
    for j in range(len(raw_ref)):
        path_mat[0, j] = DEL
    path_mat[0, 0] = EQ if ref[0] == seq[0] else SUB

    # initialize cost matrix
    cost_mat = np.zeros((len(raw_seq), len(raw_ref)))
    cost_mat[0, 0] = (raw_seq[0] - raw_ref[0]) ** 2
    for i in range(1, len(raw_seq)):
        cost_mat[i, 0] = cost_mat[i-1, 0] + (raw_seq[i] - raw_ref[0]) ** 2

    # compute entire cost matrix
    for i in range(1, len(raw_seq)):
        for j in range(1, len(raw_ref)):
            if seq[i] == ref[j]:
                diag = cost_mat[i-1, j-1]
            else:
                diag = cost_mat[i-1, j-1] + 400
            left = cost_mat[i, j-1] + 400
            top = cost_mat[i-1, j] + 400
            cost_mat[i, j] = (raw_seq[i] - raw_ref[j])**2 + min(diag, top, left)

            # fill in path matrix simultaneously
            if diag <= left:
                if diag <= top:
                    if seq[i] == ref[j]:
                        path_mat[i,j] = EQ
                    else:
                        path_mat[i,j] = SUB
                else:
                    path_mat[i,j] = INS
            else:
                if left <= top:
                    path_mat[i,j] = DEL
                else:
                    path_mat[i,j] = INS

    # return alignment path
    cigar = ""
    types = "=XID"
    i, j = len(raw_seq)-1, len(raw_ref)-1
    while i >= 0 and j >= 0:
        if path_mat[i,j] == EQ:
            cigar += types[EQ]
            i -= 1
            j -= 1
        elif path_mat[i, j] == SUB:
            cigar += types[SUB]
            i -= 1
            j -= 1
        elif path_mat[i, j] == INS:
            cigar += types[INS]
            i -= 1
        elif path_mat[i, j] == DEL:
            cigar += types[DEL]
            j -= 1
    return cigar[::-1]



class KmerModel():
    '''
    Class for nanopore k-mer based current model.
    '''

    def __init__(self, model_filename):
        self.model = {}

        try:
            with open(model_filename, 'r') as model_file:
                for line in model_file:
                    kmer, value = line.strip().split()
                    self.model[kmer] = float(value)
                self.kmer_len = len(kmer)

        except FileNotFoundError:
            print(f"Model file '{model_file}' not found.")
            exit(1)


    def convert(self, seq):
        signal = np.zeros(len(seq)-self.kmer_len+1)
        for kmer_start in range(len(signal)):
            signal[kmer_start] = self.model[ \
                    seq[kmer_start : kmer_start+self.kmer_len] 
                ]
        return signal



def argparser():

    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )

    parser.add_argument("bam")
    parser.add_argument("ref")
    parser.add_argument("vcf")
    parser.add_argument("out")

    parser.add_argument("--model_file", default="dna_6mer_model.txt")
    parser.add_argument("--contig", default="chr19")
    parser.add_argument("--min_qual", type=int, default=0)
    parser.add_argument("--window", type=int, default=25)

    return parser



def main():

    print("> parsing pore model")
    args.model = KmerModel(args.model_file)

    print("> getting DeepVariant positions")
    positions = get_positions(args.vcf, args.min_qual, args.window)

    print("> computing realignments")
    alignments = realign_bam(positions)

    print("> saving results")
    # write_results(alignments, bam, args.out)

args = None
if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main()



