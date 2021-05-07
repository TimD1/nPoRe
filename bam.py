import multiprocessing as mp
from collections import defaultdict
import numpy as np

import pysam

import cfg

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



def realign_bam(positions):
    '''
    Wrapper for multi-threaded read re-alignment at each BAM position.
    '''
    
    # with mp.Pool() as pool:
    #     read_alignments = pool.map(realign_pos, positions)

    read_alignments = realign_pos(positions[0])
    


def realign_pos(pos):
    '''
    Re-align all reads covering a single reference column (within window).
    '''

    alignments = defaultdict(list)
    bam = pysam.AlignmentFile(cfg.args.bam, 'rb')

    ref_start = pos-cfg.args.window
    ref_end = pos+cfg.args.window
    print(f'pos {pos}: {ref_start}-{ref_end}')

    for read in bam.fetch(cfg.args.contig, pos, pos+1):

        # skip read if it doesn't fully cover the window
        ref = read.get_reference_sequence().upper() \
                [ref_start-read.reference_start : ref_end-read.reference_start]
        if len(ref) < ref_end - ref_start: continue

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
        print(f'\nref:\t{ref}')
        print(f'seq:\t{seq}')


        # alignments[read.alignment.query_name].append((ref_start, ref_end, cigar))

    return alignments

def write_bam(fname, alignments, header, bam=True):
    '''
    Write a `.bam` file for a set of alignments.
    '''
    with pysam.AlignmentFile(fname, 'wb', header=header) as fh:
        for ref_id, subreads in enumerate(alignments):
            for aln in sorted(subreads, key=lambda x: x.rstart):
                a = pysam.AlignedSegment()
                a.reference_id = ref_id
                a.query_name = aln.qname
                a.query_sequence = aln.seq
                a.reference_start = aln.rstart
                a.cigarstring = aln.cigar
                a.flag = aln.flag
                a.mapping_quality = 60
                fh.write(a)
    pysam.index(fname) 
