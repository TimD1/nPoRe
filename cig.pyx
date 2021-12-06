import multiprocessing as mp
from collections import defaultdict
import numpy as np
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt
import os, re, itertools, sys

import pysam, cython

import cfg


def collapse_cigar(extended_cigar, return_groups=False):
    ''' 
    Converts extended CIGAR ops list to normal CIGAR string. 
    'DMMMII' -> '1D3M2I'
    '''
    count = 1
    last = None
    groups = []
    for op in extended_cigar:
        if last and op == last:
            count += 1
        elif last:
            groups.append((count, last))
            count = 1
        last = op

    if last:
        groups.append((count, last))

    if return_groups:
        return groups

    out = ''
    for num, op in groups:
        out += '%s%s' % (num, op)
    return out



def expand_cigar(cigar):
    ''' 
    Converts CIGAR string to list of ops. 
    '1D3M2I' -> 'DMMMII'
    '''

    cigar_str = ''
    count = 0

    for char in cigar:
        if char in '0123456789':
            count = count * 10 + ord(char) - ord('0')
        else:
            cigar_str += count * char
            count = 0
    return cigar_str



class Cigar():
    ''' Enum for pysam's cigartuples encoding.  '''
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



def get_cigars(bam):
    ''' Returns dicts of all read starts and CIGARs in BAM file. '''

    cigars = dict()
    starts = dict()
    bam = pysam.AlignmentFile(cfg.args.bam, 'rb')

    for read in bam.fetch():
        cigars[read.query_name] = read.cigartuples
        starts[read.query_name] = read.reference_start

    return cigars, starts



def extend_pysam_cigar(ops, counts):
    ''' 
    Converts split PySam cigartuples lists to extended CIGAR. 
    ['D', 'M', 'I'], [1, 3, 2] -> 'DMMMII'
    '''
    return ''.join([int(count)*cfg.cigars[op] for (count, op) in zip(counts, ops)])



@cython.boundscheck(False)
@cython.wraparound(False)
cpdef push_indels_left(char[::1] cigar, char[::1] seq, char[::1] nshifts_buf, 
        char[::1] shiftlen_buf, char push_op):
    ''' Push CIGAR indels leftwards. '''
    cdef char M = 0
    cdef char E = 7
    cdef char X = 8

    cdef int seq_ptr = 0
    cdef int cig_ptr = 0
    cdef int cig_len = len(cigar)
    cdef int i, indel_len, nshifts
    cdef char op

    while cig_ptr < cig_len:

        # get indel length
        op = cigar[cig_ptr]
        if op == push_op:
            indel_len = 1
            while cig_ptr + indel_len < cig_len and \
                    cigar[cig_ptr + indel_len] == push_op:
                indel_len += 1
        else:
            cig_ptr += 1
            if op == M or op == X or op == E:
                seq_ptr += 1
            continue

        # push indel as far left as possible (keeping seq same)
        nshifts = 0
        while cig_ptr-nshifts > 0 and seq_ptr-nshifts > 0 and \
                seq[seq_ptr-nshifts-1] == seq[seq_ptr-nshifts-1 + indel_len] and \
                (cigar[cig_ptr-nshifts-1] == E or cigar[cig_ptr-nshifts-1] == M) :
            nshifts += 1

        if nshifts:
            # fill buffers (can't fully update in-place)
            for i in range(nshifts):
                nshifts_buf[i] = cigar[cig_ptr-nshifts+i]
            for i in range(indel_len):
                shiftlen_buf[i] = cigar[cig_ptr+i]

            # update cigar
            for i in range(indel_len):
                cigar[cig_ptr-nshifts+i] = shiftlen_buf[i]
            for i in range(nshifts):
                cigar[cig_ptr-nshifts+indel_len+i] = nshifts_buf[i]

        # print(f'pushed {indel_len}{cfg.cigars[push_op]} left {nshifts} @ pos {seq_ptr}')

        # update pointers
        cig_ptr += indel_len
        if op == M or op == X or op == E:
            seq_ptr += 1
        elif op == push_op:
            seq_ptr += indel_len

    return cigar


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef push_indels_right(char[::1] cigar, char[::1] seq, 
        char[::1] nshifts_buf, char[::1] shiftlen_buf, char push_op):
    ''' Push CIGAR indels rightwards. '''
    cdef char M = 0
    cdef char E = 7
    cdef char X = 8

    cdef int cig_len = len(cigar)
    cdef int seq_len = len(seq)
    cdef int seq_ptr = seq_len-1
    cdef int cig_ptr = cig_len-1
    cdef int i, indel_len, nshifts
    cdef char op

    while cig_ptr >= 0:

        # get indel length
        op = cigar[cig_ptr]
        if op == push_op:
            indel_len = 1
            while cig_ptr - indel_len >= 0 and \
                    cigar[cig_ptr - indel_len] == push_op:
                indel_len += 1
        else:
            cig_ptr -= 1
            if op == M or op == X or op == E:
                seq_ptr -= 1
            continue

        # push indel as far right as possible (keeping seq same)
        nshifts = 0
        while cig_ptr+nshifts+1 < cig_len and seq_ptr+nshifts+1 < seq_len and \
                seq[seq_ptr+nshifts+1] == seq[seq_ptr+nshifts+1 - indel_len] and \
                (cigar[cig_ptr+nshifts+1] == E or cigar[cig_ptr+nshifts+1] == M) :
            nshifts += 1

        if nshifts:
            # fill buffers (can't fully update in-place)
            for i in range(nshifts):
                nshifts_buf[i] = cigar[cig_ptr+i+1]
            for i in range(indel_len):
                shiftlen_buf[i] = cigar[cig_ptr-indel_len+i+1]

            # update cigar
            for i in range(indel_len):
                cigar[cig_ptr-indel_len+nshifts+i+1] = shiftlen_buf[i]
            for i in range(nshifts):
                cigar[cig_ptr-indel_len+i+1] = nshifts_buf[i]

        # update pointers
        cig_ptr -= indel_len
        if op == M or op == X or op == E:
            seq_ptr -= 1
        elif op == push_op:
            seq_ptr -= indel_len

    return cigar


@cython.boundscheck(False)
@cython.wraparound(False)
cdef push_inss_thru_dels(char[::1] cigar):
    ''' Enable CIGAR insertions to be pushed leftward through deletions. '''
    cdef char I = 1
    cdef char D = 2

    cdef int cig_len = len(cigar)
    cdef int i, j, dels, del_idx, inss, ins_idx
    for i in range(cig_len-1):
        if cigar[i] == D and cigar[i+1] == I:

            # count adjacent deletions
            del_idx = i-1
            while del_idx >= 0 and cigar[del_idx] == D:
                del_idx -= 1
            dels = i - del_idx

            # count adjacent insertions
            ins_idx = i+1
            while ins_idx < cig_len and cigar[ins_idx] == I:
                ins_idx += 1
            inss = ins_idx - i - 1

            # overwrite cigar
            for j in range(inss):
                cigar[del_idx+1 + j] = I
            for j in range(dels):
                cigar[del_idx+1 + inss + j] = D

    return cigar



def seq_len(cigar):
    length = 0
    for op in cigar:
        if op in 'SXI=M':
            length += 1
    return length

def ref_len(cigar):
    length = 0
    for op in cigar:
        if op in 'XD=M':
            length += 1
    return length


def bases_to_int(seq):
    int_seq = np.zeros(len(seq), dtype=np.uint8)
    for i in range(len(seq)):
        int_seq[i] = cfg.base_dict[ seq[i] ]
    return int_seq

def int_to_bases(int_seq):
    return ''.join([cfg.bases[i] for i in int_seq])

def cig_to_int(cig):
    int_cig = np.zeros(len(cig), dtype=np.uint8)
    for i in range(len(cig)):
        int_cig[i] = cfg.cigar_dict[ cig[i] ]
    return int_cig

def int_to_cig(int_cig):
    return ''.join([cfg.cigars[i] for i in int_cig])

@cython.boundscheck(False)
@cython.wraparound(False)
cdef same_cigar(char[::1] cig1, char[::1] cig2):
    cdef int cig1_len = len(cig1)
    cdef int cig2_len = len(cig2)
    cdef int i

    if cig1_len != cig2_len:
        return False
    else:
        for i in range(cig1_len):
            if cig1[i] != cig2[i]:
                return False
        return True


cpdef standardize_cigar(read_data):
    ''' Try to force all INDELs into beginning of repetitive region. '''
    cdef char I = 1
    cdef char D = 2

    read_id, ref_name, start, stop, cigar, ref, seq, hap = read_data

    # can optionally only report INDELs, assume substitutions found already
    if cfg.args.indel_cigar:
        cigar = cigar.replace('X', 'DI').replace('=','M')
    else: 
        cigar = cigar.replace('X', 'M').replace('=','M')

    cig_len = len(cigar)
    nshifts_buf_arr = np.zeros(cig_len, dtype = np.uint8)
    cdef char[::1] nshifts_buf = nshifts_buf_arr
    shiftlen_buf_arr = np.zeros(cig_len, dtype = np.uint8)
    cdef char[::1] shiftlen_buf = shiftlen_buf_arr
    cdef char[::1] int_cig = cig_to_int(cigar)
    cdef char[::1] int_ref = bases_to_int(ref)
    cdef char[::1] int_seq = bases_to_int(seq)

    while True:
        old_cig = int_cig[:]

        # push left, through
        int_cig = push_indels_left(int_cig, int_ref, nshifts_buf, shiftlen_buf, D)
        int_cig = push_inss_thru_dels(int_cig)
        int_cig = push_indels_left(int_cig, int_seq, nshifts_buf, shiftlen_buf, I)
        int_cig = push_inss_thru_dels(int_cig)

        if same_cigar(old_cig, int_cig): break

    new_cigar = int_to_cig(int_cig).replace('ID','M')

    with cfg.counter.get_lock():
        cfg.counter.value += 1
        print(f"\r    {cfg.counter.value} CIGARs standardized.", end='', flush=True)

    return (read_id, ref_name, start, stop, new_cigar, ref, seq, hap)
