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
cdef push_indels_left(char[::1] cigar, char[::1] seq, char push_op):
    ''' Push CIGAR indels leftwards. '''
    cdef char M = 0
    cdef char E = 7
    cdef char X = 8

    cdef int diff = 0
    cdef int seq_ptr = 0
    cdef int cig_ptr = 0
    cdef int cig_len = len(cigar)
    cdef int i, seq_indel_ptr, cig_indel_ptr, indel_len, shift_len, nshifts
    cdef char op

    nshifts_buf_arr = np.zeros(cig_len, dtype = np.uint8)
    cdef char[::1] nshifts_buf = nshifts_buf_arr
    shiftlen_buf_arr = np.zeros(cig_len, dtype = np.uint8)
    cdef char[::1] shiftlen_buf = shiftlen_buf_arr

    while cig_ptr < cig_len:

        # get indel length
        op = cigar[cig_ptr]
        if op == push_op:
            indel_len = 1
            while cig_ptr + indel_len < cig_len and \
                    cigar[cig_ptr + indel_len] == push_op:
                indel_len += 1
        else:
            indel_len = 0

        # iterate, pushing shorter prefixes of indel left
        shift_len = indel_len
        seq_indel_ptr = seq_ptr
        cig_indel_ptr = cig_ptr

        # while not at CIGAR start and we still have indels to push
        while cig_indel_ptr > 0 and shift_len > 0:

            # push indel as far left as possible (keeping seq same)
            nshifts = 0
            while seq_indel_ptr-nshifts > 0 and \
                    seq[seq_indel_ptr-nshifts-1] == seq[seq_indel_ptr-nshifts-1 + shift_len] and \
                    (cigar[cig_indel_ptr-nshifts-1] == E or cigar[cig_indel_ptr-nshifts-1] == M) :
                nshifts += 1
                diff = 1
            
            # update CIGAR, try shorter prefix
            # print(cigar)

            # fill buffers (can't fully update in-place)
            for i in range(nshifts):
                nshifts_buf[i] = cigar[cig_indel_ptr-nshifts+i]
            for i in range(shift_len):
                shiftlen_buf[i] = cigar[cig_indel_ptr+i]

            # update cigar
            for i in range(shift_len):
                cigar[cig_indel_ptr-nshifts+i] = shiftlen_buf[i]
            for i in range(nshifts):
                cigar[cig_indel_ptr-nshifts+shift_len+i] = nshifts_buf[i]

            # print(' '*(cig_ptr) + '| cig_ptr')
            # print(' '*(cig_indel_ptr) + '| cig_indel_ptr')
            # print(cigar, f'{indel_len}{cfg.cigars[push_op'} total,' + \
            #              f'{shift_len}{cfg.cigars[push_op]} shifted back {nshifts}')
            # print(' '*(seq_ptr) + '| seq_ptr')
            # print(' '*(seq_indel_ptr) + '| seq_indel_ptr')
            # print(seq)
            # print(' ')

            cig_indel_ptr -= nshifts
            seq_indel_ptr -= nshifts
            shift_len -= 1

        # update pointers
        if indel_len > 0:
            cig_ptr += indel_len
        else:
            cig_ptr += 1
        if op == M or op == X or op == E:
            seq_ptr += 1
        elif op == push_op:
            seq_ptr += indel_len

    return cigar, diff



@cython.boundscheck(False)
@cython.wraparound(False)
cdef push_inss_thru_dels(char[::1] cigar):
    ''' Enable CIGAR insertions to be pushed leftward through deletions. '''
    cdef char I = 1
    cdef char D = 2

    cdef int diff = 0
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

            diff = 1
    return cigar, diff



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

def cig_to_int(cig):
    int_cig = np.zeros(len(cig), dtype=np.uint8)
    for i in range(len(cig)):
        int_cig[i] = cfg.cigar_dict[ cig[i] ]
    return int_cig

def int_to_cig(int_cig):
    return ''.join([cfg.cigars[i] for i in int_cig])


cpdef standardize_cigar(read_data):
    ''' Try to force all INDELs into beginning of repetitive region. '''
    cdef char I = 1
    cdef char D = 2

    read_id, ref_name, start, cigar, ref, seq = read_data
    cigar = expand_cigar(cigar)

    # can optionally only report INDELs, assume substitutions found already
    if cfg.args.indels_only:
        cigar = cigar.replace('X', 'DI').replace('=','M')
    else: 
        cigar = cigar.replace('X', 'M').replace('=','M')

    cdef char[::1] int_cig = cig_to_int(cigar)
    cdef char[::1] int_ref = bases_to_int(ref)
    cdef char[::1] int_seq = bases_to_int(seq)

    cdef int diff = 1
    while diff > 0: # loop until CIGAR is stable
        int_cig, diff1 = push_indels_left(int_cig, int_ref, D)
        int_cig, diff2 = push_indels_left(int_cig, int_seq, I)
        int_cig, diff3 = push_inss_thru_dels(int_cig)
        diff = diff1 + diff2 + diff3 # logical OR (if any changed)

    final_cigar = collapse_cigar(int_to_cig(int_cig).replace('ID','M'))

    # print(f'cig read:{read_id:>8}'
    #     f'\tseq:{len(seq)} {seq_len(cigar)}->{seq_len(expand_cigar(final_cigar))}'
    #     f'\tref:{len(ref)} {ref_len(cigar)}->{ref_len(expand_cigar(final_cigar))}')

    with cfg.read_count.get_lock():
        cfg.read_count.value += 1
        print(f"\r        {cfg.read_count.value} CIGARs standardized.", end='', flush=True)

    return (read_id, ref_name, start, final_cigar, ref, seq)



def change_ref(read_cig, hap_cig, ref, read, hap):

    # initialize data and pointers
    read_cig = expand_cigar(read_cig)
    hap_cig = expand_cigar(hap_cig)
    cig_ptr_read = 0
    cig_ptr_hap = 0
    ref_ptr_read = 0
    ref_ptr_hap = 0
    read_ptr = 0
    hap_ptr = 0
    new_cig = ''

    # step through both CIGARs, updating read alignment
    while cig_ptr_read < len(read_cig) and cig_ptr_hap < len(hap_cig):

        consume_read_cig, consume_hap_cig = False, False

        # hap CIGAR '='
        if hap_cig[cig_ptr_hap] == '=' and read_cig[cig_ptr_read] == '=':
            new_cig += '='
            consume_read_cig = True
            consume_hap_cig = True
        elif hap_cig[cig_ptr_hap] == '=' and read_cig[cig_ptr_read] == 'X':
            new_cig += 'X'
            consume_read_cig = True
            consume_hap_cig = True
        elif hap_cig[cig_ptr_hap] == '=' and read_cig[cig_ptr_read] == 'I':
            new_cig += 'I'
            consume_read_cig = True
            consume_hap_cig = False
        elif hap_cig[cig_ptr_hap] == '=' and read_cig[cig_ptr_read] == 'D':
            new_cig += 'D'
            consume_read_cig = True
            consume_hap_cig = True

        # hap CIGAR 'X'
        if hap_cig[cig_ptr_hap] == 'X' and read_cig[cig_ptr_read] == '=':
            new_cig += 'X'
            consume_read_cig = True
            consume_hap_cig = True
        elif hap_cig[cig_ptr_hap] == 'X' and read_cig[cig_ptr_read] == 'X':
            if read[read_ptr] == hap[hap_ptr]:
                new_cig += '='
            else:
                new_cig += 'X'
            consume_read_cig = True
            consume_hap_cig = True
        elif hap_cig[cig_ptr_hap] == 'X' and read_cig[cig_ptr_read] == 'I':
            new_cig += 'I'
            consume_read_cig = True
            consume_hap_cig = False
        elif hap_cig[cig_ptr_hap] == 'X' and read_cig[cig_ptr_read] == 'D':
            new_cig += 'D'
            consume_read_cig = True
            consume_hap_cig = True

        # hap CIGAR 'I'
        if hap_cig[cig_ptr_hap] == 'I' and read_cig[cig_ptr_read] == '=':
            new_cig += 'D'
            consume_read_cig = False
            consume_hap_cig = True
        elif hap_cig[cig_ptr_hap] == 'I' and read_cig[cig_ptr_read] == 'X':
            new_cig += 'D'
            consume_read_cig = False
            consume_hap_cig = True
        elif hap_cig[cig_ptr_hap] == 'I' and read_cig[cig_ptr_read] == 'I':
            if read[read_ptr] == hap[hap_ptr]:
                new_cig += '='
            else:
                new_cig += 'X'
            consume_read_cig = True
            consume_hap_cig = True
        elif hap_cig[cig_ptr_hap] == 'I' and read_cig[cig_ptr_read] == 'D':
            new_cig += 'D'
            consume_read_cig = False
            consume_hap_cig = True

        # hap CIGAR 'D'
        if hap_cig[cig_ptr_hap] == 'D' and read_cig[cig_ptr_read] == '=':
            new_cig += 'I'
            consume_read_cig = True
            consume_hap_cig = True
        elif hap_cig[cig_ptr_hap] == 'D' and read_cig[cig_ptr_read] == 'X':
            new_cig += 'I'
            consume_read_cig = True
            consume_hap_cig = True
        elif hap_cig[cig_ptr_hap] == 'D' and read_cig[cig_ptr_read] == 'I':
            new_cig += 'I'
            consume_read_cig = True
            consume_hap_cig = False
        elif hap_cig[cig_ptr_hap] == 'D' and read_cig[cig_ptr_read] == 'D':
            consume_read_cig = False
            consume_hap_cig = True

        # update all pointers (based on consumed CIGAR operations)
        if consume_read_cig: 
            if read_cig[cig_ptr_read] == '=':
                ref_ptr_read += 1
                read_ptr += 1
            elif read_cig[cig_ptr_read] == 'X':
                ref_ptr_read += 1
                read_ptr += 1
            elif read_cig[cig_ptr_read] == 'I':
                read_ptr += 1
            elif read_cig[cig_ptr_read] == 'D':
                ref_ptr_read += 1
            cig_ptr_read += 1
        if consume_hap_cig: 
            if hap_cig[cig_ptr_hap] == '=':
                ref_ptr_hap += 1
                hap_ptr += 1
            elif hap_cig[cig_ptr_hap] == 'X':
                ref_ptr_hap += 1
                hap_ptr += 1
            elif hap_cig[cig_ptr_hap] == 'I':
                hap_ptr += 1
            elif hap_cig[cig_ptr_hap] == 'D':
                ref_ptr_hap += 1
            cig_ptr_hap += 1


def flip_cigar_basis(cigar):
    cigar = cigar.replace('I', 'd').replace('D','i').upper()
