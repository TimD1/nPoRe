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


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef char[::1] bases_to_int(str seq):
    cdef long long seqlen = len(seq)
    int_seq_buf = np.zeros(seqlen, dtype=np.uint8)
    cdef char[::1] int_seq = int_seq_buf
    cdef long long i
    for i in range(seqlen):
        if seq[i] == 'N':
            int_seq[i] = 0
        elif seq[i] == 'A':
            int_seq[i] = 1
        elif seq[i] == 'C':
            int_seq[i] = 2
        elif seq[i] == 'G':
            int_seq[i] = 3
        elif seq[i] == 'T':
            int_seq[i] = 4
        elif seq[i] == '-':
            int_seq[i] = 5
    return int_seq

def int_to_bases(int_seq):
    return ''.join([cfg.bases[i] for i in int_seq])
