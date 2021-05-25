import multiprocessing as mp
from collections import defaultdict
import numpy as np
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt
import os, re, itertools, sys

import pysam

import cfg
from aln import *

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
    return ''.join([int(count)*cfg.cigar[op] for (count, op) in zip(counts, ops)])



def collapse_cigar(extended_cigar):
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

    out = ''
    for num, op in groups:
        out += '%s%s' % (num, op)
    return out



def push_dels_left(cigar, ref, seq):


    ref_ptr, seq_ptr, cig_ptr = 0, 0, 0
    while cig_ptr < len(cigar):

        # get deletion length
        op = cigar[cig_ptr]
        if op == 'D':
            del_len = 1
            while cig_ptr + del_len < len(cigar) and \
                    cigar[cig_ptr + del_len] == 'D':
                del_len += 1
        else:
            del_len = 0

        # iterate, pushing shorter prefixes of deletion left
        while del_len > 0:


        # update pointers
        if op == 'M' or op == 'X' or op == '=':
            ref_ptr = ref_ptr + 1
            seq_ptr = seq_ptr + 1
        elif op == 'D':
            ref_ptr = ref_ptr + 1
        elif op == 'I':
            seq_ptr = seq_ptr + 1




def standardize_cigar(cigar, ref, seq):
    ''' Try to force all INDELs into beginning of repetitive region. '''

    diff = True

    while diff: # loop until CIGAR is stable

        cigar, diff = push_dels_left(cigar, ref, seq)

