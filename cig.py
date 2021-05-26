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



def push_dels_left(cigar, ref):
    ''' Push CIGAR deletions leftwards. '''


    ref_ptr, cig_ptr = 0, 0
    diff = False
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
        shift_len = del_len
        ref_del_ptr = ref_ptr
        cig_del_ptr = cig_ptr

        # while not at CIGAR start and we still have deletions to push
        while cig_del_ptr > 0 and shift_len > 0:

            # push deletion as far left as possible (keeping ref sequence same)
            nshifts = 0
            while ref_del_ptr-nshifts > 0 and \
                    ref[ref_del_ptr-nshifts-1] == ref[ref_del_ptr-nshifts-1 + shift_len] and \
                    cigar[cig_del_ptr-nshifts-1] == '=':
                nshifts += 1
                diff = True
            
            # update CIGAR, try shorter prefix
            print(cigar)
            cigar = (
                    cigar[                      : cig_del_ptr-nshifts] +   # prefix
                    cigar[cig_del_ptr           : cig_del_ptr+shift_len] + # dels (shifted left)
                    cigar[cig_del_ptr-nshifts   : cig_del_ptr] +           # equals (shifted right)
                    cigar[cig_del_ptr+shift_len : ]                        # suffix
            )

            # print(' '*(cig_ptr) + '| cig_ptr')
            # print(' '*(cig_del_ptr) + '| cig_del_ptr')
            # print(cigar, f'{del_len}D total, {shift_len}D shifted back {nshifts}')
            # print(' '*(ref_ptr) + '| ref_ptr')
            # print(' '*(ref_del_ptr) + '| ref_del_ptr')
            # print(ref)
            # print(' ')

            cig_del_ptr -= nshifts
            ref_del_ptr -= nshifts
            shift_len -= 1

        # update pointers
        cig_ptr += max(1, del_len)
        if op == 'M' or op == 'X' or op == '=':
            ref_ptr += 1
        elif op == 'D':
            ref_ptr += del_len

    return cigar, diff




def standardize_cigar(cigar, ref, seq):
    ''' Try to force all INDELs into beginning of repetitive region. '''

    diff = True

    while diff: # loop until CIGAR is stable

        cigar, diff = push_dels_left(cigar, ref)

    return cigar

