import multiprocessing as mp
from collections import defaultdict
import numpy as np
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt
import os, re, itertools, sys

import pysam
from numba import njit

import cfg


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



@njit()
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
    return ''.join([int(count)*cfg.cigar[op] for (count, op) in zip(counts, ops)])



@njit()
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

            # push deletion as far left as possible (keeping ref same)
            nshifts = 0
            while ref_del_ptr-nshifts > 0 and \
                    ref[ref_del_ptr-nshifts-1] == ref[ref_del_ptr-nshifts-1 + shift_len] and \
                    cigar[cig_del_ptr-nshifts-1] == '=':
                nshifts += 1
                diff = True
            
            # update CIGAR, try shorter prefix
            # print(cigar)
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



@njit()
def push_inss_left(cigar, seq):
    ''' Push CIGAR insertions leftwards. '''

    seq_ptr, cig_ptr = 0, 0
    diff = False
    while cig_ptr < len(cigar):

        # get insertion length
        op = cigar[cig_ptr]
        if op == 'I':
            ins_len = 1
            while cig_ptr+ins_len < len(cigar) and \
                    cigar[cig_ptr+ins_len] == 'I':
                ins_len += 1
        else:
            ins_len = 0

        # iterate, pushing shorter prefixes of insertions left
        shift_len = ins_len
        seq_ins_ptr = seq_ptr
        cig_ins_ptr = cig_ptr

        # while not at CIGAR start and we still have insertions to push
        while cig_ins_ptr > shift_len and shift_len > 0:

            # push insertion as far left as possible (keeping seq same)
            nshifts = 0
            while seq_ins_ptr - nshifts*shift_len >= 0 and \
                    seq[seq_ins_ptr - nshifts*shift_len : seq_ins_ptr - (nshifts-1)*shift_len] == \
                    seq[seq_ins_ptr : seq_ins_ptr + shift_len] and \
                    'D' not in cigar[cig_ins_ptr - nshifts*shift_len : \
                            cig_ins_ptr - (nshifts-1)*shift_len] and \
                    'I' not in cigar[cig_ins_ptr - nshifts*shift_len : \
                            cig_ins_ptr - (nshifts-1)*shift_len]:
                nshifts += 1
                diff = True
            
            # update CIGAR, try shorter prefix
            if diff:
                # print(cigar)
                cigar = (
                        cigar[ : cig_ins_ptr-nshifts*shift_len] + # prefix
                        cigar[cig_ins_ptr : cig_ins_ptr+shift_len] + # inss (shifted left)
                        cigar[cig_ins_ptr-nshifts*shift_len : cig_ins_ptr] + # equals (shifted right)
                        cigar[cig_ins_ptr+shift_len : ] # suffix
                )

                # print(' '*(cig_ptr) + '| cig_ptr')
                # print(' '*(cig_ins_ptr) + '| cig_ins_ptr')
                # print(cigar, f'{ins_len}I total, {shift_len}I shifted back {nshifts}')
                # print(' '*(seq_ptr) + '| seq_ptr')
                # print(' '*(seq_ins_ptr) + '| seq_ins_ptr')
                # print(seq)
                # print(' ')

            cig_ins_ptr -= nshifts*shift_len
            seq_ins_ptr -= nshifts*shift_len
            shift_len -= 1

        # update pointers
        cig_ptr += max(1, ins_len)
        if op == 'M' or op == 'X' or op == '=':
            seq_ptr += 1
        elif op == 'I':
            seq_ptr += ins_len

    return cigar, diff



def subs_to_indels(cigar):
    return cigar.replace('X', 'DI')



@njit()
def push_inss_thru_dels(cigar):
    ''' Enable CIGAR insertions to be pushed leftward through deletions. '''

    diff = False
    for i in range(len(cigar)-1):
        if cigar[i] == 'D' and cigar[i+1] == 'I':

            # count adjacent deletions
            del_idx = i-1
            while del_idx >= 0 and cigar[del_idx] == 'D':
                del_idx -= 1
            dels = i - del_idx

            # count adjacent insertions
            ins_idx = i+1
            while ins_idx < len(cigar) and cigar[ins_idx] == 'I':
                ins_idx += 1
            inss = ins_idx - i - 1

            cigar = (
                    cigar[:del_idx+1] + inss*'I' + dels*'D' + cigar[ins_idx:]
            )

            diff = True
    return cigar, diff



@njit()
def seq_len(cigar):
    length = 0
    for op in cigar:
        if op in 'SXI=M':
            length += 1
    return length

@njit()
def ref_len(cigar):
    length = 0
    for op in cigar:
        if op in 'XD=M':
            length += 1
    return length



def standardize_cigar(read_data):
    ''' Try to force all INDELs into beginning of repetitive region. '''

    read_id, ref_name, start, cigar, ref, seq = read_data

    cigar0 = subs_to_indels(cigar)

    diff = True
    while diff: # loop until CIGAR is stable

        cigar1, diff1 = push_dels_left(cigar0, ref)

        cigar2, diff2 = push_inss_left(cigar1, seq)

        cigar3, diff3 = push_inss_thru_dels(cigar2)

        diff = diff1 or diff2 or diff3
        cigar0 = cigar3

    final_cigar = collapse_cigar(cigar0)

    return (read_id, ref_name, start, final_cigar, ref, seq)

