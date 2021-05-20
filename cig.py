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

    cigars = dict()
    starts = dict()
    bam = pysam.AlignmentFile(cfg.args.bam, 'rb')

    for read in bam.fetch():
        cigars[read.query_name] = read.cigartuples
        starts[read.query_name] = read.reference_start

    return cigars, starts



def extend_pysam_cigar(ops, counts):
    return ''.join([int(count)*cfg.cigar[op] for (count, op) in zip(counts, ops)])



def collapse_cigar(extended_cigar):
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
