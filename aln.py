#!/usr/bin/env python
'''
Simple Smith-Waterman aligner
'''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from numba import njit
import scipy.ndimage as ndimage

import cfg
from cig import *


@njit()
def fix_matrix_properties(scores, delta = 0.01):
    ''' Modify matrix so scores follow expected pattern. '''

    # don't penalize diagonals
    l = scores.shape[0]
    for i in range(1, l):
        scores[i,i] = min(scores[i,i], scores[i-1, i-1])

    # more insertions should be more penalized
    for j in range(1, l):
        for i in range(j-1, -1, -1):
            scores[i,j] = max(
                    scores[i,j], 
                    scores[i+1,j] + delta, 
                    scores[i,j-1] + delta
            )

    # more deletions should be more penalized
    for i in range(1,l):
        for j in range(i-1, -1, -1):
            scores[i,j] = max(
                    scores[i,j], 
                    scores[i,j+1] + delta, 
                    scores[i-1,j] + delta
            )

    # prefer insertions from longer homopolymers
    best = np.ones(l) * 1000
    for j in range(1,l):
        for i in range(j-1, -1, -1):
            ins_len = j - i
            if scores[i,j] < best[ins_len]:
                best[ins_len] = scores[i,j]
                for total_ins_len in range(ins_len+1, l):
                    best[total_ins_len] = min(
                            best[total_ins_len], 
                            best[ins_len] + best[total_ins_len-ins_len]
                    )
            else:
                scores[i,j] = min(
                        scores[i, j], 
                        best[ins_len] - delta
                )

    # prefer deletions from longer homopolymers
    best = np.ones(l) * 1000
    for i in range(1,l):
        for j in range(i-1, -1, -1):
            del_len = i - j
            if scores[i,j] < best[del_len]:
                best[del_len] = scores[i,j]
                for total_del_len in range(del_len+1, l):
                    best[total_del_len] = min(
                            best[total_del_len], 
                            best[del_len] + best[total_del_len-del_len]
                    )
            else:
                scores[i,j] = min(
                        scores[i, j], 
                        best[del_len] - delta
                )

    return scores



def calc_score_matrices(subs, hps):

    # calculate homopolymer scores matrix
    hps = np.sum(hps, axis=0)
    hp_scores = np.zeros_like(hps)
    for ref_len in range(cfg.args.max_hp):
        total = np.sum(hps[ref_len])
        for call_len in range(cfg.args.max_hp):
            count = int(hps[ref_len, call_len])
            bias = 10
            frac = (count + 0.01 + int(ref_len==call_len)*bias) / (total + 0.01*cfg.args.max_hp + bias)
            hp_scores[ref_len, call_len] = -np.log(frac)
    # hp_scores = ndimage.gaussian_filter(hp_scores, sigma=(2,2))
    hp_scores = fix_matrix_properties(hp_scores)

    # calculate substitution scores matrix
    sub_scores = np.zeros_like(subs)
    for i in range(len(cfg.bases)):
        total = np.sum(subs[i])
        for j in range(len(cfg.bases)):
            if i != j:
                sub_scores[i, j] = -np.log( (subs[i,j]+0.1) / total )
            else:
                sub_scores[i, j] = 0

    return sub_scores, hp_scores



def plot_hp_score_matrix(hps, prefix="score_mat"):

    # surface plot
    x, y = np.meshgrid(range(cfg.args.max_hp), range(cfg.args.max_hp))
    fig = plt.figure(figsize=(20,10))
    ax1 = fig.add_subplot(1,2,1, projection='3d')
    ax2 = fig.add_subplot(1,2,2, projection='3d')
    ax1.plot_trisurf(x.flatten(), -y.flatten(), hps.flatten(), cmap='RdYlGn_r', edgecolor='none')
    ax2.plot_trisurf(-x.flatten(), -y.flatten(), hps.flatten(), cmap='RdYlGn_r', edgecolor='none')
    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/{prefix}_surface.png', dpi=200)

    # contour plot
    plt.figure(figsize=(15,15))
    plot = plt.contour(hps, levels=10)
    plt.xlabel('Homopolymer Length Called')
    plt.ylabel('Actual Homopolymer Length')
    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/{prefix}_contour.png', dpi=200)



@njit()
def get_hp_lengths(seq):
    ''' Calculate HP length of substring starting at each index. '''

    # TODO: make this more efficient
    hp_lens = [0]*(len(seq))
    for start in range(len(seq)):
        for stop in range(start+1, len(seq)):
            if seq[stop] != seq[start]:
                hp_lens[start] = stop - start
                break
    if len(seq):
        hp_lens[-1] += 1
    return hp_lens



@njit()
def hp_indel_score(ref_hp_len, indel_len, hp_scores):

    # error, don't allow
    if ref_hp_len <= 0:
        return 100
    elif ref_hp_len + indel_len < 0:
        return 100

    # force lengths to fit in matrix
    ref_hp_len = min(ref_hp_len, len(hp_scores)-1)
    call_hp_len = min(ref_hp_len+indel_len, len(hp_scores)-1)
    return hp_scores[ref_hp_len, call_hp_len]


@njit()
def base_idx(base):
    if base == 'A':
        return 0
    elif base == 'C':
        return 1
    elif base == 'G':
        return 2
    elif base == 'T':
        return 3
    else:
        return -1



@njit()
def align(ref, query, orig_ref, sub_scores, hp_scores, indel_start=5, indel_extend=2, verbose=False):
    ''' Perform alignment. '''

    ref_hp_lens = get_hp_lengths(ref)

    rows = len(query) + 1
    cols = len(ref) + 1

    dims = 3
    VALUE = 0
    TYPE = 1
    RUNLEN = 2

    typs = 5# types
    SUB = 0 # substitution
    INS = 1 # insertion
    LHP = 2 # lengthen homopolymer
    DEL = 3 # deletion
    SHP = 4 # shorten homopolymer

    # initialize first row/col of matrix
    matrix = np.zeros((typs, rows, cols, dims))
    for typ in range(typs):
        for row in range(1, rows):
            ins_val = indel_start if row == 1 else \
                    matrix[typ, row-1, 0, VALUE] + indel_extend
            matrix[typ, row, 0, :] = [ins_val, INS, 0]

    for typ in range(typs):
        for col in range(1, cols):
            del_val = indel_start if col == 1 else \
                    matrix[typ, 0, col-1, VALUE] + indel_extend
            matrix[typ, 0, col, :] = [del_val, DEL, 0]

    # calculate matrix
    for row in range(1, rows):
        for col in range(1, cols):
            ref_idx = col - 1
            query_idx = row - 1

            # UPDATE INS MATRIX
            # continue INS
            ins_run = matrix[INS, row-1, col, RUNLEN] + 1
            ins_val = matrix[INS, row-1, col, VALUE] + indel_extend
            matrix[INS, row, col, :] = [ins_val, INS, ins_run]

            # start INS
            min_val = ins_val
            for typ in [SUB, LHP, SHP]:
                start_val = matrix[typ, row-1, col, VALUE] + indel_start
                if start_val < min_val:
                    min_val = start_val
                    matrix[INS, row, col, :] = [min_val, typ, 1]


            # UPDATE DEL MATRIX
            # continue DEL
            del_run = matrix[DEL, row, col-1, RUNLEN] + 1
            del_val = matrix[DEL, row, col-1, VALUE] + indel_extend
            matrix[DEL, row, col, :] = [del_val, DEL, del_run]

            # start DEL
            min_val = del_val
            for typ in [SUB, LHP, SHP]:
                start_val = matrix[typ, row, col-1, VALUE] + indel_start
                if start_val < min_val:
                    min_val = start_val
                    matrix[DEL, row, col, :] = [min_val, typ, 1]


            # UPDATE LHP MATRIX
            if query_idx+2 < len(query) and ref_idx+1 < len(ref) and \
                    query[query_idx+1] == query[query_idx+2]: # only insert same base

                # continue LHP
                lhp_run = int(matrix[LHP, row-1, col, RUNLEN] + 1)
                lhp_val = matrix[LHP, row-lhp_run, col, VALUE] + \
                        hp_indel_score(ref_hp_lens[ref_idx+1], lhp_run, hp_scores)
                matrix[LHP, row, col, :] = [lhp_val, LHP, lhp_run]

                # start LHP
                # if query[query_idx+1] != query[query_idx]:
                min_val = lhp_val
                for typ in [SUB, DEL, INS, SHP]:
                    start_val = matrix[typ, row-1, col, VALUE] + \
                            hp_indel_score(ref_hp_lens[ref_idx+1], 1, hp_scores)
                    if start_val < min_val:
                        min_val = start_val
                        matrix[LHP, row, col, :] = [min_val, typ, 1]

            else: # don't allow insertion of different base
                lhp_val = max(matrix[:, row-1, col, VALUE]) + 100
                matrix[LHP, row, col, :] = [lhp_val, SUB, 0]


            # UPDATE SHP MATRIX
            if ref_idx+1 < len(ref) and \
                    ref[ref_idx+1] == ref[ref_idx]: # only delete same base

                # continue SHP
                shp_run = int(matrix[SHP, row, col-1, RUNLEN] + 1)
                shp_val = matrix[SHP, row, col-shp_run, VALUE] + \
                        hp_indel_score(ref_hp_lens[ref_idx], -shp_run, hp_scores)
                matrix[SHP, row, col, :] = [shp_val, SHP, shp_run]

                # start SHP
                # if ref_idx and ref[ref_idx] != ref[ref_idx-1]:
                min_val = shp_val
                for typ in [SUB, DEL, INS, LHP]:
                    start_val = matrix[typ, row, col-1, VALUE] + \
                        hp_indel_score(ref_hp_lens[ref_idx], -1, hp_scores)
                    if start_val < min_val:
                        min_val = start_val
                        matrix[SHP, row, col, :] = [min_val, typ, 1]

            else: # don't allow deleting different base
                shp_val = max(matrix[:, row, col-1, VALUE]) + 100
                matrix[SHP, row, col, :] = [shp_val, SUB, 0]


            # UPDATE SUB MATRIX
            # simple SUB lookup
            sub_val = matrix[SUB, row-1, col-1, VALUE] + \
                    sub_scores[ 
                            base_idx( query[query_idx] ), 
                            base_idx( ref[ref_idx] ) 
                    ]
            min_val = sub_val
            matrix[SUB, row, col, :] = [min_val, SUB, 0]

            # end INDEL
            for typ in [INS, LHP, DEL, SHP]:
                end_val = matrix[typ, row, col, VALUE]
                if end_val < min_val:
                    min_val = end_val
                    matrix[SUB, row, col, :] = [min_val, typ, 0]


    # initialize backtracking from last cell
    aln = ''
    path = []
    row, col = rows-1, cols-1
    old_typ = np.argmin(matrix[:, row, col, VALUE])

    # backtrack
    while row > 0 or col > 0:
        path.append((int(old_typ), row, col))
        val, new_typ, runlen = matrix[int(old_typ), row, col, :]
        op = ''
        if new_typ == LHP or new_typ == INS:   # each move is an insertion
            op = 'I'
            row -= 1
        elif new_typ == SHP or new_typ == DEL: # each move is a deletion
            op = 'D'
            col -= 1
        elif new_typ == SUB: # only sub if stay in same matrix
            if old_typ == SUB:
                row -= 1
                col -= 1
                op = '=' if orig_ref[col] == query[row] else 'X'
        else:
            print("ERROR: unknown alignment matrix type '" + str(typ) + "'.")
        aln += op
        old_typ = new_typ

    # debug print matrices
    if verbose:
        types = ['SUB', 'INS', 'LHP', 'DEL', 'SHP']
        ops = 'MILDS'
        for typ, name in enumerate(types):
            print('\n\n', name)
            print('  -    ')
            print('    '.join(ref))
            print('\n')
            for row in range(rows):
                if row == 0:
                    print('-')
                else:
                    print(query[row-1])

                for col in range(cols):
                    print(" " + str(int(matrix[typ, row, col, RUNLEN])) + \
                            ops[int(matrix[typ, row, col, TYPE])] + \
                            '$' if (typ, row, col) in path else ' ')
                print('\n')

    # we backtracked, so get forward alignment
    return aln[::-1]



def dump(ref, seq, cigar):
    ''' Pretty print full alignment result. '''

    ref_str = ''
    cig_str = ''
    seq_str = ''

    ref_idx = 0
    seq_idx = 0

    for idx, op in enumerate(cigar):
        if op == '=':
            ref_str += ref[ref_idx]
            ref_idx += 1
            seq_str += seq[seq_idx]
            seq_idx += 1
            cig_str += '|'

        elif op == 'X':
            ref_str += ref[ref_idx]
            ref_idx += 1
            seq_str += seq[seq_idx]
            seq_idx += 1
            cig_str += 'X'

        elif op == 'D':
            ref_str += ref[ref_idx]
            ref_idx += 1
            seq_str += '-'
            cig_str += ' '

        elif op == 'I':
            ref_str += '-'
            seq_str += seq[seq_idx]
            seq_idx += 1
            cig_str += ' '

        else:
            print(f"ERROR: unrecognized CIGAR operation '{op}' at cigar index {len(cig_str)}.")
            exit(1)

    print(f"REF: {ref_str} len: {len(ref)} ciglen: {sum([op in 'XD=M' for op in cigar])}\n"
          f"     {cig_str}\n"
          f"SEQ: {seq_str} len: {len(seq)} ciglen: {sum([op in 'SXI=M' for op in cigar])}\n"
          f"Cigar: {collapse_cigar(cigar)}\n\n")
