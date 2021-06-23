import numpy as np
np.set_printoptions(linewidth=200)
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

    # confusion matrix
    max_hp = 20
    plt.figure(figsize=(15,15))
    plt.matshow(hps[:max_hp,:max_hp], cmap=plt.cm.Reds, alpha=0.5)
    for i in range(max_hp):
        for j in range(max_hp):
            plt.text(x=j, y=i, s=f'{hps[j,i]:.1f}', fontsize=7, va='center', ha='center')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title('Homopolymer Score Matrix')
    plt.savefig(f'{cfg.args.stats_dir}/{prefix}_scores.png', dpi=300)
    plt.close()

    # surface plot
    x, y = np.meshgrid(range(cfg.args.max_hp), range(cfg.args.max_hp))
    fig = plt.figure(figsize=(20,10))
    ax1 = fig.add_subplot(1,2,1, projection='3d')
    ax2 = fig.add_subplot(1,2,2, projection='3d')
    ax1.plot_trisurf(x.flatten(), -y.flatten(), hps.flatten(), cmap='RdYlGn_r', edgecolor='none')
    ax2.plot_trisurf(-x.flatten(), -y.flatten(), hps.flatten(), cmap='RdYlGn_r', edgecolor='none')
    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/{prefix}_surface.png', dpi=200)
    plt.close()

    # contour plot
    plt.figure(figsize=(15,15))
    plot = plt.contour(hps, levels=10)
    plt.xlabel('Homopolymer Length Called')
    plt.ylabel('Actual Homopolymer Length')
    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/{prefix}_contour.png', dpi=200)
    plt.close()



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
def get_inss_dels(cigar):
    ''' CIGAR must contain only "I" and "D". '''
    inss = np.zeros(len(cigar) + 1)
    dels = np.zeros(len(cigar) + 1)

    for i in range(len(cigar)):

        # verify CIGAR contains only INSs and DELs
        is_ins = cigar[i] == 'I'
        is_del = cigar[i] == 'D'
        if not is_ins and not is_del:
            print("ERROR: unexpected CIGAR type during 'get_inss_dels()'.")

        # create DP arrays for fast lookups
        if is_ins:
            inss[i+1] = inss[i] + 1
        else:
            inss[i+1] = inss[i]

        if is_del:
            dels[i+1] = dels[i] + 1
        else:
            dels[i+1] = dels[i]
    return inss, dels



@njit()
def a_to_b(a_row, a_col, inss, dels, r):
    b_row = int(a_row + a_col)
    b_col = int(inss[b_row] - a_row + r)
    return b_row, b_col



@njit()
def b_to_a(b_row, b_col, inss, dels, r):
    return (int(inss[int(b_row)] + r - b_col), int(dels[int(b_row)] - r + b_col))



@njit()
def align(ref, seq, orig_ref, cigar, sub_scores, hp_scores, 
        indel_start=5, indel_extend=2, r = 40, verbose=False):
    ''' Perform alignment. '''

    # convert CIGAR so that each movement is row+1 or col+1, enables easy banding
    cigar = expand_cigar(cigar)
    cigar = cigar.replace('X','DI').replace('=','DI') \
            .replace('M','DI').replace('S','')
    inss, dels = get_inss_dels(cigar)

    # set helpful constants
    a_rows = len(seq) + 1
    a_cols = len(ref) + 1
    b_rows = len(seq) + len(ref) + 1
    b_cols = 2*r + 1

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

    # precompute hompolymers
    ref_hp_lens = get_hp_lengths(ref)

    matrix = np.zeros((typs, b_rows, b_cols, dims))
    # calculate matrix
    for b_row in range(b_rows):
        for b_col in range(b_cols):

            # precompute useful positions
            a_row, a_col = b_to_a(b_row, b_col, inss, dels, r)
            ref_idx = a_col - 1
            seq_idx = a_row - 1
            b_top = a_to_b(a_row-1, a_col, inss, dels, r)
            b_left = a_to_b(a_row, a_col-1, inss, dels, r)
            b_diag = a_to_b(a_row-1, a_col-1, inss, dels, r)

            # skip cells out of range of original "A" matrix
            if a_row < 0 or a_col < 0 or a_row >= a_rows or a_col >= a_cols:
                continue

            # very first cell has no score (necessary to skip other elifs)
            elif a_row == 0 and a_col == 0:
                continue

            # initialize first row/col of matrix A
            elif a_row == 0:
                for typ in range(typs):
                    if typ == INS or typ == LHP:
                        del_val = indel_start + (100+indel_extend)*(a_col-1)
                    else:
                        del_val = indel_start + indel_extend*(a_col-1)
                    matrix[typ, b_row, b_col, :] = [del_val, DEL, a_col]
                continue
            elif a_col == 0:
                for typ in range(typs):
                    if typ == DEL or typ == SHP:
                        ins_val = indel_start + (100+indel_extend)*(a_row-1)
                    else:
                        ins_val = indel_start + indel_extend*(a_row-1)
                    matrix[typ, b_row, b_col, :] = [ins_val, INS, a_row]
                continue

            # shouldn't be possible
            elif b_row == 0:
                print("ERROR: b_row == 0")

            # enforce new path remains within r cells of original path
            elif b_col == 0 or b_col == 2*r:
                val = max(matrix[:, b_row-1, b_col, VALUE]) + 100
                for typ in range(typs):
                    matrix[typ, b_row, b_col, :] = [val, SUB, 0]
                continue


            # UPDATE INS MATRIX
            # continue INS
            if a_row == 1:
                ins_run = 1
            else:
                ins_run = matrix[INS, b_top[0], b_top[1], RUNLEN] + 1
            ins_val = matrix[INS, b_top[0], b_top[1], VALUE] + indel_extend
            matrix[INS, b_row, b_col, :] = [ins_val, INS, ins_run]

            # start INS
            start_val = matrix[SUB, b_top[0], b_top[1], VALUE] + indel_start
            if start_val < ins_val:
                matrix[INS, b_row, b_col, :] = [start_val, INS, 1]


            # UPDATE DEL MATRIX
            # continue DEL
            if a_col == 1:
                del_run = 1
            else:
                del_run = matrix[DEL, b_left[0], b_left[1], RUNLEN] + 1
            del_val = matrix[DEL, b_left[0], b_left[1], VALUE] + indel_extend
            matrix[DEL, b_row, b_col, :] = [del_val, DEL, del_run]

            # start DEL
            start_val = matrix[SUB, b_left[0], b_left[1], VALUE] + indel_start
            if start_val < del_val:
                matrix[DEL, b_row, b_col, :] = [start_val, DEL, 1]


            # UPDATE LHP MATRIX
            if seq_idx+2 < len(seq) and ref_idx+1 < len(ref) and \
                    seq[seq_idx+1] == seq[seq_idx+2]: # only insert same base

                # continue LHP
                if a_row == 1:
                    lhp_run = 1
                else:
                    lhp_run = int(matrix[LHP, b_top[0], b_top[1], RUNLEN] + 1)
                lhp_val = matrix[LHP, a_to_b(a_row-lhp_run, a_col, inss, dels, r)[0], 
                        a_to_b(a_row-lhp_run, a_col, inss, dels, r)[1], VALUE] + \
                        hp_indel_score(ref_hp_lens[ref_idx+1], lhp_run, hp_scores)
                matrix[LHP, b_row, b_col, :] = [lhp_val, LHP, lhp_run]

                # start LHP
                if seq[seq_idx+1] != seq[seq_idx]:
                    start_val = matrix[SUB, b_top[0], b_top[1], VALUE] + \
                            hp_indel_score(ref_hp_lens[ref_idx+1], 1, hp_scores)
                    if start_val < lhp_val:
                        matrix[LHP, b_row, b_col, :] = [start_val, LHP, 1]

            else: # don't allow insertion of different base
                lhp_val = matrix[SUB, b_diag[0], b_diag[1], VALUE] + 100
                matrix[LHP, b_row, b_col, :] = [lhp_val, LHP, 0]


            # UPDATE SHP MATRIX
            if ref_idx+1 < len(ref) and \
                    ref[ref_idx+1] == ref[ref_idx]: # only delete same base

                # continue SHP
                if a_col == 1:
                    shp_run = 1
                else:
                    shp_run = int(matrix[SHP, b_left[0], b_left[1], RUNLEN] + 1)
                shp_val = matrix[SHP, a_to_b(a_row, a_col-shp_run, inss, dels, r)[0], \
                        a_to_b(a_row, a_col-shp_run, inss, dels, r)[1], VALUE] + \
                        hp_indel_score(ref_hp_lens[ref_idx], -shp_run, hp_scores)
                matrix[SHP, b_row, b_col, :] = [shp_val, SHP, shp_run]

                # start SHP
                if ref_idx and ref[ref_idx] != ref[ref_idx-1]:
                    start_val = matrix[SUB, b_left[0], b_left[1], VALUE] + \
                        hp_indel_score(ref_hp_lens[ref_idx], -1, hp_scores)
                    if start_val < shp_val:
                        matrix[SHP, b_row, b_col, :] = [start_val, SHP, 1]

            else: # don't allow deleting different base
                shp_val = matrix[SUB, b_diag[0], b_diag[1], VALUE] + 100
                matrix[SHP, b_row, b_col, :] = [shp_val, SHP, 0]


            # UPDATE SUB MATRIX
            if matrix[SUB, b_diag[0], b_diag[1], TYPE] == SUB:
                runlen = matrix[SUB, b_diag[0], b_diag[1], RUNLEN] + 1
            else:
                runlen = 1
            sub_val = matrix[SUB, b_diag[0], b_diag[1], VALUE] + \
                    sub_scores[ 
                            base_idx( seq[seq_idx] ), 
                            base_idx( ref[ref_idx] ) 
                    ]
            min_val = sub_val
            matrix[SUB, b_row, b_col, :] = [min_val, SUB, runlen]

            # end INDEL
            for typ in [INS, LHP, DEL, SHP]:
                end_val = matrix[typ, b_row, b_col, VALUE]
                if end_val < min_val:
                    min_val = end_val
                    runlen = matrix[typ, b_row, b_col, RUNLEN]
                    matrix[SUB, b_row, b_col, :] = [min_val, typ, runlen]


    # initialize backtracking from last cell
    aln = ''
    path = []
    row, col = a_rows-1, a_cols-1
    b_pos = a_to_b(row, col, inss, dels, r)
    typ = matrix[SUB, b_pos[0], b_pos[1], VALUE]

    # backtrack
    while row > 0 or col > 0:
        if row < 0:
            print("\nERROR: row < 0 @  A: (", row, ",", col, "), B: (", b_pos[0], ",", b_pos[1], ").")
            print(matrix[SUB,:10,:,VALUE])
            print(matrix[SUB,:10,:,TYPE])
            print(matrix[SUB,:10,:,RUNLEN])
            break
        if col < 0:
            print("\nERROR: col < 0 @  A: (", row, ",", col, "), B: (", b_pos[0], ",", b_pos[1], ").")
            print(matrix[SUB,:10,:,VALUE])
            print(matrix[SUB,:10,:,TYPE])
            print(matrix[SUB,:10,:,RUNLEN])
            break
        path.append((SUB, row, col))
        b_pos = a_to_b(row, col, inss, dels, r)
        val, typ, runlen = matrix[SUB, b_pos[0], b_pos[1], :]
        # Invalid SHPs and LHPs are runlen 0. They should never be reached during 
        # backtracking (due to high penalties), but leaving this here just-in-case
        runlen = max(1, int(runlen))
        if runlen == 0:
            print("\nERROR: RUNLEN == 0 @  A: (", row, ",", col, "), B: (", b_pos[0], ",", b_pos[1], ").")
            print("REF:", orig_ref[col-20:col], orig_ref[col], orig_ref[col+1:col+20])
            print("SEQ:", seq[row-20:row], seq[row], seq[row+1:row+20])
            print(matrix[SUB,b_pos[0]-5:b_pos[0]+5,:,VALUE])
            print(matrix[SUB,b_pos[0]-5:b_pos[0]+5,:,TYPE])
            print(matrix[SUB,b_pos[0]-5:b_pos[0]+5,:,RUNLEN])
            break
        op = ''
        if typ == LHP or typ == INS:   # each move is an insertion
            for i in range(runlen):
                op += 'I'
            row -= runlen
        elif typ == SHP or typ == DEL: # each move is a deletion
            for i in range(runlen):
                op += 'D'
            col -= runlen
        elif typ == SUB: # only sub if stay in same matrix
            i = 0
            while i < runlen:
                row -= 1
                col -= 1
                op += '=' if orig_ref[col] == seq[row] else 'X'
                i += 1
        else:
            print("ERROR: unknown alignment matrix type:", typ)
            break
        aln += op

    # debug print matrices
    if verbose:
        types = ['SUB', 'INS', 'LHP', 'DEL', 'SHP']
        ops = 'MILDS'

        # PRINT A
        for typ, name in enumerate(types):
            print('\nA:', name)
            s = '    -'
            for base in ref:
                s += '    ' + base
            print(s)

            for row in range(a_rows):

                if row == 0:
                    s = '-'
                else:
                    s = seq[row-1]

                for col in range(a_cols):
                    b_pos = a_to_b(row, col, inss, dels, r)
                    op = ops[int(matrix[typ, b_pos[0], b_pos[1], TYPE])]
                    if (typ, row, col) in path:
                        mark = '$'
                    elif col == 0 or row == 0:
                        mark = '*'
                    elif b_pos[1] == 0 or b_pos[1] == 2*r:
                        mark = '.'
                    else:
                        mark = ' '
                    runlen = int(matrix[typ, b_pos[0], b_pos[1], RUNLEN])
                    if runlen < 10:
                        s += "  " + str(runlen) + op + mark
                    else:
                        s += " " + str(runlen) + op + mark
                print(s)

        # PRINT B
            print('\nB:', name)

            # print columns
            s = ' '
            for col in range(b_cols):
                if col < 10:
                    s = s + '    ' + str(col)
                else:
                    s = s + '   ' + str(col)
            print(s)

            for row in range(b_rows):

                # align labels
                if row < 10:
                    s = str(row) + ' '
                else:
                    s = str(row)

                for col in range(b_cols):
                    a_row, a_col = b_to_a(row, col, inss, dels, r)
                    op = ops[int(matrix[typ, row, col, TYPE])]
                    if (typ, row, col) in path:
                        mark = '$'
                    elif a_row < 0 or a_col < 0 or a_row >= a_rows or a_col >= a_cols:
                        mark = '~'
                    elif a_col == 0 or a_row == 0:
                        mark = '*'
                    elif col == 0 or col == 2*r:
                        mark = '.'
                    else:
                        mark = ' '
                    runlen = int(matrix[typ, row, col, RUNLEN])
                    if runlen < 10:
                        s += "  " + str(runlen) + op + mark
                    else:
                        s += " " + str(runlen) + op + mark
                print(s)
            print(' ')


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

    print(f"REF: len: {len(ref)} ciglen: {sum([op in 'XD=M' for op in cigar])}\n"
          f"SEQ: len: {len(seq)} ciglen: {sum([op in 'SXI=M' for op in cigar])}\n"
          f"Cigar: {cigar}\n\n")

    for x in range(0, len(cig_str), 80):
        print(ref_str[x : x+80])
        print(cig_str[x : x+80])
        print(seq_str[x : x+80])
        print(' ')

