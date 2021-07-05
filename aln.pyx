import numpy as np
np.set_printoptions(linewidth=200)
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import scipy.ndimage as ndimage

import cython

import cfg

from cig import *


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
    hp_scores = np.zeros_like(hps, dtype=np.float32)
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
    sub_scores = np.zeros((len(cfg.bases),len(cfg.bases)), dtype=np.float32)
    for i in range(1, len(cfg.bases)):
        total = np.sum(subs[i-1])
        for j in range(1, len(cfg.bases)):
            if i != j:
                sub_scores[i, j] = -np.log( (subs[i-1,j-1]+0.1) / total )
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



@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[:] get_hp_lengths(char[:] seq):
    ''' Calculate HP length of substring starting at each index. '''

    cdef int seq_len = len(seq)
    hp_lens_buf = np.zeros(seq_len, dtype=np.intc)
    cdef int[:] hp_lens = hp_lens_buf
    cdef int start, stop

    for start in range(seq_len):
        for stop in range(start+1, seq_len):
            if seq[stop] != seq[start]:
                hp_lens[start] = stop - start
                break
    if seq_len:
        hp_lens[seq_len-1] += 1
    return hp_lens



@cython.boundscheck(False)
@cython.wraparound(False)
cdef float hp_indel_score(int ref_hp_len, int indel_len, float[:,:] hp_scores):

    # error, don't allow
    if ref_hp_len <= 0:
        return 100
    elif ref_hp_len + indel_len < 0:
        return 100

    # force lengths to fit in matrix
    cdef int hp_scores_len = len(hp_scores)
    if ref_hp_len > hp_scores_len-1:
        ref_hp_len = hp_scores_len-1

    cdef int call_hp_len = ref_hp_len + indel_len
    if call_hp_len > hp_scores_len-1:
        call_hp_len = hp_scores_len-1
    return hp_scores[ref_hp_len, call_hp_len]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[:] get_inss(str cigar):
    ''' CIGAR must contain only "I" and "D". '''

    cdef int cig_len = len(cigar)
    inss_buf = np.zeros(cig_len+1, dtype=np.intc)
    cdef int[:] inss = inss_buf
    cdef int i

    for i in range(cig_len):
        if cigar[i] == 'I':
            inss[i+1] = inss[i] + 1
        else:
            inss[i+1] = inss[i]
    return inss



@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[:] get_dels(str cigar):
    ''' CIGAR must contain only "I" and "D". '''

    cdef int cig_len = len(cigar)
    dels_buf = np.zeros(cig_len+1, dtype=np.intc)
    cdef int[:] dels = dels_buf
    cdef int i

    for i in range(cig_len):
        if cigar[i] == 'D':
            dels[i+1] = dels[i] + 1
        else:
            dels[i+1] = dels[i]
    return dels



@cython.boundscheck(False)
@cython.wraparound(False)
cdef int a_to_b0(int a_row, int a_col, int[:] inss, int[:] dels, int r):
    return a_row + a_col

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int a_to_b1(int a_row, int a_col, int[:] inss, int[:] dels, int r):
    cdef int b_row, b_col
    b_row = a_row + a_col
    b_col = inss[b_row] - a_row + r
    return b_col



@cython.boundscheck(False)
@cython.wraparound(False)
cdef int b_to_a0(int b_row, int b_col, int[:] inss, int[:] dels, int r):
    return inss[b_row] + r - b_col

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int b_to_a1(int b_row, int b_col, int[:] inss, int[:] dels, int r):
    return dels[b_row] - r + b_col



@cython.boundscheck(False)
@cython.wraparound(False)
cdef float get_max(float[:,:,:,:] matrix, int row, int col, int typs, int VALUE):
    cdef int typ
    cdef float max_val = matrix[0, row, col, VALUE]
    for typ in range(1, typs):
        if matrix[typ, row, col, VALUE] > max_val:
            max_val = matrix[typ, row, col, VALUE]
    return max_val


# @profile
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef align(char[::1] ref, char[::1] seq, char[::1] orig_ref, str cigar, 
        float[:,::1] sub_scores, float[:,::1] hp_scores, 
        float indel_start=5, float indel_extend=2, 
        int r = 30, int verbosity=0):
    ''' Perform alignment. 
    verbosity:
     0 -> no printing
     1 -> print final matrix
     2 -> print each cell (with pointers and scores)
    '''

    # convert CIGAR so that each movement is row+1 or col+1, enables easy banding
    cigar = expand_cigar(cigar)
    cigar = cigar.replace('X','DI').replace('=','DI') \
            .replace('M','DI').replace('S','')

    cdef int[:] inss = get_inss(cigar)
    cdef int[:] dels = get_dels(cigar)

    # set helpful constants
    cdef int a_rows = len(seq) + 1
    cdef int a_cols = len(ref) + 1
    cdef int b_rows = len(seq) + len(ref) + 1
    cdef int max_rows = 110000
    cdef int b_cols = 2*r + 1

    cdef int dims = 3
    cdef int VALUE = 0
    cdef int TYPE = 1
    cdef int RUNLEN = 2

    cdef int typs = 5# types
    cdef int MAT = 0 # substitution
    cdef int INS = 1 # insertion
    cdef int LHP = 2 # lengthen homopolymer
    cdef int DEL = 3 # deletion
    cdef int SHP = 4 # shorten homopolymer

    # precompute hompolymers
    cdef int[:] ref_hp_lens = get_hp_lengths(ref)

    matrix_buf = np.zeros((typs, max_rows, b_cols, dims), dtype=np.float32)
    cdef float[:,:,:,::1] matrix = matrix_buf

    # calculate matrix
    cdef int b_row, b_col, a_row, a_col, ref_idx, seq_idx
    cdef int b_top0, b_top1, b_left0, b_left1, b_diag0, b_diag1
    cdef int typ, i
    cdef float start_val, end_val
    cdef float min_val, max_val, ins_val, del_val, shp_val, lhp_val, sub_val
    cdef int ins_run, del_run, shp_run, lhp_run, runlen
    for b_row in range(b_rows):
        for b_col in range(b_cols):

            # precompute useful positions
            a_row = b_to_a0(b_row, b_col, inss, dels, r)
            a_col = b_to_a1(b_row, b_col, inss, dels, r)
            ref_idx = a_col - 1
            seq_idx = a_row - 1
            b_top0 = a_to_b0(a_row-1, a_col, inss, dels, r)
            b_top1 = a_to_b1(a_row-1, a_col, inss, dels, r)
            b_left0 = a_to_b0(a_row, a_col-1, inss, dels, r)
            b_left1 = a_to_b1(a_row, a_col-1, inss, dels, r)
            b_diag0 = a_to_b0(a_row-1, a_col-1, inss, dels, r)
            b_diag1 = a_to_b1(a_row-1, a_col-1, inss, dels, r)

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
                    matrix[typ, b_row, b_col, VALUE] = del_val
                    matrix[typ, b_row, b_col, TYPE] = DEL
                    matrix[typ, b_row, b_col, RUNLEN] = a_col
                continue
            elif a_col == 0:
                for typ in range(typs):
                    if typ == DEL or typ == SHP:
                        ins_val = indel_start + (100+indel_extend)*(a_row-1)
                    else:
                        ins_val = indel_start + indel_extend*(a_row-1)
                    matrix[typ, b_row, b_col, VALUE] = ins_val
                    matrix[typ, b_row, b_col, TYPE] = INS
                    matrix[typ, b_row, b_col, RUNLEN] = a_row
                continue

            # shouldn't be possible
            elif b_row == 0:
                print("ERROR: b_row == 0")

            # enforce new path remains within r cells of original path
            elif b_col == 0 or b_col == 2*r:
                max_val = get_max(matrix, b_row-1, b_col, typs, VALUE) + 100
                for typ in range(typs):
                    matrix[typ, b_row, b_col, VALUE] = max_val
                    matrix[typ, b_row, b_col, TYPE] = MAT
                    matrix[typ, b_row, b_col, RUNLEN] = 0
                continue


            # UPDATE INS MATRIX
            # continue INS
            if a_row == 1:
                ins_run = 1
            else:
                ins_run = <int>(matrix[INS, b_top0, b_top1, RUNLEN]) + 1
            ins_val = matrix[INS, b_top0, b_top1, VALUE] + indel_extend
            matrix[INS, b_row, b_col, VALUE] = ins_val
            matrix[INS, b_row, b_col, TYPE] = INS
            matrix[INS, b_row, b_col, RUNLEN] = ins_run

            # start INS
            start_val = matrix[MAT, b_top0, b_top1, VALUE] + indel_start
            if start_val < ins_val:
                matrix[INS, b_row, b_col, VALUE] = start_val
                matrix[INS, b_row, b_col, TYPE] = INS
                matrix[INS, b_row, b_col, RUNLEN] = 1


            # UPDATE DEL MATRIX
            # continue DEL
            if a_col == 1:
                del_run = 1
            else:
                del_run = <int>(matrix[DEL, b_left0, b_left1, RUNLEN]) + 1
            del_val = matrix[DEL, b_left0, b_left1, VALUE] + indel_extend
            matrix[DEL, b_row, b_col, VALUE] = del_val
            matrix[DEL, b_row, b_col, TYPE] = DEL
            matrix[DEL, b_row, b_col, RUNLEN] = del_run

            # start DEL
            start_val = matrix[MAT, b_left0, b_left1, VALUE] + indel_start
            if start_val < del_val:
                matrix[DEL, b_row, b_col, VALUE] = start_val
                matrix[DEL, b_row, b_col, TYPE] = DEL
                matrix[DEL, b_row, b_col, RUNLEN] = 1


            # UPDATE LHP MATRIX
            if seq_idx+2 < a_rows-1 and ref_idx+1 < a_cols-1 and \
                    seq[seq_idx+1] == seq[seq_idx+2]: # only insert same base

                # continue LHP
                if a_row == 1:
                    lhp_run = 1
                else:
                    lhp_run = <int>(matrix[LHP, b_top0, b_top1, RUNLEN]) + 1
                lhp_val = matrix[LHP, 
                        a_to_b0(a_row-lhp_run, a_col, inss, dels, r), 
                        a_to_b1(a_row-lhp_run, a_col, inss, dels, r), VALUE] + \
                        hp_indel_score(ref_hp_lens[ref_idx+1], lhp_run, hp_scores)
                matrix[LHP, b_row, b_col, VALUE] = lhp_val
                matrix[LHP, b_row, b_col, TYPE] = LHP
                matrix[LHP, b_row, b_col, RUNLEN] = lhp_run

                # start LHP
                start_val = matrix[MAT, b_top0, b_top1, VALUE] + \
                        hp_indel_score(ref_hp_lens[ref_idx+1], 1, hp_scores)
                if start_val < lhp_val:
                    matrix[LHP, b_row, b_col, VALUE] = start_val
                    matrix[LHP, b_row, b_col, TYPE] = LHP
                    matrix[LHP, b_row, b_col, RUNLEN] = 1

            else: # don't allow insertion of different base
                lhp_val = matrix[MAT, b_diag0, b_diag1, VALUE] + 100
                matrix[LHP, b_row, b_col, VALUE] = lhp_val
                matrix[LHP, b_row, b_col, TYPE] = LHP
                matrix[LHP, b_row, b_col, RUNLEN] = 0


            # UPDATE SHP MATRIX
            if ref_idx+1 < a_cols-1 and \
                    ref[ref_idx+1] == ref[ref_idx]: # only delete same base

                # continue SHP
                if a_col == 1:
                    shp_run = 1
                else:
                    shp_run = <int>(matrix[SHP, b_left0, b_left1, RUNLEN] + 1)
                shp_val = matrix[SHP, 
                        a_to_b0(a_row, a_col-shp_run, inss, dels, r), 
                        a_to_b1(a_row, a_col-shp_run, inss, dels, r), VALUE] + \
                        hp_indel_score(ref_hp_lens[ref_idx], -shp_run, hp_scores)
                matrix[SHP, b_row, b_col, VALUE] = shp_val
                matrix[SHP, b_row, b_col, TYPE] = SHP
                matrix[SHP, b_row, b_col, RUNLEN] = shp_run

                # start SHP
                start_val = matrix[MAT, b_left0, b_left1, VALUE] + \
                    hp_indel_score(ref_hp_lens[ref_idx], -1, hp_scores)
                if start_val < shp_val:
                    matrix[SHP, b_row, b_col, VALUE] = start_val
                    matrix[SHP, b_row, b_col, TYPE] = SHP
                    matrix[SHP, b_row, b_col, RUNLEN] = 1

            else: # don't allow deleting different base
                shp_val = matrix[MAT, b_diag0, b_diag1, VALUE] + 100
                matrix[SHP, b_row, b_col, VALUE] = shp_val
                matrix[SHP, b_row, b_col, TYPE] = SHP
                matrix[SHP, b_row, b_col, RUNLEN] = 0


            # UPDATE MAT MATRIX
            if matrix[MAT, b_diag0, b_diag1, TYPE] == MAT:
                runlen = <int>(matrix[MAT, b_diag0, b_diag1, RUNLEN]) + 1
            else:
                runlen = 1
            sub_val = matrix[MAT, b_diag0, b_diag1, VALUE] + \
                    sub_scores[ seq[seq_idx], ref[ref_idx] ]
            min_val = sub_val
            matrix[MAT, b_row, b_col, VALUE] = min_val
            matrix[MAT, b_row, b_col, TYPE] = MAT
            matrix[MAT, b_row, b_col, RUNLEN] = runlen

            # end INDEL
            if verbosity >= 2:
                print("REF:", ref)
                s = "     "
                for x in range(ref_idx):
                    s += " "
                s += "^"
                print(s)

                print("SEQ:", seq)
                s = "     "
                for x in range(seq_idx):
                    s += " "
                s += "^"
                print(s)
                print("M", sub_val)
            for typ in [INS, LHP, DEL, SHP]:

                end_val = matrix[typ, b_row, b_col, VALUE]
                if verbosity >= 2:
                    chars = "MILDS"
                    print(chars[typ], end_val)
                if end_val < min_val:
                    min_val = end_val
                    runlen = <int>(matrix[typ, b_row, b_col, RUNLEN])
                    matrix[MAT, b_row, b_col, VALUE] = min_val
                    matrix[MAT, b_row, b_col, TYPE] = typ
                    matrix[MAT, b_row, b_col, RUNLEN] = runlen
            if verbosity >= 2:
                print("\n")


    # initialize backtracking from last cell
    cdef str aln = ''
    cdef str op
    path = []
    cdef int row = a_rows - 1
    cdef int col = a_cols - 1
    cdef float val
    b_pos0 = a_to_b0(row, col, inss, dels, r)
    b_pos1 = a_to_b1(row, col, inss, dels, r)
    typ = <int>(matrix[MAT, b_pos0, b_pos1, VALUE])

    # backtrack
    while row > 0 or col > 0:
        if row < 0:
            print("\nERROR: row < 0 @  A: (", row, 
                    ",", col, "), B: (", b_pos0, ",", b_pos1, ").")
            break
        if col < 0:
            print("\nERROR: col < 0 @  A: (", row, 
                    ",", col, "), B: (", b_pos0, ",", b_pos1, ").")
            break
        path.append((MAT, row, col))
        b_pos0 = a_to_b0(row, col, inss, dels, r)
        b_pos1 = a_to_b1(row, col, inss, dels, r)
        val = matrix[MAT, b_pos0, b_pos1, VALUE]
        typ = <int>(matrix[MAT, b_pos0, b_pos1, TYPE])
        runlen = <int>(matrix[MAT, b_pos0, b_pos1, RUNLEN])

        # Invalid SHPs and LHPs are runlen 0. They should never be reached during 
        # backtracking (due to high penalties), but leaving this here just-in-case
        if runlen < 1: runlen = 1
        op = ''
        if typ == LHP or typ == INS:   # each move is an insertion
            for i in range(runlen):
                op += 'I'
            row -= runlen
        elif typ == SHP or typ == DEL: # each move is a deletion
            for i in range(runlen):
                op += 'D'
            col -= runlen
        elif typ == MAT: # only sub if stay in same matrix
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
    if verbosity >= 1:
        types = ['MAT', 'INS', 'LHP', 'DEL', 'SHP']
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
                    b_pos0 = a_to_b0(row, col, inss, dels, r)
                    b_pos1 = a_to_b1(row, col, inss, dels, r)
                    op = ops[<int>(matrix[typ, b_pos0, b_pos1, TYPE])]
                    if (typ, row, col) in path:
                        mark = '$'
                    elif col == 0 or row == 0:
                        mark = '*'
                    elif b_pos1 == 0 or b_pos1 == 2*r:
                        mark = '.'
                    else:
                        mark = ' '
                    runlen = <int>(matrix[typ, b_pos0, b_pos1, RUNLEN])
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
                    a_row = b_to_a0(row, col, inss, dels, r)
                    a_col = b_to_a1(row, col, inss, dels, r)
                    op = ops[<int>(matrix[typ, row, col, TYPE])]
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
                    runlen = <int>(matrix[typ, row, col, RUNLEN])
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

