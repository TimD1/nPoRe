import numpy as np
import math
np.set_printoptions(linewidth=200)
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import scipy.ndimage as ndimage

import cython

import cfg

from cig import *


def fix_matrix_properties(scores, delta = 0.01):
    ''' Modify score matrix to adhere to the following properties:
        - diagonals should all have same low penalty (correct call)
        - a longer INDEL from same N-polymer should be penalized more
        - a similar INDEL from a longer N-polymer should be penalized less
    '''

    ns = scores.shape[0]
    l = scores.shape[1]
    INF = 10000

    for n in range(ns):

        # don't penalize diagonals
        for i in range(1, l):
            scores[n,i,i] = min(scores[n,i,i], scores[n,i-1,i-1])

        # more insertions should be more penalized
        for j in range(1, l):
            for i in range(j-1, -1, -1):
                scores[n,i,j] = max(
                        scores[n,i,j], 
                        scores[n,i+1,j] + delta, 
                        scores[n,i,j-1] + delta
                )

        # more deletions should be more penalized
        for i in range(1,l):
            for j in range(i-1, -1, -1):
                scores[n,i,j] = max(
                        scores[n,i,j], 
                        scores[n,i,j+1] + delta, 
                        scores[n,i-1,j] + delta
                )

        # prefer insertions from longer homopolymers
        best = np.ones(l) * INF
        for j in range(1,l):
            for i in range(j-1, -1, -1):
                ins_len = j - i
                if scores[n,i,j] < best[ins_len]:
                    best[ins_len] = scores[n,i,j]
                    for total_ins_len in range(ins_len+1, l):
                        best[total_ins_len] = min(
                                best[total_ins_len], 
                                best[ins_len] + best[total_ins_len-ins_len]
                        )
                else:
                    scores[n,i,j] = min(
                            scores[n,i,j], 
                            best[ins_len] - delta
                    )

        # prefer deletions from longer homopolymers
        best = np.ones(l) * INF
        for i in range(1,l):
            for j in range(i-1, -1, -1):
                del_len = i - j
                if scores[n,i,j] < best[del_len]:
                    best[del_len] = scores[n,i,j]
                    for total_del_len in range(del_len+1, l):
                        best[total_del_len] = min(
                                best[total_del_len], 
                                best[del_len] + best[total_del_len-del_len]
                        )
                else:
                    scores[n,i,j] = min(
                            scores[n,i,j], 
                            best[del_len] - delta
                    )

    return scores



def calc_score_matrices(subs, nps, inss, dels):

    # calculate homopolymer scores matrix
    np_scores = np.zeros_like(nps, dtype=np.float32)
    for n in range(cfg.args.max_np):
        for ref_len in range(cfg.args.max_np_len):
            total = np.sum(nps[n, ref_len])
            for call_len in range(cfg.args.max_np_len):
                count = int(nps[n, ref_len, call_len])
                bias = 10
                frac = (count + 0.01 + int(ref_len==call_len)*bias) / (total + 0.01*cfg.args.max_np_len + bias)
                np_scores[n, ref_len, call_len] = -np.log(frac)
    np_scores = fix_matrix_properties(np_scores)

    # calculate substitution scores matrix
    sub_scores = np.zeros((cfg.nbases,cfg.nbases), dtype=np.float32)
    for i in range(1, cfg.nbases):
        for j in range(1, cfg.nbases):
            if i != j:
                sub_scores[i, j] = -np.log( (subs[i,j]+0.01) / (np.sum(subs[i])+0.1) )
            else:
                sub_scores[i, j] = 0

    ins_scores = np.zeros_like(inss, dtype=np.float32)
    total = np.sum(inss)
    for l in range(cfg.args.max_np_len):
        frac = (inss[l] + 0.01) / (total + 0.01*cfg.args.max_np_len)
        ins_scores[l] = -np.log(frac)

    del_scores = np.zeros_like(dels, dtype=np.float32)
    total = np.sum(dels)
    for l in range(cfg.args.max_np_len):
        frac = (dels[l] + 0.01) / (total + 0.01*cfg.args.max_np_len)
        del_scores[l] = -np.log(frac)

    return sub_scores, np_scores, ins_scores, del_scores



def plot_np_score_matrices(nps, max_np_len = 40):
    for n in range(cfg.args.max_np):
        x, y = np.meshgrid(range(max_np_len), range(max_np_len))
        fig = plt.figure(figsize=(20,10))
        ax1 = fig.add_subplot(1,2,1, projection='3d')
        ax2 = fig.add_subplot(1,2,2, projection='3d')
        ax1.plot_trisurf(x.flatten(), -y.flatten(), 
                nps[n,:max_np_len,:max_np_len].flatten(), 
                cmap='RdYlGn_r', edgecolor='none')
        ax2.plot_trisurf(-x.flatten(), -y.flatten(), 
                nps[n,:max_np_len,:max_np_len].flatten(), 
                cmap='RdYlGn_r', edgecolor='none')
        plt.tight_layout()
        plt.savefig(f'{cfg.args.stats_dir}/{n+1}-polymer_scores.png', dpi=200)
        plt.close()



@cython.boundscheck(False)
@cython.wraparound(False)
cpdef int[:] get_ns(char[:] seq):
    ''' Calculate N-polymer type and length of substring starting at each index. 
         seq:     A T A T A T T T T T T T A A A
         ns:      2 2 2 0 0 1 1 1 1 1 1 0 1 1 0
         np_lens: 3 2 2 0 0 7 6 5 4 3 2 0 3 2 0
    '''

    cdef int seq_len = len(seq)
    np_lens_buf = np.zeros(seq_len, dtype=np.intc)
    ns_buf = np.zeros(seq_len, dtype=np.intc)
    cdef int[:] np_lens = np_lens_buf
    cdef int[:] ns = ns_buf
    cdef int n, n_idx, np_repeat_len
    cdef int seq_idx, seq_ptr

    for seq_idx in range(seq_len): # iterate over sequence

        best_len = 0
        for n in range(1, cfg.args.max_np+1): # check each length N-polymer

            # get np repeat length at this position
            np_repeat_len = 0
            seq_ptr = seq_idx
            while seq_ptr+n < seq_len and seq[seq_ptr] == seq[seq_ptr+n]:
                seq_ptr += 1
                if (seq_ptr-seq_idx) % n == 0: # finished n-polymer
                    np_repeat_len += 1
            if np_repeat_len: np_repeat_len += 1 # count first

            # save best n and n-polymer repeat lengths
            if n * np_repeat_len > ns[seq_idx] * np_lens[seq_idx]:
                ns[seq_idx] = n
                np_lens[seq_idx] = np_repeat_len

    return ns



@cython.boundscheck(False)
@cython.wraparound(False)
cpdef int[:] get_np_lens(char[:] seq):
    ''' Calculate N-polymer type and length of substring starting at each index. 
         seq:     A T A T A T T T T T T T A A A
         ns:      2 2 2 0 0 1 1 1 1 1 1 0 1 1 0
         np_lens: 3 2 2 0 0 7 6 5 4 3 2 0 3 2 0
    '''

    cdef int seq_len = len(seq)
    np_lens_buf = np.zeros(seq_len, dtype=np.intc)
    ns_buf = np.zeros(seq_len, dtype=np.intc)
    cdef int[:] np_lens = np_lens_buf
    cdef int[:] ns = ns_buf
    cdef int n, n_idx, np_repeat_len
    cdef int seq_idx, seq_ptr

    for seq_idx in range(seq_len): # iterate over sequence

        best_len = 0
        for n in range(1, cfg.args.max_np+1): # check each length N-polymer

            # get np repeat length at this position
            np_repeat_len = 0
            seq_ptr = seq_idx
            while seq_ptr+n < seq_len and seq[seq_ptr] == seq[seq_ptr+n]:
                seq_ptr += 1
                if (seq_ptr-seq_idx) % n == 0: # finished n-polymer
                    np_repeat_len += 1
            if np_repeat_len: np_repeat_len += 1 # count first

            # save best n and n-polymer repeat lengths
            if n * np_repeat_len > ns[seq_idx] * np_lens[seq_idx]:
                ns[seq_idx] = n
                np_lens[seq_idx] = np_repeat_len

    return np_lens



@cython.boundscheck(False)
@cython.wraparound(False)
cdef float np_score(int ref_np_len, int indel_len, float[:,:] np_scores):

    # error, don't allow
    if ref_np_len <= 0:
        return 100
    elif ref_np_len + indel_len < 0:
        return 100

    # force lengths to fit in matrix
    cdef int np_scores_len = len(np_scores)
    if ref_np_len > np_scores_len-1:
        ref_np_len = np_scores_len-1

    cdef int call_hp_len = ref_np_len + indel_len
    if call_hp_len > np_scores_len-1:
        call_hp_len = np_scores_len-1
    return np_scores[ref_np_len, call_hp_len]


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
cdef int a_to_b_row(int a_row, int a_col, int[:] inss, int[:] dels, int r):
    return a_row + a_col

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int a_to_b_col(int a_row, int a_col, int[:] inss, int[:] dels, int r):
    cdef int b_row, b_col
    b_row = a_row + a_col
    b_col = inss[b_row] - a_row + r
    return b_col



@cython.boundscheck(False)
@cython.wraparound(False)
cdef int b_to_a_row(int b_row, int b_col, int[:] inss, int[:] dels, int r):
    return inss[b_row] + r - b_col

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int b_to_a_col(int b_row, int b_col, int[:] inss, int[:] dels, int r):
    return dels[b_row] - r + b_col



@cython.boundscheck(False)
@cython.wraparound(False)
cdef float get_max(float[:,:,:,:] matrix, int row, int col, int typs, int VAL):
    cdef int typ
    cdef float max_val = matrix[0, row, col, VAL]
    for typ in range(1, typs):
        if matrix[typ, row, col, VAL] > max_val:
            max_val = matrix[typ, row, col, VAL]
    return max_val



@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[:] get_breaks(int chunk_size, int array_size, int[:] inss, int[:] dels):
    cdef int buf_len = 1 + math.ceil( (array_size-1) / (chunk_size-1) )
    breaks_buf = np.zeros(buf_len, dtype=np.intc)
    cdef int[:] breaks = breaks_buf
    cdef int i
    for i in range(buf_len-1):
        breaks[i] = i * (chunk_size-1)

        # don't split on 'DI', since that may have originally been '='
        if i > 0 and inss[breaks[i]+1] == inss[breaks[i]]+1 and \
                dels[breaks[i]] == dels[breaks[i]-1]+1:
            breaks[i] -= 1

    breaks[buf_len-1] = array_size-1
    return breaks



@cython.boundscheck(False)
@cython.wraparound(False)
cpdef align(char[::1] ref, char[::1] seq, str cigar, 
        float[:,::1] sub_scores, float[:,::1] np_scores, 
        float indel_start=5, float indel_extend=2, int max_b_rows = 20001,
        int r = 30):
    ''' Perform alignment.  '''

    # convert CIGAR so that each movement is row+1 or col+1, enables easy banding
    cigar = cigar.replace('X','DI').replace('=','DI') \
            .replace('M','DI').replace('S','').replace('H','')

    # precompute offsets, breakpoints, and homopolymers
    cdef int[:] inss = get_inss(cigar)
    cdef int[:] dels = get_dels(cigar)
    cdef int[:] breaks = get_breaks(max_b_rows, len(seq) + len(ref) + 1, inss, dels)
    cdef int[:] ref_nps = get_ns(ref)
    cdef int[:] ref_np_lens = get_np_lens(ref)

    # define useful constants
    cdef int a_rows = len(seq) + 1
    cdef int a_cols = len(ref) + 1
    cdef int b_cols = 2*r + 1
    cdef int b_rows = -1

    cdef int dims = 3 # dimensions in matrix
    cdef int VAL = 0  # value (alignment score)
    cdef int TYP = 1  # type (predecessor matrix dimension)
    cdef int RUN = 2  # run length (for indels)

    cdef int typs = 5 # types
    cdef int MAT = 0  # match/substitution
    cdef int INS = 1  # insertion
    cdef int NPI = 2  # lengthen homopolymer
    cdef int DEL = 3  # deletion
    cdef int NPD = 4  # shorten homopolymer

    cdef str aln = ''
    cdef str full_aln = ''
    cdef str op
    matrix_buf = np.zeros((typs, max_b_rows, b_cols, dims), dtype=np.float32)
    cdef float[:,:,:,::1] matrix = matrix_buf

    cdef int b_row, b_col, a_row, a_col, ref_idx, seq_idx
    cdef int b_top_row, b_top_col, b_left_row, b_left_col, b_diag_row, b_diag_col
    cdef int run, typ, i, brk, next_brk, brk_idx
    cdef float val1, val2

    # iterate over b matrix in chunks set by breakpoints
    for brk_idx in range(len(breaks)-1):
        brk = breaks[brk_idx]
        next_brk = breaks[brk_idx+1]
        b_rows = next_brk - brk + 1
        matrix_buf.fill(0)

        # calculate matrix chunk
        for b_row in range(b_rows):
            for b_col in range(b_cols):

                # precompute useful positions
                a_row = b_to_a_row(b_row + brk, b_col, inss, dels, r)
                a_col = b_to_a_col(b_row + brk, b_col, inss, dels, r)
                b_top_row = a_to_b_row(a_row-1, a_col, inss, dels, r) - brk
                b_top_col = a_to_b_col(a_row-1, a_col, inss, dels, r)
                b_left_row = a_to_b_row(a_row, a_col-1, inss, dels, r) - brk
                b_left_col = a_to_b_col(a_row, a_col-1, inss, dels, r)
                b_diag_row = a_to_b_row(a_row-1, a_col-1, inss, dels, r) - brk
                b_diag_col = a_to_b_col(a_row-1, a_col-1, inss, dels, r)
                ref_idx = a_col - 1
                seq_idx = a_row - 1

                # skip cells out of range of this chunk of original "A" matrix
                if a_row < inss[brk] or a_col < dels[brk] or \
                        a_row > inss[next_brk] or a_col > dels[next_brk]:
                    continue

                # very first cell has no score (necessary to skip other elifs)
                elif a_row == inss[brk] and a_col == dels[brk]:
                    continue

                # initialize first row/col of matrix A
                elif a_row == inss[brk]:
                    for typ in range(typs):
                        if typ == INS or typ == NPI: # don't allow
                            val1 = indel_start + (100+indel_extend)*(a_col-dels[brk]-1)
                        else:
                            val1 = indel_start + indel_extend*(a_col-dels[brk]-1)
                        matrix[typ, b_row, b_col, VAL] = val1
                        matrix[typ, b_row, b_col, TYP] = DEL
                        matrix[typ, b_row, b_col, RUN] = a_col - dels[brk]
                    continue
                elif a_col == dels[brk]:
                    for typ in range(typs):
                        if typ == DEL or typ == NPD: # don't allow
                            val1 = indel_start + (100+indel_extend)*(a_row-inss[brk]-1)
                        else:
                            val1 = indel_start + indel_extend*(a_row-inss[brk]-1)
                        matrix[typ, b_row, b_col, VAL] = val1
                        matrix[typ, b_row, b_col, TYP] = INS
                        matrix[typ, b_row, b_col, RUN] = a_row - inss[brk]
                    continue

                # enforce new path remains within r cells of original path
                elif b_col == 0 or b_col == 2*r:
                    val1 = get_max(matrix, b_row-1, b_col, typs, VAL) + 100
                    for typ in range(typs):
                        matrix[typ, b_row, b_col, VAL] = val1
                        matrix[typ, b_row, b_col, TYP] = MAT
                        matrix[typ, b_row, b_col, RUN] = 0
                    continue


                # UPDATE INS MATRIX
                val1 = matrix[MAT, b_top_row, b_top_col, VAL] + indel_start
                val2 = matrix[INS, b_top_row, b_top_col, VAL] + indel_extend
                if val1 < val2: # start insertion
                    matrix[INS, b_row, b_col, VAL] = val1
                    matrix[INS, b_row, b_col, TYP] = INS
                    matrix[INS, b_row, b_col, RUN] = 1
                else: # continue insertion
                    if a_row == inss[brk] + 1:
                        run = 1
                    else:
                        run = <int>(matrix[INS, b_top_row, b_top_col, RUN]) + 1
                    matrix[INS, b_row, b_col, VAL] = val2
                    matrix[INS, b_row, b_col, TYP] = INS
                    matrix[INS, b_row, b_col, RUN] = run


                # UPDATE DEL MATRIX
                val1 = matrix[MAT, b_left_row, b_left_col, VAL] + indel_start
                val2 = matrix[DEL, b_left_row, b_left_col, VAL] + indel_extend
                if val1 < val2: # start deletion
                    matrix[DEL, b_row, b_col, VAL] = val1
                    matrix[DEL, b_row, b_col, TYP] = DEL
                    matrix[DEL, b_row, b_col, RUN] = 1
                else: # continue deletion
                    if a_col == dels[brk] + 1:
                        run = 1
                    else:
                        run = <int>(matrix[DEL, b_left_row, b_left_col, RUN]) + 1
                    matrix[DEL, b_row, b_col, VAL] = val2
                    matrix[DEL, b_row, b_col, TYP] = DEL
                    matrix[DEL, b_row, b_col, RUN] = run


                # UPDATE NPI MATRIX
                np_start = ref_idx == 0 and ref_nps[0] != 0 or \
                        ref_nps[ref_idx] != ref_nps[ref_idx-1] or \
                        ref_np_lens[ref_idx] > ref_np_lens[ref_idx-1]

                # # UPDATE NPI MATRIX
                # if seq_idx+2 < a_rows-1 and ref_idx+1 < a_cols-1 and \
                #         seq[seq_idx+1] == seq[seq_idx+2]: # only insert same base

                #     val1 = matrix[MAT, b_top_row, b_top_col, VAL] + \
                #             np_score(ref_np_lens[ref_idx+1], 1, np_scores)
                #     if a_row == inss[brk] + 1:
                #         run = 1
                #     else:
                #         run = <int>(matrix[NPI, b_top_row, b_top_col, RUN]) + 1
                #     val2 = matrix[NPI, 
                #             a_to_b_row(a_row-run, a_col, inss, dels, r) - brk, 
                #             a_to_b_col(a_row-run, a_col, inss, dels, r), VAL] + \
                #             np_score(ref_np_lens[ref_idx+1], run, np_scores)
                #     if val1 < val2: # start lengthen
                #         matrix[NPI, b_row, b_col, VAL] = val1
                #         matrix[NPI, b_row, b_col, TYP] = NPI
                #         matrix[NPI, b_row, b_col, RUN] = 1
                #     else: # continue lengthen
                #         matrix[NPI, b_row, b_col, VAL] = val2
                #         matrix[NPI, b_row, b_col, TYP] = NPI
                #         matrix[NPI, b_row, b_col, RUN] = run

                # else: # don't allow insertion of different base
                #     val2 = matrix[MAT, b_diag_row, b_diag_col, VAL] + 100
                #     matrix[NPI, b_row, b_col, VAL] = val2
                #     matrix[NPI, b_row, b_col, TYP] = NPI
                #     matrix[NPI, b_row, b_col, RUN] = 0


                # UPDATE NPD MATRIX
                if ref_idx+1 < a_cols-1 and \
                        ref[ref_idx+1] == ref[ref_idx]: # only delete same base

                    val1 = matrix[MAT, b_left_row, b_left_col, VAL] + \
                        np_score(ref_np_lens[ref_idx], -1, np_scores)
                    if a_col == dels[brk] + 1:
                        run = 1
                    else:
                        run = <int>(matrix[NPD, b_left_row, b_left_col, RUN] + 1)
                    val2 = matrix[NPD, 
                            a_to_b_row(a_row, a_col-run, inss, dels, r) - brk, 
                            a_to_b_col(a_row, a_col-run, inss, dels, r), VAL] + \
                            np_score(ref_np_lens[ref_idx], -run, np_scores)
                    if val1 < val2: # start shorten
                        matrix[NPD, b_row, b_col, VAL] = val1
                        matrix[NPD, b_row, b_col, TYP] = NPD
                        matrix[NPD, b_row, b_col, RUN] = 1
                    else: # continue shorten
                        matrix[NPD, b_row, b_col, VAL] = val2
                        matrix[NPD, b_row, b_col, TYP] = NPD
                        matrix[NPD, b_row, b_col, RUN] = run

                else: # don't allow deleting different base
                    val2 = matrix[MAT, b_diag_row, b_diag_col, VAL] + 100
                    matrix[NPD, b_row, b_col, VAL] = val2
                    matrix[NPD, b_row, b_col, TYP] = NPD
                    matrix[NPD, b_row, b_col, RUN] = 0


                # UPDATE MAT MATRIX
                if matrix[MAT, b_diag_row, b_diag_col, TYP] == MAT:
                    run = <int>(matrix[MAT, b_diag_row, b_diag_col, RUN]) + 1
                else:
                    run = 1
                val1 = matrix[MAT, b_diag_row, b_diag_col, VAL] + \
                        sub_scores[ seq[seq_idx], ref[ref_idx] ]
                matrix[MAT, b_row, b_col, VAL] = val1
                matrix[MAT, b_row, b_col, TYP] = MAT
                matrix[MAT, b_row, b_col, RUN] = run

                # end INDEL
                for typ in [INS, NPI, DEL, NPD]:
                    val2 = matrix[typ, b_row, b_col, VAL]
                    if val2 < val1:
                        val1 = val2
                        run = <int>(matrix[typ, b_row, b_col, RUN])
                        matrix[MAT, b_row, b_col, VAL] = val2
                        matrix[MAT, b_row, b_col, TYP] = typ
                        matrix[MAT, b_row, b_col, RUN] = run


        # initialize backtracking from last cell
        a_row = inss[next_brk]
        a_col = dels[next_brk]
        aln = ''
        b_row = a_to_b_row(a_row, a_col, inss, dels, r) - brk
        b_col = a_to_b_col(a_row, a_col, inss, dels, r)
        typ = <int>(matrix[MAT, b_row, b_col, VAL])
        run = 0

        # backtrack
        while a_row > inss[brk] or a_col > dels[brk]:
            b_row = a_to_b_row(a_row, a_col, inss, dels, r) - brk
            b_col = a_to_b_col(a_row, a_col, inss, dels, r)
            val = matrix[MAT, b_row, b_col, VAL]
            typ = <int>(matrix[MAT, b_row, b_col, TYP])
            run = <int>(matrix[MAT, b_row, b_col, RUN])

            # Invalid NPDs and NPIs are run 0. They should never be reached during 
            # backtracking (due to high penalties), but leaving this here just-in-case
            if run < 1: run = 1

            op = ''
            if typ == NPI or typ == INS:   # each move is an insertion
                for i in range(run):
                    op += 'I'
                a_row -= run
            elif typ == NPD or typ == DEL: # each move is a deletion
                for i in range(run):
                    op += 'D'
                a_col -= run
            elif typ == MAT: # only sub if stay in same matrix
                i = 0
                while i < run:
                    a_row -= 1
                    a_col -= 1
                    if ref[a_col] == seq[a_row]:
                        op += '='
                    else:
                        op += 'X'
                    i += 1
            else:
                print("ERROR: unknown alignment matrix type:", typ)
                break
            aln += op

        full_aln += aln[::-1]

    return full_aln



def dump(ref, seq, cigar):
    ''' Pretty print full alignment result. '''

    ref_str = ''
    cig_str = ''
    seq_str = ''

    ref_idx = 0
    seq_idx = 0

    for idx, op in enumerate(cigar):
        if op == '=' or op == 'M':
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

