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



def plot_np_score_matrices(nps, max_np_len = 20):
    for n in range(cfg.args.max_np):

        # score matrix
        plt.figure(figsize=(max_np_len,max_np_len))
        plt.matshow(nps[n,:max_np_len, :max_np_len], cmap='RdYlGn_r',)
        for i in range(max_np_len):
            for j in range(max_np_len):
                plt.text(x=j, y=i, s=f'{nps[n,i,j]:.1f}', fontsize=5, 
                        va='center', ha='center')
        plt.xlabel('Called')
        plt.ylabel('Actual')
        plt.xticks(range(max_np_len))
        plt.yticks(range(max_np_len))
        plt.title(f'{n+1}-Polymer Score Matrix')
        plt.savefig(f'{cfg.args.stats_dir}/{n+1}-polymer_scores.png', dpi=300)
        plt.close()

        # surface plot
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
        plt.savefig(f'{cfg.args.stats_dir}/{n+1}-polymer_scores_surface.png', dpi=200)
        plt.close()



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int[:,::1] get_np_info(char[::1] seq):
    ''' Calculate N-polymer information. 

         seq:     A T A T A T T T T T T T A A A G C

         np_info:
         RPTS:    3 3 3 3 3 7 7 7 7 7 7 7 3 3 3 0 0
         RPT:     0 0 1 1 2 0 1 2 3 4 5 6 0 1 2 0 0
         N:       2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 0 0
         IDX:     0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
         
         RPTS = number of times the sequence is repeated
         RPT = 0-based index of the current repeat
         N = number of bases in repeated sequence
         IDX = 0-based index of current base within repeat

         Explanation:
         3(AT), 7(T), 3(A) with some overlap. A sequence must repeat at least 
         twice to be considered an n-polymer. Bases are considered part of 
         the longest repeat in which they're included (N * RPTS).

    '''

    cdef int seq_len = len(seq)
    np_info_buf = np.zeros((4, seq_len), dtype=np.intc)
    cdef int[:,::1] np_info = np_info_buf
    cdef int n, np_repeat_len, pos, rpt, idx
    cdef int seq_idx, seq_ptr

    # define constant values for indexing into `np_info` array
    cdef int RPTS = 0
    cdef int RPT = 1
    cdef int N = 2
    cdef int IDX = 3

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

            # save n-polymer info
            if np_repeat_len > 2 and \
                    n * np_repeat_len > np_info[N, seq_idx] * np_info[RPTS, seq_idx]:
                for rpt in range(np_repeat_len):
                    for idx in range(n):
                        pos = seq_idx + rpt*n + idx
                        np_info[RPTS, pos] = np_repeat_len
                        np_info[RPT, pos] = rpt
                        np_info[N, pos] = n
                        np_info[IDX, pos] = idx

    return np_info



@cython.boundscheck(False)
@cython.wraparound(False)
cdef float np_score(int n, int ref_np_len, int indel_len, float[:,:,::1] np_scores, int max_np):

    # error, don't allow
    if ref_np_len <= 0:
        return 100
    elif ref_np_len + indel_len < 0:
        return 100
    elif n < 1 or n > max_np:
        return 100

    # force lengths to fit in matrix
    cdef int call_np_len = ref_np_len + indel_len
    if ref_np_len > max_np-1:
        ref_np_len = max_np-1
    if call_np_len > max_np-1:
        call_np_len = max_np-1

    return np_scores[n-1, ref_np_len, call_np_len]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[::1] get_inss(str cigar):
    ''' CIGAR must contain only "I" and "D". '''

    cdef int cig_len = len(cigar)
    inss_buf = np.zeros(cig_len+1, dtype=np.intc)
    cdef int[::1] inss = inss_buf
    cdef int i

    for i in range(cig_len):
        if cigar[i] == 'I':
            inss[i+1] = inss[i] + 1
        else:
            inss[i+1] = inss[i]
    return inss



@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[::1] get_dels(str cigar):
    ''' CIGAR must contain only "I" and "D". '''

    cdef int cig_len = len(cigar)
    dels_buf = np.zeros(cig_len+1, dtype=np.intc)
    cdef int[::1] dels = dels_buf
    cdef int i

    for i in range(cig_len):
        if cigar[i] == 'D':
            dels[i+1] = dels[i] + 1
        else:
            dels[i+1] = dels[i]
    return dels



@cython.boundscheck(False)
@cython.wraparound(False)
cdef int a_to_b_row(int a_row, int a_col, int[::1] inss, int[::1] dels, int r):
    return a_row + a_col

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int a_to_b_col(int a_row, int a_col, int[::1] inss, int[::1] dels, int r):
    cdef int b_row, b_col
    b_row = a_row + a_col
    b_col = inss[b_row] - a_row + r
    return b_col



@cython.boundscheck(False)
@cython.wraparound(False)
cdef int b_to_a_row(int b_row, int b_col, int[::1] inss, int[::1] dels, int r):
    return inss[b_row] + r - b_col

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int b_to_a_col(int b_row, int b_col, int[::1] inss, int[::1] dels, int r):
    return dels[b_row] - r + b_col



@cython.boundscheck(False)
@cython.wraparound(False)
cdef int[::1] get_breaks(int chunk_size, int array_size, int[::1] inss, int[::1] dels):
    cdef int buf_len = 1 + math.ceil( (array_size-1) / (chunk_size-1) )
    breaks_buf = np.zeros(buf_len, dtype=np.intc)
    cdef int[::1] breaks = breaks_buf
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
cdef int match(char[::1] A, char[::1] B):
    cdef int i
    if len(A) != len(B): 
        return 0
    else:
        for i in range(len(A)):
            if A[i] != B[i]:
                return 0
    return 1



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef align(char[::1] full_ref, char[::1] seq, str cigar, 
        float[:,::1] sub_scores, float[:,:,::1] np_scores, 
        float indel_start=5, float indel_extend=2, int max_b_rows = 20000,
        int r = 30, int verbose=0):
    ''' Perform alignment.  '''

    # convert CIGAR so that each movement is row+1 or col+1, enables easy banding
    cigar = cigar.replace('X','DI').replace('=','DI').replace('M','DI')

    # precompute offsets, breakpoints, and homopolymers
    cdef int[::1] inss = get_inss(cigar)
    cdef int[::1] dels = get_dels(cigar)
    cdef int[::1] breaks = get_breaks(max_b_rows, len(seq) + len(full_ref) + 1, inss, dels)
    cdef int[:,::1] np_info
    cdef char[::1] ref

    # define useful constants
    cdef int a_rows = len(seq) + 1
    cdef int a_cols = len(full_ref) + 1
    cdef int b_cols = 2*r + 1
    cdef int b_rows = -1

    # n-polymer info indices
    cdef int RPTS = 0
    cdef int RPT = 1
    cdef int N = 2
    cdef int IDX = 3

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
    cdef str op, s
    # for alignment offset purposes, M -> DI. Breaks are sometimes shifted by 1, 
    # so that this DI doesn't cross the boundary of alignment chunks. 
    # As a result, we must increase the matrix buffer size. (hence max_b_rows+1)
    matrix_buf = np.zeros((typs, max_b_rows+1, b_cols, dims), dtype=np.float32)
    cdef float[:,:,:,::1] matrix = matrix_buf
    # max single-move penalty is < 10, so a path with this penalty will never 
    # be chosen. Still kept small enough so that len(seq)*INF < INT_MAX
    cdef int INF = 100

    cdef int b_row, b_col, a_row, a_col, ref_idx, seq_idx, row, col
    cdef int b_top_row, b_top_col, b_left_row, b_left_col, b_diag_row, b_diag_col
    cdef int b_runleft_row, b_runleft_col, b_runup_row, b_runup_col
    cdef int b_ndown_row, b_ndown_col, b_nright_row, b_nright_col
    cdef int run, typ, i, brk, next_brk, brk_idx
    cdef int rpts, rpt, n, idx
    cdef int max_np_len = cfg.args.max_np_len
    cdef float val1, val2

    # iterate over b matrix in chunks set by breakpoints
    for brk_idx in range(len(breaks)-1):

        brk = breaks[brk_idx]
        next_brk = breaks[brk_idx+1]
        b_rows = next_brk - brk + 1
        matrix_buf.fill(0)
        ref = full_ref[ dels[brk] : dels[next_brk]+1 ]
        np_info = get_np_info(ref)

        # initialize N-polymer matrices with invalid states
        for b_row in range(b_rows):
            for b_col in range(b_cols):
                a_row = b_to_a_row(b_row + brk, b_col, inss, dels, r)
                a_col = b_to_a_col(b_row + brk, b_col, inss, dels, r)
                if a_row < inss[brk] or a_col < dels[brk] or \
                        a_row > inss[next_brk] or a_col > dels[next_brk] or \
                        b_col == 0 or b_col == 2*r:
                    continue
                matrix[NPI, b_row, b_col, VAL] = INF * (a_row-inss[brk] + a_col-dels[brk])
                matrix[NPI, b_row, b_col, TYP] = MAT
                matrix[NPI, b_row, b_col, RUN] = 0
                matrix[NPD, b_row, b_col, VAL] = INF * (a_row-inss[brk] + a_col-dels[brk])
                matrix[NPD, b_row, b_col, TYP] = MAT
                matrix[NPD, b_row, b_col, RUN] = 0

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
                ref_idx = a_col - dels[brk] - 1
                seq_idx = a_row - 1

                # skip cells out of range of this chunk of original "A" matrix
                if a_row < inss[brk] or a_col < dels[brk] or \
                        a_row > inss[next_brk] or a_col > dels[next_brk]:
                    continue

                # enforce new path remains within r cells of original path
                elif b_col == 0 or b_col == 2*r:
                    for typ in range(typs):
                        matrix[typ, b_row, b_col, VAL] = INF * (b_row+1)
                        matrix[typ, b_row, b_col, TYP] = MAT
                        matrix[typ, b_row, b_col, RUN] = 0
                    continue

                # get n-polymer info
                if a_col >= a_cols - 1:
                    rpts = rpt = n = idx = 0
                else:
                    rpts = np_info[RPTS, ref_idx+1]
                    rpt = np_info[RPT, ref_idx+1]
                    n = np_info[N, ref_idx+1]
                    idx = np_info[IDX, ref_idx+1]


                # UPDATE INS MATRIX
                if a_row == inss[brk]: # first row
                    matrix[INS, b_row, b_col, VAL] = INF * (a_col-dels[brk]+1)
                    matrix[INS, b_row, b_col, TYP] = DEL
                    matrix[INS, b_row, b_col, RUN] = a_col - dels[brk]
                else:
                    val1 = matrix[MAT, b_top_row, b_top_col, VAL] + indel_start
                    matrix[INS, b_row, b_col, VAL] = val1
                    matrix[INS, b_row, b_col, TYP] = INS
                    matrix[INS, b_row, b_col, RUN] = 1

                    val2 = matrix[INS, b_top_row, b_top_col, VAL] + indel_extend
                    if val2 < val1:
                        if a_row == inss[brk] + 1:
                            run = 1
                        else:
                            run = <int>(matrix[INS, b_top_row, b_top_col, RUN]) + 1
                        matrix[INS, b_row, b_col, VAL] = val2
                        matrix[INS, b_row, b_col, TYP] = INS
                        matrix[INS, b_row, b_col, RUN] = run

                    val3 = matrix[NPI, b_top_row, b_top_col, VAL] + indel_extend
                    if val3 < val1 and val3 < val2:
                        if a_row == inss[brk] + 1:
                            run = 1
                        else:
                            run = <int>(matrix[NPI, b_top_row, b_top_col, RUN]) + 1
                        matrix[INS, b_row, b_col, VAL] = val3
                        matrix[INS, b_row, b_col, TYP] = INS
                        matrix[INS, b_row, b_col, RUN] = run


                # UPDATE DEL MATRIX
                if a_col == dels[brk]: # first col
                    matrix[DEL, b_row, b_col, VAL] = INF * (a_row-inss[brk]+1)
                    matrix[DEL, b_row, b_col, TYP] = INS
                    matrix[DEL, b_row, b_col, RUN] = a_row - inss[brk]
                else:
                    val1 = matrix[MAT, b_left_row, b_left_col, VAL] + indel_start
                    matrix[DEL, b_row, b_col, VAL] = val1
                    matrix[DEL, b_row, b_col, TYP] = DEL
                    matrix[DEL, b_row, b_col, RUN] = 1

                    val2 = matrix[DEL, b_left_row, b_left_col, VAL] + indel_extend
                    if val2 < val1:
                        if a_col == dels[brk] + 1:
                            run = 1
                        else:
                            run = <int>(matrix[DEL, b_left_row, b_left_col, RUN]) + 1
                        matrix[DEL, b_row, b_col, VAL] = val2
                        matrix[DEL, b_row, b_col, TYP] = DEL
                        matrix[DEL, b_row, b_col, RUN] = run

                    val3 = matrix[NPD, b_left_row, b_left_col, VAL] + indel_extend
                    if val3 < val1 and val3 < val2:
                        if a_col == dels[brk] + 1:
                            run = 1
                        else:
                            run = <int>(matrix[NPD, b_left_row, b_left_col, RUN]) + 1
                        matrix[DEL, b_row, b_col, VAL] = val3
                        matrix[DEL, b_row, b_col, TYP] = DEL
                        matrix[DEL, b_row, b_col, RUN] = run


                # UPDATE MAT MATRIX
                if a_row > inss[brk] and a_col > dels[brk]: # can move diag
                    if matrix[MAT, b_diag_row, b_diag_col, TYP] == MAT:
                        run = <int>(matrix[MAT, b_diag_row, b_diag_col, RUN]) + 1
                    else:
                        run = 1
                    val1 = matrix[MAT, b_diag_row, b_diag_col, VAL] + \
                            sub_scores[ seq[seq_idx], ref[ref_idx] ]
                    matrix[MAT, b_row, b_col, VAL] = val1
                    matrix[MAT, b_row, b_col, TYP] = MAT
                    matrix[MAT, b_row, b_col, RUN] = run

                else:
                    # ensure val1 isn't chosen
                    val1 = matrix[DEL, b_row, b_col, VAL] + INF

                # end INDEL
                for typ in range(1,typs): # [INS, NPI, DEL, NPD]
                    val2 = matrix[typ, b_row, b_col, VAL]
                    if val2 < val1:
                        val1 = val2
                        run = <int>(matrix[typ, b_row, b_col, RUN])
                        matrix[MAT, b_row, b_col, VAL] = val2
                        matrix[MAT, b_row, b_col, TYP] = typ
                        matrix[MAT, b_row, b_col, RUN] = run


                # UPDATE NPI MATRIX
                if a_row == inss[brk]: # first row
                    matrix[NPI, b_row, b_col, VAL] = INF * (a_col-dels[brk])
                    matrix[NPI, b_row, b_col, TYP] = DEL
                    matrix[NPI, b_row, b_col, RUN] = a_col - dels[brk]

                # start insertion
                b_ndown_row = a_to_b_row(a_row+n, a_col, inss, dels, r) - brk
                b_ndown_col = a_to_b_col(a_row+n, a_col, inss, dels, r)
                if a_row+n <= inss[next_brk] and b_ndown_col > 0: # np spans breakpoint
                    if n > 0 and idx == 0 and rpt == 0 and \
                            match(ref[ref_idx+1:ref_idx+n+1], 
                                    seq[seq_idx+1:seq_idx+1+n]):
                        val1 = matrix[MAT, b_row, b_col, VAL] + \
                                np_score(n, rpts, 1, np_scores, max_np_len)
                        matrix[NPI, b_ndown_row, b_ndown_col, VAL] = val1
                        matrix[NPI, b_ndown_row, b_ndown_col, TYP] = NPI
                        matrix[NPI, b_ndown_row, b_ndown_col, RUN] = n

                    # continue insertion
                    elif n > 0 and idx == 0 and \
                            match(ref[ref_idx+1:ref_idx+n+1], 
                                    seq[seq_idx+1:seq_idx+1+n]): 
                        run = <int>(matrix[NPI, b_row, b_col, RUN]) + n
                        b_runup_row = a_to_b_row(a_row+n-run, a_col, inss, dels, r) - brk
                        b_runup_col = a_to_b_col(a_row+n-run, a_col, inss, dels, r)
                        if run > n and a_row+n-run >= inss[brk] and b_runup_col < 2*r:
                            val2 = matrix[NPI, b_runup_row, b_runup_col, VAL] + \
                                np_score(n, rpts, <int>(run/n), np_scores, max_np_len)
                            matrix[NPI, b_ndown_row, b_ndown_col, VAL] = val2
                            matrix[NPI, b_ndown_row, b_ndown_col, TYP] = NPI
                            matrix[NPI, b_ndown_row, b_ndown_col, RUN] = run


                # UPDATE NPD MATRIX
                if a_col == dels[brk]: # first col
                    matrix[NPD, b_row, b_col, VAL] = INF * (a_row-inss[brk])
                    matrix[NPD, b_row, b_col, TYP] = INS
                    matrix[NPD, b_row, b_col, RUN] = a_row - inss[brk]

                b_nright_row = a_to_b_row(a_row, a_col+n, inss, dels, r) - brk
                b_nright_col = a_to_b_col(a_row, a_col+n, inss, dels, r)
                if a_col+n <= dels[next_brk] and b_nright_col < 2*r: # np spans breakpoint
                    if n > 0 and idx == 0 and rpt == 0: # start deletion
                        val1 = matrix[MAT, b_row, b_col, VAL] + \
                            np_score(n, rpts, -1, np_scores, max_np_len)
                        matrix[NPD, b_nright_row, b_nright_col, VAL] = val1
                        matrix[NPD, b_nright_row, b_nright_col, TYP] = NPD
                        matrix[NPD, b_nright_row, b_nright_col, RUN] = n

                    elif n > 0 and idx == 0: # continue deletion
                        run = <int>(matrix[NPD, b_row, b_col, RUN]) + n
                        b_runleft_row = a_to_b_row(a_row, a_col+n-run, inss, dels, r) - brk
                        b_runleft_col = a_to_b_col(a_row, a_col+n-run, inss, dels, r)
                        if run > n and a_col+n-run >= dels[brk] and b_runleft_col > 0:
                            val2 = matrix[NPD, b_runleft_row, b_runleft_col, VAL] + \
                                np_score(n, rpts, <int>(-run/n), np_scores, max_np_len)
                            matrix[NPD, b_nright_row, b_nright_col, VAL] = val2
                            matrix[NPD, b_nright_row, b_nright_col, TYP] = NPD
                            matrix[NPD, b_nright_row, b_nright_col, RUN] = run


        # initialize backtracking from last cell
        a_row = inss[next_brk]
        a_col = dels[next_brk]
        aln = ''
        b_row = a_to_b_row(a_row, a_col, inss, dels, r) - brk
        b_col = a_to_b_col(a_row, a_col, inss, dels, r)
        run = 0
        path = []

        # backtrack
        while a_row > inss[brk] or a_col > dels[brk]:
            b_row = a_to_b_row(a_row, a_col, inss, dels, r) - brk
            b_col = a_to_b_col(a_row, a_col, inss, dels, r)
            val = matrix[MAT, b_row, b_col, VAL]
            typ = <int>(matrix[MAT, b_row, b_col, TYP])
            run = <int>(matrix[MAT, b_row, b_col, RUN])

            # if verbose: # do some error-checking
            path.append((MAT, a_row, a_col))
            if a_row < 0:
                print(f"ERROR: row < 0 @ A:({a_row},{a_col}), B:({b_row},{b_col})")
                break
            if a_col < 0:
                print(f"ERROR: col < 0 @ A:({a_row},{a_col}), B:({b_row},{b_col})")
                break
            if run < 1:
                print(f"\nERROR: run 0 @ A:({a_row},{a_col}), B:({b_row},{b_col}),  type {typ}, val {val}")
                break

            op = ''
            if typ == NPI or typ == INS:   # each move is an insertion
                for i in range(run):
                    op += 'I' #if typ == INS else 'L'
                a_row -= run
            elif typ == NPD or typ == DEL: # each move is a deletion
                for i in range(run):
                    op += 'D' #if typ == DEL else 'S'
                a_col -= run
            elif typ == MAT: # only sub if stay in same matrix
                i = 0
                while i < run:
                    a_row -= 1
                    a_col -= 1
                    if ref[a_col-dels[brk]] == seq[a_row]:
                        op += '='
                    else:
                        op += 'X'
                    i += 1
            else:
                print("ERROR: unknown alignment matrix type:", typ)
                break
            aln += op

        full_aln += aln[::-1]

        if verbose:
            # debug print matrices
            types = ['MAT', 'INS', 'NPI', 'DEL', 'NPD']
            ops = 'MILDS'
            bases = 'NACGT'

            # PRINT A
            for typ, name in enumerate(types):
                print(f'\nA: {name}, chunk {brk_idx}')
                s = '    ~'
                for base in ref:
                    s += '        ' + bases[base]
                print(s)

                for a_row in range(inss[brk], inss[next_brk]):

                    if a_row == inss[brk]:
                        s = '~'
                    else:
                        s = bases[seq[a_row-1]]

                    for a_col in range(dels[brk], dels[next_brk]):
                        b_row = a_to_b_row(a_row, a_col, inss, dels, r) - brk
                        b_col = a_to_b_col(a_row, a_col, inss, dels, r)
                        op = ops[<int>(matrix[typ, b_row, b_col, TYP])]
                        val = f"{<int>(matrix[typ, b_row, b_col, VAL]):04}"
                        if (typ, a_row, a_col) in path:
                            mark = '$'
                        elif a_col == dels[brk] or a_row == inss[brk]:
                            mark = '*'
                        elif b_row == 0 or b_col == 2*r:
                            mark = '.'
                        else:
                            mark = ' '
                        run = <int>(matrix[typ, b_row, b_col, RUN])
                        if run < 10:
                            s += "  " + str(run) + op + mark + val
                        else:
                            s += " " + str(run) + op + mark + val
                    print(s)

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

        elif op == 'D' or op == 'S':
            ref_str += ref[ref_idx]
            ref_idx += 1
            seq_str += '-'
            cig_str += ' '

        elif op == 'I' or op == 'L':
            ref_str += '-'
            seq_str += seq[seq_idx]
            seq_idx += 1
            cig_str += ' '

        else:
            print(f"ERROR: unrecognized CIGAR operation '{op}' at cigar index {len(cig_str)}.")

    print(f"REF: len: {len(ref)} ciglen: {sum([op in 'XD=M' for op in cigar])}\n"
          f"SEQ: len: {len(seq)} ciglen: {sum([op in 'SXI=M' for op in cigar])}\n"
          f"Cigar: {cigar}\n\n")

    for x in range(0, len(cig_str), 80):
        print(ref_str[x : x+80])
        print(cig_str[x : x+80])
        print(seq_str[x : x+80])
        print(' ')

