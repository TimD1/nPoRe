import numpy as np
import math
np.set_printoptions(linewidth=200)
import matplotlib.pyplot as plt
from scipy import ndimage
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

        # don't penalize diagonals, don't allow 
        for i in range(1, l):
            scores[n,0,i] = 20
            scores[n,1,i] = 20
            scores[n,2,i] = 20
            scores[n,i,i] = 0

        # more insertions should be more penalized
        for j in range(1, l):
            for i in range(j-1, -1, -1):
                scores[n,i,j] = max(
                        scores[n,i,j], 
                        scores[n,i+1,j] + delta, 
                        scores[n,i,j-1] + delta
                )

        # more deletions should be more penalized
        for i in range(4,l):
            for j in range(i-1, -1, -1):
                scores[n,i,j] = max(
                        scores[n,i,j], 
                        scores[n,i,j+1] + delta, 
                        scores[n,i-1,j] + delta
                )

        # prefer INDELs in longer npolymers
        for i in range(4,l):
            for j in range(1,l):
                if i != j:
                    scores[n,i,j] = min(
                            scores[n,i,j],
                            scores[n,i-1,j-1]-delta,
                    )

    return scores



def calc_score_matrices(subs, nps, inss, dels, eps=0.01):

    # calculate homopolymer scores matrix
    np_scores = np.zeros_like(nps, dtype=np.float32)
    for n in range(cfg.args.max_n):
        for ref_len in range(cfg.args.max_l):
            total = np.sum(nps[n, ref_len])
            for call_len in range(cfg.args.max_l):
                count = int(nps[n, ref_len, call_len])
                frac = (count + eps) / (total + eps)
                np_scores[n, ref_len, call_len] = -np.log(frac)
        # np_scores[n] = ndimage.gaussian_filter(np_scores[n], sigma=1)
    np_scores = fix_matrix_properties(np_scores)

    # calculate substitution scores matrix
    sub_scores = np.zeros((cfg.nbases,cfg.nbases), dtype=np.float32)
    for i in range(1, cfg.nbases):
        for j in range(1, cfg.nbases):
            if i != j:
                sub_scores[i, j] = -np.log( (subs[i,j]+eps) / (np.sum(subs[i])+eps) )
            else:
                sub_scores[i, j] = 0

    ins_scores = np.zeros_like(inss, dtype=np.float32)
    total = np.sum(inss)
    for l in range(cfg.args.max_l):
        frac = (inss[l] + eps) / (total + eps)
        ins_scores[l] = -np.log(frac)

    del_scores = np.zeros_like(dels, dtype=np.float32)
    total = np.sum(dels)
    for l in range(cfg.args.max_l):
        frac = (dels[l] + eps) / (total + eps)
        del_scores[l] = -np.log(frac)

    return sub_scores, np_scores, ins_scores, del_scores



def plot_np_score_matrices(nps, max_l = 50):
    for n in range(cfg.args.max_n):

        # score matrix
        med_np_len = 20
        plt.figure(figsize=(med_np_len,med_np_len))
        plt.matshow(nps[n,:med_np_len, :med_np_len], cmap='RdYlGn_r',)
        for i in range(med_np_len):
            for j in range(med_np_len):
                plt.text(x=j, y=i, s=f'{nps[n,i,j]:.1f}', fontsize=5, 
                        va='center', ha='center')
        plt.xlabel('Called')
        plt.ylabel('Actual')
        plt.xticks(range(med_np_len))
        plt.yticks(range(med_np_len))
        plt.title(f'{n+1}-Polymer Score Matrix')
        plt.savefig(f'{cfg.args.stats_dir}/{n+1}-polymer_scores.png', dpi=300)
        plt.close()

        # surface plot
        x, y = np.meshgrid(range(max_l), range(3,max_l))
        fig = plt.figure(figsize=(20,10))
        ax1 = fig.add_subplot(1,2,1, projection='3d')
        ax2 = fig.add_subplot(1,2,2, projection='3d')
        ax1.plot_trisurf(x.flatten(), -y.flatten(), 
                nps[n,3:max_l,:max_l].flatten(), 
                cmap='RdYlGn_r', edgecolor='none')
        ax2.plot_trisurf(-x.flatten(), -y.flatten(), 
                nps[n,3:max_l,:max_l].flatten(), 
                cmap='RdYlGn_r', edgecolor='none')
        plt.tight_layout()
        plt.savefig(f'{cfg.args.stats_dir}/{n+1}-polymer_scores_surface.png', dpi=200)
        plt.close()

        # best-fit plot
        fig, ax = plt.subplots(1,2, figsize=(20,10))

        for i in range(3,med_np_len):
            xs = []
            ys = []
            for j in range(i, med_np_len):
                xs.append(j-i)
                ys.append(nps[n,i,j])
            ax[0].plot(xs, ys)

        for i in range(med_np_len):
            xs = []
            ys = []
            for j in range(i, -1, -1):
                xs.append(i-j)
                ys.append(nps[n,i,j])
            ax[1].plot(xs, ys)

        ax[0].set_title('INSs')
        ax[0].set_xlabel('INS Length')
        ax[0].set_ylabel('Score')

        ax[1].legend([f'HP len {x}' for x in range(3,med_np_len)])
        ax[1].set_title('DELs')
        ax[1].set_xlabel('DEL Length')
        ax[1].set_ylabel('Score')

        plt.tight_layout()
        plt.savefig(f'{cfg.args.stats_dir}/{n+1}-polymer_scores_plot.png', dpi=200)
        plt.close()

        # contour plot
        fig = plt.figure(figsize=(20,20))
        plt.contour(x, y, nps[n, 3:max_l, :max_l], 
                levels=list(range(20)), cmap='RdYlGn_r')
        plt.tight_layout()
        plt.savefig(f'{cfg.args.stats_dir}/{n+1}-polymer_scores_contour.png', dpi=200)
        plt.close()



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef int[:,:,::1] get_np_info(char[::1] seq):
    ''' Calculate N-polymer information. 

         seq:     A T A T A T A T T T T T T A A A G C G C G C

         N = 1
         L:       0 0 0 0 0 0 0 6 6 6 6 6 6 3 3 3 0 0 0 0 0 0
         L_IDX:   0 0 0 0 0 0 0 0 1 2 3 4 5 0 1 2 0 0 0 0 0 0

         N = 2
         L:       4 3 4 3 4 3 4 0 0 0 0 0 0 0 0 0 3 0 3 0 3 0
         L_IDX:   0 0 1 1 2 2 3 0 0 0 0 0 0 0 0 0 0 0 1 0 2 0

         N = 3
         L:       0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         L_IDX:   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
         
         L = number of times the sequence is repeated
         L_IDX = 0-based index of the current repeat
         N = number of bases in repeated sequence

         Explanation: shape (len(seq), 2, N)
         4(AT), 3(TA), 6(T), 3(A), 3(GC) with some overlap.  A sequence must 
         repeat at least three times to be considered an n-polymer. Note that
         6(T) is not considered to also be 3(TT).
    '''

    cdef int seq_len = len(seq)
    cdef int max_l = cfg.args.max_l
    cdef int max_n = cfg.args.max_n
    np_info_buf = np.zeros((seq_len, 2, max_n), dtype=np.intc)
    cdef int[:,:,::1] np_info = np_info_buf
    cdef int n, l, pos, l_idx, n2, longest
    cdef int seq_idx, seq_ptr

    # define constant values for indexing into `np_info` array
    cdef int L = 0
    cdef int L_IDX = 1

    for seq_idx in range(seq_len): # iterate over sequence

        # if base 'N' (encoded as 0), skip
        if not seq[seq_idx]:
            continue

        for n in range(1, max_n+1): # check each length N-polymer
            n_idx = n-1

            # get np repeat length at this position
            l = 0
            seq_ptr = seq_idx
            while seq_ptr+n < seq_len and seq[seq_ptr] == seq[seq_ptr+n]:
                seq_ptr += 1
                if (seq_ptr-seq_idx) % n == 0: # finished n-polymer
                    l += 1
            if l: l += 1 # count first

            # save n-polymer info, if long enough
            if l > 2:
                # don't overwrite equivalent: 6T with 3(TT)
                longest = True
                for n2 in range(1, n):
                    if l*n <= np_info[seq_idx, L, n2-1]*n2: 
                        longest = False

                # don't overwrite previous: 6T with 5T with 4T...
                for l_idx in range(l):
                    pos = seq_idx + l_idx*n
                    if longest and l > np_info[pos, L, n_idx]:
                        np_info[pos, L, n_idx] = min(max_l, l)
                        np_info[pos, L_IDX, n_idx] = l_idx

    return np_info



@cython.boundscheck(False)
@cython.wraparound(False)
cdef float np_score(int n, int ref_np_len, int indel_len, float[:,:,::1] np_scores, int max_n):

    # error, don't allow
    if ref_np_len <= 0:
        return 100
    elif ref_np_len + indel_len < 0:
        return 100
    elif n < 1 or n > max_n:
        return 100

    # force lengths to fit in matrix
    cdef int call_np_len = ref_np_len + indel_len
    if ref_np_len > max_n-1:
        ref_np_len = max_n-1
    if call_np_len > max_n-1:
        call_np_len = max_n-1

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
cpdef align(char[::1] full_ref, char[::1] full_seq, str cigar, 
        float[:,::1] sub_scores, float[:,:,::1] np_scores, 
        float indel_start=5, float indel_extend=1, int max_b_rows = 20000,
        int r = 30, int verbose=0):
    ''' Perform alignment.  '''

    # convert CIGAR so that each movement is row+1 or col+1, enables easy banding
    cigar = cigar.replace('X','DI').replace('=','DI').replace('M','DI')

    # precompute offsets, breakpoints, and homopolymers
    cdef int[::1] inss = get_inss(cigar)
    cdef int[::1] dels = get_dels(cigar)
    cdef int[::1] breaks = get_breaks(max_b_rows, 
            len(full_seq) + len(full_ref) + 1, inss, dels)
    cdef int[:,:,::1] np_info, np_info_seq

    # define useful constants
    cdef int a_rows = len(full_seq) + 1
    cdef int a_cols = len(full_ref) + 1
    cdef int b_cols = 2*r + 1
    cdef int b_rows = -1

    # n-polymer info indices
    cdef int L = 0
    cdef int L_IDX = 1
    cdef int N = 2

    cdef int dims = 3 # dimensions in matrix
    cdef int VAL = 0  # value (alignment score)
    cdef int TYP = 1  # type (predecessor matrix dimension)
    cdef int RUN = 2  # run length (for indels)

    cdef int typs = 5 # types
    cdef int MAT = 0  # match/substitution
    cdef int INS = 1  # insertion
    cdef int LEN = 2  # lengthen homopolymer
    cdef int DEL = 3  # deletion
    cdef int SHR = 4  # shorten homopolymer

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
    cdef int run, typ, i, brk, next_brk, brk_idx, n, n_seq, end_idx
    cdef int[::1] l, l_idx, l_seq, l_idx_seq
    cdef int max_l = cfg.args.max_l
    cdef int max_n = cfg.args.max_n
    cdef float val1, val2

    zeros_buf = np.zeros(max_n, dtype=np.intc)
    cdef int[::1] zeros = zeros_buf
    cdef char[::1] ref, seq

    # iterate over b matrix in chunks set by breakpoints
    for brk_idx in range(len(breaks)-1):

        brk = breaks[brk_idx]
        next_brk = breaks[brk_idx+1]
        b_rows = next_brk - brk + 1
        matrix_buf.fill(0)

        # only current chunk so np_info isn't too large
        ref = full_ref[dels[brk] : dels[next_brk]+1]
        seq = full_seq[inss[brk] : inss[next_brk]+1]
        np_info = get_np_info(ref)
        np_info_seq = get_np_info(seq)

        if verbose:
            print("REF:")
            print_np_info(ref)
            print("SEQ:")
            print_np_info(seq)

        # initialize N-polymer matrices with invalid states
        for b_row in range(b_rows):
            for b_col in range(b_cols):
                a_row = b_to_a_row(b_row + brk, b_col, inss, dels, r)
                a_col = b_to_a_col(b_row + brk, b_col, inss, dels, r)
                if a_row < inss[brk] or a_col < dels[brk] or \
                        a_row > inss[next_brk] or a_col > dels[next_brk] or \
                        b_col == 0 or b_col == 2*r:
                    continue
                matrix[LEN, b_row, b_col, VAL] = INF * (a_row-inss[brk] + a_col-dels[brk])
                matrix[LEN, b_row, b_col, TYP] = MAT
                matrix[LEN, b_row, b_col, RUN] = 0
                matrix[SHR, b_row, b_col, VAL] = INF * (a_row-inss[brk] + a_col-dels[brk])
                matrix[SHR, b_row, b_col, TYP] = MAT
                matrix[SHR, b_row, b_col, RUN] = 0

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
                seq_idx = a_row - inss[brk] - 1

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
                    l = zeros_buf
                    l_idx = zeros_buf
                else:
                    l = np_info[ref_idx+1, L, :]
                    l_idx = np_info[ref_idx+1, L_IDX, :]
                if a_row >= a_rows - 1:
                    l_seq = zeros_buf
                    l_idx_seq = zeros_buf
                else:
                    l_seq = np_info_seq[seq_idx+1, L, :]
                    l_idx_seq = np_info_seq[seq_idx+1, L_IDX, :]


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
                for typ in range(1,typs): # [INS, LEN, DEL, SHR]
                    val2 = matrix[typ, b_row, b_col, VAL]
                    if val2 < val1:
                        val1 = val2
                        run = <int>(matrix[typ, b_row, b_col, RUN])
                        matrix[MAT, b_row, b_col, VAL] = val2
                        matrix[MAT, b_row, b_col, TYP] = typ
                        matrix[MAT, b_row, b_col, RUN] = run


                # UPDATE LEN MATRIX
                if a_row == inss[brk]: # first row
                    matrix[LEN, b_row, b_col, VAL] = INF * (a_col-dels[brk])
                    matrix[LEN, b_row, b_col, TYP] = DEL
                    matrix[LEN, b_row, b_col, RUN] = a_col - dels[brk]

                # start insertion
                for n in range(1,max_n+1):
                    n_idx = n - 1
                    if l[n_idx] == 0 or l_seq[n_idx] == 0 or \
                            l_idx[n_idx] != 0 or not \
                            match(seq[seq_idx+1:seq_idx+1+n], 
                                    ref[ref_idx+1:ref_idx+1+n]):
                        continue
                    b_ndown_row = a_to_b_row(a_row+n, a_col, inss, dels, r) - brk
                    b_ndown_col = a_to_b_col(a_row+n, a_col, inss, dels, r)
                    if a_row+n <= inss[next_brk] and b_ndown_col > 0: # np stays in chunk

                        if l_idx_seq[n_idx] == 0: # start insertion
                            val1 = matrix[MAT, b_row, b_col, VAL] + \
                                    np_score(n, l[n_idx], 1, np_scores, max_l)
                            if val1 < matrix[LEN, b_ndown_row, b_ndown_col, VAL]:
                                matrix[LEN, b_ndown_row, b_ndown_col, VAL] = val1
                                matrix[LEN, b_ndown_row, b_ndown_col, TYP] = LEN
                                matrix[LEN, b_ndown_row, b_ndown_col, RUN] = n

                        else: # continue insertion
                            # get matrix position of ins start
                            run = <int>(matrix[LEN, b_row, b_col, RUN])
                            b_runup_row = a_to_b_row(a_row-run, a_col, inss, dels, r) - brk
                            b_runup_col = a_to_b_col(a_row-run, a_col, inss, dels, r)

                            if run > 0 and a_row-run >= inss[brk] and b_runup_col < 2*r:
                                val1 = matrix[MAT, b_runup_row, b_runup_col, VAL] + \
                                    np_score(n, l[n_idx], <int>(run/n)+1, np_scores, max_l)
                                if val1 < matrix[LEN, b_ndown_row, b_ndown_col, VAL]:
                                    matrix[LEN, b_ndown_row, b_ndown_col, VAL] = val1
                                    matrix[LEN, b_ndown_row, b_ndown_col, TYP] = LEN
                                    matrix[LEN, b_ndown_row, b_ndown_col, RUN] = run + n


                # UPDATE SHR MATRIX
                if a_col == dels[brk]: # first col
                    matrix[SHR, b_row, b_col, VAL] = INF * (a_row-inss[brk])
                    matrix[SHR, b_row, b_col, TYP] = INS
                    matrix[SHR, b_row, b_col, RUN] = a_row - inss[brk]

                for n in range(1, max_n+1):
                    n_idx = n-1
                    if l[n_idx] == 0: continue
                    b_nright_row = a_to_b_row(a_row, a_col+n, inss, dels, r) - brk
                    b_nright_col = a_to_b_col(a_row, a_col+n, inss, dels, r)
                    if a_col+n <= dels[next_brk] and b_nright_col < 2*r: # np stays in chunk
                        if l_idx[n_idx] == 0: # start deletion
                            val1 = matrix[MAT, b_row, b_col, VAL] + \
                                np_score(n, l[n_idx], -1, np_scores, max_l)
                            if val1 < matrix[SHR, b_nright_row, b_nright_col, VAL]:
                                matrix[SHR, b_nright_row, b_nright_col, VAL] = val1
                                matrix[SHR, b_nright_row, b_nright_col, TYP] = SHR
                                matrix[SHR, b_nright_row, b_nright_col, RUN] = n

                        else: # continue deletion
                            run = <int>(matrix[SHR, b_row, b_col, RUN])
                            b_runleft_row = a_to_b_row(a_row, a_col-run, inss, dels, r) - brk
                            b_runleft_col = a_to_b_col(a_row, a_col-run, inss, dels, r)

                            if run > 0 and a_col-run >= dels[brk] and b_runleft_col > 0:
                                val1 = matrix[MAT, b_runleft_row, b_runleft_col, VAL] + \
                                    np_score(n, l[n_idx], <int>(-run/n)-1, np_scores, max_l)
                                if val1 < matrix[SHR, b_nright_row, b_nright_col, VAL]:
                                    matrix[SHR, b_nright_row, b_nright_col, VAL] = val1
                                    matrix[SHR, b_nright_row, b_nright_col, TYP] = SHR
                                    matrix[SHR, b_nright_row, b_nright_col, RUN] = run + n


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
                print(f"\nERROR: row < 0 @ A:({a_row},{a_col}), B:({b_row},{b_col})")
                with open(f'{cfg.args.out_prefix}.log', 'a+') as fh:
                    print(f'val: {val}, typ: {typ}, run: {run}, a_row: {a_row}, a_col: {a_col}, b_row: {b_row}, b_col: {b_col}', file=fh)
                    print(f'path: {path}', file=fh)
                    print(f'aln: {aln}', file=fh)
                    print(f'ref: {ref}', file=fh)
                    print(f'seq: {seq}\n', file=fh)

                break
            if a_col < 0:
                print(f"\nERROR: col < 0 @ A:({a_row},{a_col}), B:({b_row},{b_col})")
                with open(f'{cfg.args.out_prefix}.log', 'a+') as fh:
                    print(f'val: {val}, typ: {typ}, run: {run}, a_row: {a_row}, a_col: {a_col}, b_row: {b_row}, b_col: {b_col}', file=fh)
                    print(f'path: {path}', file=fh)
                    print(f'aln: {aln}', file=fh)
                    print(f'ref: {ref}', file=fh)
                    print(f'seq: {seq}\n', file=fh)
                break
            if run < 1:
                print(f"\nERROR: run 0 @ A:({a_row},{a_col}), B:({b_row},{b_col}),  type {typ}, val {val}")
                with open(f'{cfg.args.out_prefix}.log', 'a+') as fh:
                    print(f'val: {val}, typ: {typ}, run: {run}, a_row: {a_row}, a_col: {a_col}, b_row: {b_row}, b_col: {b_col}', file=fh)
                    print(f'path: {path}', file=fh)
                    print(f'aln: {aln}', file=fh)
                    print(f'ref: {ref}', file=fh)
                    print(f'seq: {seq}\n', file=fh)
                break

            op = ''
            if typ == LEN or typ == INS:   # each move is an insertion
                for i in range(run):
                    op += 'I' #if typ == INS else 'L'
                a_row -= run
            elif typ == SHR or typ == DEL: # each move is a deletion
                for i in range(run):
                    op += 'D' #if typ == DEL else 'S'
                a_col -= run
            elif typ == MAT: # only sub if stay in same matrix
                i = 0
                while i < run:
                    a_row -= 1
                    a_col -= 1
                    if ref[a_col-dels[brk]] == seq[a_row-inss[brk]]:
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
            types = ['MAT', 'INS', 'LEN', 'DEL', 'SHR']
            ops = 'MILDS'
            bases = 'NACGT'

            # PRINT A
            for typ, name in enumerate(types):
                print(f'\nA: {name}, chunk {brk_idx} ({brk}, {next_brk})')
                s = '    ~'

                for a_col in range(dels[brk], dels[next_brk]):
                    s += '        ' + bases[ref[a_col-dels[brk]]]
                print(s)

                for a_row in range(inss[brk], inss[next_brk]+1):

                    if a_row == inss[brk]:
                        s = '~'
                    else:
                        s = bases[seq[a_row-inss[brk]-1]]

                    for a_col in range(dels[brk], dels[next_brk]+1):
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
            print(" ")

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



def print_np_info(char[::1] seq):
    ''' Print sequence n-polymer annotations. '''

    np_info = get_np_info(seq)
    L, L_IDX = 0, 1

    print('bases: ', end='')
    for c in seq:
        print(f'{"NACGT"[c]} ', end='')
    print('')

    for n in range(1, cfg.args.max_n+1):
        n_idx = n-1
        print(f'n={n} l: ', end='')
        for l in np_info[:, L, n_idx]:
            print(f'{l} ', end='')
        print('')

        print('l_idx: ', end='')
        for l_idx in np_info[:, L_IDX, n_idx]:
            print(f'{l_idx} ', end='')
        print('')
    print('')
