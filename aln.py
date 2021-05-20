#!/usr/bin/env python
'''
Simple Smith-Waterman aligner
'''
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import scipy.ndimage as ndimage
import cfg


def fix_matrix_properties(scores, delta = 0.01):
    ''' Modify matrix so scores follow expected pattern. '''

    # don't penalize diagonals
    l = scores.shape[0]
    for i in range(1, l):
        scores[i,i] = min(scores[i,i], scores[i-1, i-1])

    # more insertions should be more penalized
    for j in range(1, l):
        for i in reversed(range(0, j)):
            scores[i,j] = max(
                    scores[i,j], 
                    scores[i+1,j] + delta, 
                    scores[i,j-1] + delta
            )

    # more deletions should be more penalized
    for i in range(1,l):
        for j in reversed(range(0, i)):
            scores[i,j] = max(
                    scores[i,j], 
                    scores[i,j+1] + delta, 
                    scores[i-1,j] + delta
            )

    # prefer insertions from longer homopolymers
    best = np.ones(l) * 1000
    for j in range(1,l):
        for i in reversed(range(0,j)):
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
        for j in reversed(range(0,i)):
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




class Matrix(object):
    ''' Wrapper for matrices. '''
    def __init__(self, rows, cols, init=None):
        self.rows = rows
        self.cols = cols
        self.values = [init, ] * rows * cols

    def get(self, row, col):
        return self.values[(row * self.cols) + col]

    def set(self, row, col, val):
        self.values[(row * self.cols) + col] = val


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



class Aligner(object):
    ''' Class for performing local alignment. '''
    def __init__(self, sub_scores, hp_scores, indel_start=5, indel_extend=2, verbose=False):
        ''' Set parameters for local alignment. '''
        self.sub_scores = sub_scores
        self.hp_scores = hp_scores
        self.indel_start = indel_start
        self.indel_extend = indel_extend
        self.verbose = verbose


    def hp_indel_score(self, ref_hp_len, indel_len):

        # error, don't allow
        if ref_hp_len <= 0:
            return 100
        elif ref_hp_len + indel_len < 0:
            return 100

        # force lengths to fit in matrix
        ref_hp_len = min(ref_hp_len, cfg.args.max_hp-1)
        call_hp_len = min(ref_hp_len+indel_len, cfg.args.max_hp-1)
        return self.hp_scores[ref_hp_len, call_hp_len]


    def sub_score(self, query_base, ref_base):
        return self.sub_scores[cfg.bases[ref_base], cfg.bases[query_base]]


    def align(self, ref, query, ref_name='', query_name='', rc=False):
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
                ins_val = self.indel_start if row == 1 else \
                        matrix[typ, row-1, 0, VALUE] + self.indel_extend
                matrix[typ, row, 0, :] = [ins_val, INS, 0]

        for typ in range(typs):
            for col in range(1, cols):
                del_val = self.indel_start if col == 1 else \
                        matrix[typ, 0, col-1, VALUE] + self.indel_extend
                matrix[typ, 0, col, :] = [del_val, DEL, 0]

        # calculate matrix
        for row in range(1, rows):
            for col in range(1, cols):
                ref_idx = col - 1
                query_idx = row - 1

                # UPDATE INS MATRIX
                # continue INS
                ins_run = matrix[INS, row-1, col, RUNLEN] + 1
                ins_val = matrix[INS, row-1, col, VALUE] + self.indel_extend
                matrix[INS, row, col, :] = [ins_val, INS, ins_run]

                # start INS
                min_val = ins_val
                for typ in [SUB, LHP, SHP]:
                    start_val = matrix[typ, row-1, col, VALUE] + self.indel_start
                    if start_val < min_val:
                        min_val = start_val
                        matrix[INS, row, col, :] = [min_val, typ, 1]


                # UPDATE DEL MATRIX
                # continue DEL
                del_run = matrix[DEL, row, col-1, RUNLEN] + 1
                del_val = matrix[DEL, row, col-1, VALUE] + self.indel_extend
                matrix[DEL, row, col, :] = [del_val, DEL, del_run]

                # start DEL
                min_val = del_val
                for typ in [SUB, LHP, SHP]:
                    start_val = matrix[typ, row, col-1, VALUE] + self.indel_start
                    if start_val < min_val:
                        min_val = start_val
                        matrix[DEL, row, col, :] = [min_val, typ, 1]


                # UPDATE LHP MATRIX
                if query_idx+2 < len(query) and ref_idx+1 < len(ref) and \
                        query[query_idx+1] == query[query_idx+2]: # only insert same base

                    # continue LHP
                    lhp_run = int(matrix[LHP, row-1, col, RUNLEN] + 1)
                    lhp_val = matrix[LHP, row-lhp_run, col, VALUE] + \
                            self.hp_indel_score(ref_hp_lens[ref_idx+1], lhp_run)
                    matrix[LHP, row, col, :] = [lhp_val, LHP, lhp_run]

                    # start LHP
                    # if query[query_idx+1] != query[query_idx]:
                    min_val = lhp_val
                    for typ in [SUB, DEL, INS, SHP]:
                        start_val = matrix[typ, row-1, col, VALUE] + \
                                self.hp_indel_score(ref_hp_lens[ref_idx+1], 1)
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
                            self.hp_indel_score(ref_hp_lens[ref_idx], -shp_run)
                    matrix[SHP, row, col, :] = [shp_val, SHP, shp_run]

                    # start SHP
                    # if ref_idx and ref[ref_idx] != ref[ref_idx-1]:
                    min_val = shp_val
                    for typ in [SUB, DEL, INS, LHP]:
                        start_val = matrix[typ, row, col-1, VALUE] + \
                            self.hp_indel_score(ref_hp_lens[ref_idx], -1)
                        if start_val < min_val:
                            min_val = start_val
                            matrix[SHP, row, col, :] = [min_val, typ, 1]

                else: # don't allow deleting different base
                    shp_val = max(matrix[:, row, col-1, VALUE]) + 100
                    matrix[SHP, row, col, :] = [shp_val, SUB, 0]


                # UPDATE SUB MATRIX
                # simple SUB lookup
                sub_val = matrix[SUB, row-1, col-1, VALUE] + \
                        self.sub_score(query[query_idx], ref[ref_idx])
                min_val = sub_val
                matrix[SUB, row, col, :] = [min_val, SUB, 0]

                # end INDEL
                for typ in [INS, LHP, DEL, SHP]:
                    end_val = matrix[typ, row, col, VALUE]
                    if end_val < min_val:
                        min_val = end_val
                        matrix[SUB, row, col, :] = [min_val, typ, 0]


        # initialize backtracking from last cell
        aln, path = [], []
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
                    op = 'M'
            else:
                print(f"ERROR: unknown alignment matrix type '{typ}'.")
                exit(1)

            aln.append(op)

            old_typ = new_typ

        final_row = row
        final_col = col

        # we backtracked, so get forward alignment
        aln.reverse()

        if self.verbose:
            types = ['SUB', 'INS', 'LHP', 'DEL', 'SHP']
            ops = 'MILDS'
            for typ, name in enumerate(types):
                print('\n\n', name)
                sys.stdout.write('  -    ')
                sys.stdout.write('    '.join(ref))
                sys.stdout.write('\n')
                for row in range(rows):
                    if row == 0:
                        sys.stdout.write('-')
                    else:
                        sys.stdout.write(query[row-1])

                    for col in range(cols):
                        sys.stdout.write(' %2d%s%s' % 
                                (int(matrix[typ, row, col, RUNLEN]), ops[int(matrix[typ, row, col, TYPE])], 
                                    '$' if (typ, row, col) in path else ' '))
                    sys.stdout.write('\n')

        # return alignment
        return Alignment(query, ref, final_row, final_col, _reduce_cigar(aln), 
                matrix[SUB, rows-1, cols-1, VALUE], 
                ref_name, query_name, rc)



def _reduce_cigar(operations):
    ''' Count adjacent CIGAR operations. 
        - reduces 'MMIIIMMMMDDD' to [(2, 'M'), (3, 'I'), (4, 'M'), (3, 'D')]
    '''
    count = 1
    last = None
    ret = []
    for op in operations:
        if last and op == last:
            count += 1
        elif last:
            ret.append((count, last))
            count = 1
        last = op

    if last:
        ret.append((count, last))
    return ret



def _cigar_str(cigar):
    ''' Convert CIGAR string from list of tuples to string. '''
    out = ''
    for num, op in cigar:
        out += '%s%s' % (num, op)
    return out



class Alignment(object):
    ''' Class for working with alignment results and printing. '''
    def __init__(self, query, ref, q_pos, r_pos, cigar, score, ref_name='', 
            query_name='', rc=False):
        self.query = query
        self.ref = ref
        self.q_pos = q_pos
        self.r_pos = r_pos
        self.cigar = cigar
        self.score = score
        self.r_name = ref_name
        self.q_name = query_name
        self.rc = rc

        self.r_offset = 0
        self.r_region = None

        self.query = query
        self.ref = ref

        q_len = 0
        r_len = 0

        self.matches = 0
        self.mismatches = 0

        i = self.r_pos
        j = self.q_pos

        for count, op in self.cigar:
            if op == 'M':
                q_len += count
                r_len += count
                for k in range(count):
                    if self.query[j] == self.ref[i]:
                        self.matches += 1
                    else:
                        self.mismatches += 1
                    i += 1
                    j += 1

            elif op == 'I':
                q_len += count
                j += count
                self.mismatches += count
            elif op == 'D':
                r_len += count
                i += count
                self.mismatches += count

        self.q_end = q_pos + q_len
        self.r_end = r_pos + r_len
        if self.mismatches + self.matches > 0:
            self.identity = float(self.matches) / (self.mismatches + self.matches)
        else:
            self.identity = 0


    @property
    def extended_cigar_str(self):
        ''' Return =XID CIGAR string. '''
        qpos = 0
        rpos = 0
        ext_cigar_str = ''
        working = []
        for count, op in self.cigar:
            if op == 'M':
                for k in range(count):
                    if self.query[self.q_pos + qpos + k] == self.ref[self.r_pos + rpos + k]:
                        ext_cigar_str += '='
                    else:
                        ext_cigar_str += 'X'
                qpos += count
                rpos += count

            elif op == 'I':
                qpos += count
                ext_cigar_str += 'I' * count
            elif op == 'D':
                rpos += count
                ext_cigar_str += 'D' * count

            working = _reduce_cigar(ext_cigar_str)

        out = ''
        for num, op in working:
            out += '%s%s' % (num, op)
        return out


    @property
    def cigar_str(self):
        ''' Return MID CIGAR string. '''
        return _cigar_str(self.cigar)


    def dump(self, wrap=None, out=sys.stdout):
        ''' Pretty print full alignment result. '''
        i = self.r_pos
        j = self.q_pos

        q = ''
        m = ''
        r = ''
        qlen = 0
        rlen = 0

        for count, op in self.cigar:
            if op == 'M':
                qlen += count
                rlen += count
                for k in range(count):
                    q += self.query[j]
                    r += self.ref[i]
                    if self.query[j] == self.ref[i]:
                        m += '|'
                    else:
                        m += '.'

                    i += 1
                    j += 1
            elif op == 'D':
                rlen += count
                for k in range(count):
                    q += '-'
                    r += self.ref[i]
                    m += ' '
                    i += 1
            elif op == 'I':
                qlen += count
                for k in range(count):
                    q += self.query[j]
                    r += '-'
                    m += ' '
                    j += 1

            elif op == 'N':
                q += '-//-'
                r += '-//-'
                m += '    '

        if self.q_name:
            out.write('\n\nQuery: %s%s (%s nt)\n' % (self.q_name, \
                    ' (reverse-complement)' if self.rc else '', len(self.query)))
        if self.r_name:
            if self.r_region:
                out.write('Ref  : %s (%s)\n\n' % (self.r_name, self.r_region))
            else:
                out.write('Ref  : %s (%s nt)\n\n' % (self.r_name, len(self.ref)))

        poslens = [self.q_pos + 1, self.q_end + 1, 
                self.r_pos + self.r_offset + 1, self.r_end + self.r_offset + 1]
        maxlen = max([len(str(x)) for x in poslens])

        q_pre = '\n\nQuery: %%%ss ' % maxlen
        r_pre = 'Ref  : %%%ss ' % maxlen
        m_pre = ' ' * (8 + maxlen)

        rpos = self.r_pos
        if not self.rc:
            qpos = self.q_pos
        else:
            qpos = self.q_end

        while q and r and m:
            if not self.rc:
                out.write(q_pre % (qpos + 1))  # pos is displayed as 1-based
            else:
                out.write(q_pre % (qpos))  # revcomp is 1-based on the 3' end

            if wrap:
                qfragment = q[:wrap]
                mfragment = m[:wrap]
                rfragment = r[:wrap]

                q = q[wrap:]
                m = m[wrap:]
                r = r[wrap:]
            else:
                qfragment = q
                mfragment = m
                rfragment = r

                q = ''
                m = ''
                r = ''

            out.write(qfragment)
            if not self.rc:
                for base in qfragment:
                    if base != '-':
                        qpos += 1
            else:
                for base in qfragment:
                    if base != '-':
                        qpos -= 1

            if not self.rc:
                out.write(' %s\n' % qpos)
            else:
                out.write(' %s\n' % (qpos + 1))

            out.write(m_pre)
            out.write(mfragment)
            out.write('\n')
            out.write(r_pre % (rpos + self.r_offset + 1))
            out.write(rfragment)
            for base in rfragment:
                if base != '-':
                    rpos += 1
            out.write(' %s\n\n' % (rpos + self.r_offset))

        out.write("Score: %s\n" % self.score)
        out.write("Matches: %s (%.1f%%)\n" % (self.matches, self.identity * 100))
        out.write("Mismatches: %s\n" % (self.mismatches,))
        out.write("CIGAR: %s\n" % self.cigar_str)
