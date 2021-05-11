#!/usr/bin/env python
'''
Simple Smith-Waterman aligner
'''
import sys
from io import StringIO


def calc_score_matrices(subs, hps):

    # calculate homopolymer scores matrix
    hp_scores = np.zeros_like(hps)
    for base in cfg.bases.values():
        for ref_len in range(cfg.args.max_hp):
            total = np.sum(hps[base, ref_len])
            for call_len in range(cfg.args.max_hp):
                count = int(hps[base, ref_len, call_len])
                bias = 10
                frac = (count + 0.1 + int(i==j)*bias) / (total + 0.1*cfg.args.max_hp + bias)
                hp_scores[base, ref_len, call_len] = -np.log(frac)

    # calculate substitution scores matrix
    sub_scores = np.zeros_like(subs)
    for i in range(len(cfg.bases)):
        total = np.sum(subs[i])
        for j in range(len(cfg.bases)):
            count = int(subs[i,j])
            bias = 10
            frac = (count + 0.1 + int(i==j)*bias) / (total + 0.1*cfg.args.max_hp + bias)
            sub_scores[i, j] = -np.log(frac)

    return sub_scores, hp_scores



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



class Alignment(object):
    ''' Class for performing local alignment. '''
    def __init__(self, sub_scores, hp_scores, verbose=False):
        ''' Set parameters for local alignment. '''
        self.sub_scores = sub_scores
        self.hp_scores = hp_scores
        self.verbose = verbose


    def align(self, ref, query, ref_name='', query_name='', rc=False):
        ''' Perform alignment. '''

        # initialize first row/col of matrix
        matrix = Matrix( len(query)+1, len(ref)+1, (0, ' ', 0))
        for row in range(1, matrix.rows):
            val = self.gap_penalty + (row-1) * self.gap_extension_penalty
            matrix.set(row, 0, (val, 'I', row))
        for col in range(1, matrix.cols):
            val = self.gap_penalty + (row-1) * self.gap_extension_penalty
            matrix.set(0, col, (val, 'D', col))

        # calculate matrix
        VALUE = 0
        TYPE = 1
        RUNLEN = 2
        for row in range(1, matrix.rows):
            for col in range(1, matrix.cols):

                # look up score for match/sub
                sub_val = matrix.get(row-1, col-1)[VALUE] + \
                        self.scoring_matrix.score(query[row-1], ref[col-1])

                # calculate insertion score
                ins_run = 0
                if matrix.get(row - 1, col)[TYPE] == 'I':
                    ins_run = matrix.get(row-1, col)[RUNLEN]
                    ins_val = matrix.get(row-1, col)[VALUE] + \
                            self.gap_extension_penalty
                else:
                    ins_val = matrix.get(row-1, col)[VALUE] + self.gap_penalty

                # calculate deletion score
                del_run = 0
                if matrix.get(row, col-1)[TYPE] == 'D':
                    del_run = matrix.get(row, col - 1)[RUNLEN]
                    del_val = matrix.get(row, col - 1)[VALUE] + \
                            self.gap_extension_penalty
                else:
                    del_val = matrix.get(row, col-1)[VALUE] + self.gap_penalty

                # determine optimal alignment for this cell
                cell_val = min(sub_val, del_val, ins_val)
                if del_run and cell_val == del_val: # continue deletion
                    cell = (cell_val, 'D', del_run+1)
                elif ins_run and cell_val == ins_val: # continue insertion
                    cell = (cell_val, 'I', ins_run+1)
                elif cell_val == sub_val: # match/sub
                    cell = (cell_val, 'M', 0)
                elif cell_val == del_val: # start deletion
                    cell = (cell_val, 'D', 1)
                elif cell_val == ins_val: # start insertion
                    cell = (cell_val, 'I', 1)
                else:
                    cell = (0, 'X', 0) # error

                matrix.set(row, col, cell)

        # initialize backtracking from last cell
        row = matrix.rows - 1
        col = matrix.cols - 1
        val = matrix.get(row, col)[VALUE]
        op, aln, path = '', [], []

        # backtrack
        while row > 0 or col > 0:

            # update path and alignment
            val, op, runlen = matrix.get(row, col)
            path.append((row, col))
            aln.append(op)

            # update pointers for next cell
            if op == 'M':
                row -= 1
                col -= 1
            elif op == 'I':
                row -= 1
            elif op == 'D':
                col -= 1
            else:
                break
        # we backtracked, so get forward alignment
        aln.reverse()

        # debug printing
        if self.verbose:
            self.dump_matrix(ref, query, matrix, path)
            print(aln)

        # return alignment
        return Alignment(query, ref, row, col, _reduce_cigar(aln), 
                matrix.get(matrix.rows-1, matrix.cols-1)[VALUE], 
                ref_name, query_name, rc)


    def dump_matrix(self, ref, query, matrix, path, show_row=-1, show_col=-1):
        ''' Pretty print alignment matrix. '''
        sys.stdout.write('      -      ')
        sys.stdout.write('       '.join(ref))
        sys.stdout.write('\n')
        for row in range(matrix.rows):
            if row == 0:
                sys.stdout.write('-')
            else:
                sys.stdout.write(query[row - 1])

            for col in range(matrix.cols):
                if show_row == row and show_col == col:
                    sys.stdout.write('       *')
                else:
                    sys.stdout.write(' %5s%s%s' % 
                            (matrix.get(row, col)[0], matrix.get(row, col)[1], 
                                '$' if (row, col) in path else ' '))
            sys.stdout.write('\n')



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
