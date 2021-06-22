import os, sys
import argparse

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from aln import b_to_a, a_to_b, get_inss_dels

cigar = 'DIDIIDIDDIDIDIDD'
inss, dels = get_inss_dels(cigar)
r = 2
b_rows = 17
a_rows = 8
a_cols = 10
b_cols = 5

for b_row in range(b_rows):
    for b_col in range(b_cols):
        a_row, a_col = b_to_a(b_row, b_col, inss, dels, r)
        if a_row < 0 or a_col < 0 or a_row >= a_rows or a_col >= a_cols:
            char = "~"
        elif a_row == 0 or a_col == 0:
            char = '*'
        else:
            char = '.'
        b_row2, b_col2 = a_to_b(a_row, a_col, inss, dels, r)
        print(f'({b_row2},{b_col2}) {char} '.ljust(10), end='')
    print(' ')


