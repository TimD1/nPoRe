from aln import btoa, atob

cigar = 'DIDIIDIDDIDIDIDD'
r = 2
brows = 17
arows = 8
acols = 10
bcols = 5

for brow in range(brows):
    for bcol in range(bcols):
        arow, acol = btoa(brow, bcol, cigar, r)
        if arow < 0 or acol < 0 or arow >= arows or acol >= acols:
            char = "~"
        elif arow == 0 or acol == 0:
            char = '*'
        else:
            char = '.'
        brow2, bcol2 = atob(arow, acol, cigar, r)
        print(f'({brow2},{bcol2}) {char} '.ljust(10), end='')
    print(' ')


