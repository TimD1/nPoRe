import numpy as np
import subprocess
import os

def collapse_cigar(extended_cigar, return_groups=False):
    ''' 
    Converts extended CIGAR ops list to normal CIGAR string.
    'DMMMII' -> '1D3M2I'
    ''' 
    count = 1 
    last = None 
    groups = [] 
    for op in extended_cigar: 
        if last and op == last: 
            count += 1
        elif last:
            groups.append((count, last))
            count = 1 
        last = op 

    if last:
        groups.append((count, last))

    if return_groups: 
        return groups 

    out = ''
    for num, op in groups:
        out += '%s%s' % (num, op) 
    return out



reflen = 1000
nreads = 1000
readlen_min = 300
readlen_max = 700

# generate reference fasta
ref_fn = "ref.fasta"
ref_fasta = open(ref_fn, "w")
int_bases = np.random.randint(4, size=reflen)
ref = ''.join(["ACGT"[i] for i in int_bases])
print(">ref", file=ref_fasta)
print(ref, file=ref_fasta)
ref_fasta.close()

# write sam header
reads_sam = open("reads.sam", 'w')
print("@HD\tVN:1.6\tSO:coordinate", file=reads_sam)
print(f"@SQ\tSN:ref\tLN:{reflen+1}", file=reads_sam)
print("@PG\tID:generate_bam.py", file=reads_sam)
reads_sam.close()
reads_sam = open("reads.sam", 'a')

# generate reads
reads_fn = "reads.fastq"
reads_fastq = open(reads_fn, "w")
for idx in range(nreads):
    readlen = np.random.randint(readlen_min, readlen_max)
    readstart = np.random.randint(reflen-readlen)
    read = ref[readstart : readstart+readlen]
    qscores = ''.join([chr(q) for q in np.random.randint(33, 127, size=readlen)])

    # save reads to fastq
    print(f"@read{idx}", file=reads_fastq)
    print(read, file=reads_fastq)
    print("+", file=reads_fastq)
    print(qscores, file=reads_fastq)

    # convert to sam
    new_read = ""
    new_qscores = ""
    cigar = ""
    for base, q in zip(read, qscores):

        # randomly mutate some positions
        sub_base = np.random.randint(100) < 3
        ins_base = np.random.randint(100) < 5
        del_base = np.random.randint(100) < 3

        if ins_base:
            new_read += "ACGT"[np.random.randint(4)]
            new_qscores += chr(np.random.randint(33,127))
            cigar += 'I'
            readlen += 1

        if sub_base: # force new base, keep same Q
            new_read += "ACGTACGT"["ACGT".find(base)+np.random.randint(1,4)]
            new_qscores += q 
            cigar += 'X'
        elif del_base: # don't allow both sub and del
            cigar += 'D'
            readlen -= 1
        else: # not sub or del, add original base
            new_read += base
            new_qscores += q
            cigar += '='
    hap = np.random.randint(2) + 1

    print(f'read{idx}\t0\tref\t{readstart+1}\t60\t{collapse_cigar(cigar)}\t*\t0\t{readlen}\t{new_read}\t{new_qscores}\tHP:i:{hap}', file=reads_sam)
reads_fastq.close()
reads_sam.close()

# convert sam to bam
path = os.getcwd()
subprocess.call([
    "../scripts/align.sh", f"{path}/reads.sam", f"{path}/ref.fasta"
])
