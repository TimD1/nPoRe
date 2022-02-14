import multiprocessing as mp
import os, re, itertools, sys, gc, psutil, io, subprocess
from bisect import bisect_left, bisect_right
from time import perf_counter
import cython

import pysam
import numpy as np
import matplotlib.pyplot as plt

import cfg
from aln import *
from cig import *
from vcf import *
from util import *


def get_read_data(bam_fn):

    try: # verify bam exists
        bam = pysam.AlignmentFile(bam_fn, 'rb')
    except FileNotFoundError:
        print(f"\nERROR: BAM file '{bam_fn}' not found.")
        exit(1)

    kept = 0
    for ctg, start, stop in cfg.args.regions:
        for read in bam.fetch(ctg, start, stop):
            if cfg.args.max_reads and kept >= cfg.args.max_reads:
                return
            if not read.is_secondary and not read.is_supplementary \
                    and not read.is_unmapped:
                kept += 1
                yield (
                    read.query_name,
                    read.flag,
                    read.reference_name,
                    read.reference_start,
                    read.mapping_quality,
                    read.cigarstring,
                    read.reference_start + read.reference_length,
                    read.query_alignment_sequence.upper(),
                    ''.join([chr(33+x) for x in read.query_alignment_qualities]),
                    read.get_reference_sequence().upper(),
                    0 if not read.has_tag('HP') else int(read.get_tag('HP'))
                )



def realign_read(read_data):
    '''
    Re-align reads using better estimates of SNP frequencies in new alignment
    algorithm, taking n-polymer indels into account.
    '''
    ### ALIGN ###
    read_id, flag, ref_name, start, mapq, cigar, stop, seq, quals, ref, hap = \
            read_data
    cigar = expand_cigar(cigar).replace('S','').replace('H','')
    int_ref = bases_to_int(ref)
    int_seq = bases_to_int(seq)
    cigar = align(int_ref, int_seq, cigar, cfg.args.sub_scores, cfg.args.np_scores)

    ### STANDARDIZE ###
    cigar = cigar.replace('X', 'M').replace('=','M')
    cig_len = len(cigar)
    nshifts_buf = np.zeros(cig_len, dtype = np.uint8)
    shiftlen_buf = np.zeros(cig_len, dtype = np.uint8)
    int_cig = cig_to_int(cigar)
    while True:
        old_cig = int_cig[:]
        I, D = 1, 2
        int_cig = push_indels_left(int_cig, int_ref, nshifts_buf, shiftlen_buf, D)
        int_cig = push_inss_thru_dels(int_cig)
        int_cig = push_indels_left(int_cig, int_seq, nshifts_buf, shiftlen_buf, I)
        int_cig = push_inss_thru_dels(int_cig)
        if same_cigar(old_cig, int_cig): break
    cigar = int_to_cig(int_cig).replace('ID','M')

    ### WRITE ###
    with cfg.counter.get_lock():
        out_bam_fh = open(f'{cfg.args.out_prefix}.sam', 'a')
        print(f"{read_id}\t{flag}\t{ref_name}\t{start+1}\t{mapq}\t{collapse_cigar(cigar)}\t*\t0\t{stop-start}\t{seq}\t{quals}\tHP:i:{hap}", file=out_bam_fh)
        out_bam_fh.close()

    # free unused RAM
    del cigar, seq, quals, int_ref, int_seq, int_cig, nshifts_buf, shiftlen_buf
    if psutil.virtual_memory().percent > 90:
        gc.collect()
    


def realign_hap(hap_data):
    '''
    Re-align reads using better estimates of SNP frequencies in new alignment
    algorithm, taking n-polymer indels into account.
    '''
    ### ALIGN ###
    contig, hap, seq, ref, cigar = hap_data
    int_ref = bases_to_int(ref)
    int_seq = bases_to_int(seq)
    cigar = align(int_ref, int_seq, cigar, cfg.args.sub_scores, cfg.args.np_scores)

    ### STANDARDIZE ###
    cigar = cigar.replace('X', 'M').replace('=','M')
    cig_len = len(cigar)
    nshifts_buf = np.zeros(cig_len, dtype = np.uint8)
    shiftlen_buf = np.zeros(cig_len, dtype = np.uint8)
    int_cig = cig_to_int(cigar)
    while True:
        old_cig = int_cig[:]
        I, D = 1, 2
        int_cig = push_indels_left(int_cig, int_ref, nshifts_buf, shiftlen_buf, D)
        int_cig = push_inss_thru_dels(int_cig)
        int_cig = push_indels_left(int_cig, int_seq, nshifts_buf, shiftlen_buf, I)
        int_cig = push_inss_thru_dels(int_cig)
        if same_cigar(old_cig, int_cig): break
    cigar = int_to_cig(int_cig).replace('ID','M')

    with cfg.counter.get_lock():
        cfg.counter.value += 1
        print(f"\r    {cfg.counter.value} reads processed.", end='', flush=True)
    return contig, hap, seq, ref, cigar

    

def create_header(outfile):
    bam = pysam.AlignmentFile(cfg.args.bam, 'rb')
    header = { 'HD': {
                   'VN': '1.6', 
                   'SO': 'coordinate'
               },
               'SQ': [{'LN': l, 'SN': ctg} for l, ctg in \
                       zip(bam.lengths, bam.references)],
               'PG': [{
                   'PN': 'realigner',
                   'ID': 'realigner',
                   'VN': cfg.__version__,
                   'CL': ' '.join(sys.argv)
               }]
             }
    fh = pysam.Samfile(outfile, 'w', header=header)
    fh.close()



def get_ranges(regions):
    ''' Split (start, stop) into `n` even chunks. '''

    contigs = []
    starts = []
    stops = []
    for contig, start, stop in regions:
        contig_starts = list(range(start, stop, cfg.args.chunk_width))
        contig_stops = [ min(stop, st + cfg.args.chunk_width) 
                for st in contig_starts ]
        starts.extend(contig_starts)
        stops.extend(contig_stops)
        contigs.extend([contig]*len(contig_starts))
    return list(zip(contigs, starts, stops))



def get_confusion_matrices():
    ''' Load cached SUB/INDEL confusion matrices if they exist. 
        Otherwise, calculate them from the provided BAM.
    '''

    if not cfg.args.recalc_cms and \
            os.path.isfile(f'{cfg.args.stats_dir}/subs_cm.npy') and \
            os.path.isfile(f'{cfg.args.stats_dir}/nps_cm.npy') and \
            os.path.isfile(f'{cfg.args.stats_dir}/inss_cm.npy') and \
            os.path.isfile(f'{cfg.args.stats_dir}/dels_cm.npy'):

        print("> loading confusion matrices")
        return np.load(f'{cfg.args.stats_dir}/subs_cm.npy'), \
                np.load(f'{cfg.args.stats_dir}/nps_cm.npy'), \
                np.load(f'{cfg.args.stats_dir}/inss_cm.npy'), \
                np.load(f'{cfg.args.stats_dir}/dels_cm.npy')

    else:
        print("> calculating confusion matrices")
        num_chunks = count_chunks(cfg.args.regions)
        print(f"    0 of {num_chunks} chunks processed", end='', flush=True)
        ranges = get_ranges(cfg.args.regions)
        with mp.Pool() as pool: # multi-threaded
            results = list(pool.map(calc_confusion_matrices, ranges))
        print(" ")

        # sum results
        sub_cms, np_cms, ins_cms, del_cms = list(map(np.array, zip(*results)))
        subs = np.sum(sub_cms, axis=0)
        nps = np.sum(np_cms, axis=0)
        inss = np.sum(ins_cms, axis=0)
        dels = np.sum(del_cms, axis=0)

        # cache results
        np.save(f'{cfg.args.stats_dir}/subs_cm', subs)
        np.save(f'{cfg.args.stats_dir}/nps_cm', nps)
        np.save(f'{cfg.args.stats_dir}/inss_cm', inss)
        np.save(f'{cfg.args.stats_dir}/dels_cm', dels)

        if cfg.args.recalc_exit:
            exit(0)

        return subs, nps, inss, dels



def plot_confusion_matrices(subs, nps, inss, dels, max_l = 10, eps=0.01):
    ''' Generate confusion matrix plots for each n-polymer, substitutions, 
        and INDELs.
    '''

    # plot homopolymer confusion matrices
    for n in range(cfg.args.max_n):
        fig, ax = plt.subplots(figsize=(max_l,max_l))
        ax.matshow(nps[n,:max_l,:max_l] / \
                (1 + np.sum(nps[n,:max_l,:max_l],axis=1)[:, np.newaxis]), 
                cmap=plt.cm.Blues, alpha=0.5)

        # add labels
        for i in range(max_l):
            total = np.sum(nps[n, i, :max_l])
            for j in range(max_l):
                count = int(nps[n, i, j])
                frac = (count + eps) / (total + eps)
                ax.text(x=j, y=i, 
                        s=f'{count}\n{frac*100:.1f}%\n{-np.log(frac):.2f}', 
                        va='center', ha='center')

        # formatting
        plt.ylabel('Actual')
        plt.xlabel('Predicted')
        plt.title(f'{n+1}-Polymer Confusion Matrix')
        ax.set_xticks(range(max_l))
        ax.set_yticks(range(max_l))
        plt.tight_layout()
        plt.savefig(f'{cfg.args.stats_dir}/{n+1}-polymer_cm.png', dpi=300)
        plt.close()

    # plot substitution confusion matrix
    fig, ax = plt.subplots(figsize=(cfg.nbases,cfg.nbases))
    ax.matshow(subs, cmap=plt.cm.Greys, alpha=0.5)

    # add labels
    for i in range(cfg.nbases):
        total = np.sum(subs[i])
        for j in range(cfg.nbases):
            count = int(subs[i,j])
            frac = (count + 0.1 + int(i==j)*10) / (total + 10 + max_l*0.1)
            ax.text(x = j, y = i, 
                    s = f'{count}\n{frac*100:.1f}%\n{-np.log(frac):.2f}', 
                    va = 'center', ha = 'center')

    # formatting
    plt.ylabel('Actual')
    plt.xlabel('Predicted')
    ax.set_xticks(range(cfg.nbases))
    ax.set_xticklabels(cfg.bases)
    ax.set_yticks(range(cfg.nbases))
    ax.set_yticklabels(cfg.bases)
    plt.title(f'Substitutions Confusion Matrix')
    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/subs_cm.png', dpi=300)
    plt.close()


    fig, ax = plt.subplots(2, 1, figsize=(max_l,5))
    ax[0].matshow(inss[np.newaxis,:max_l], cmap=plt.cm.Greens, alpha=0.5)
    ax[1].matshow(dels[np.newaxis,:max_l], cmap=plt.cm.Reds, alpha=0.5)

    # add labels
    total = np.sum(inss)
    for i in range(max_l):
        count = int(inss[i])
        frac = (count + 0.1 + int(i == j)*10) / (total + 10 + max_l*0.1)
        ax[0].text(x=i, y=0, 
                s=f'{count}\n{frac*100:.1f}%\n{-np.log(frac):.2f}', 
                va='center', ha='center')
    total = np.sum(dels)
    for i in range(max_l):
        count = int(dels[i])
        frac = (count + 0.1 + int(i == j)*10) / (total + 10 + max_l*0.1)
        ax[1].text(x=i, y=0, 
                s=f'{count}\n{frac*100:.1f}%\n{-np.log(frac):.2f}', 
                va='center', ha='center')

    # formatting
    ax[0].set_ylabel('INSs')
    ax[1].set_ylabel('DELs')
    ax[0].set_xticks(range(max_l))
    ax[1].set_xticks(range(max_l))
    ax[0].set_yticks([])
    ax[1].set_yticks([])
    plt.suptitle(f'INDEL Confusion Matrices')
    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/indels_cm.png', dpi=300)
    plt.close()



def get_pileups(bam, ctg, start, end):

    # get samtools mpileup
    pileups = subprocess.Popen(["samtools", "mpileup", 
        "-r", f'{ctg}:{start+1}-{end}', bam],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    cut = subprocess.Popen(["cut", "-f5"], 
            stdin=pileups.stdout, stdout=subprocess.PIPE)

    # return iterator
    while True:
        reads = cut.stdout.readline()
        if not reads: 
            break
        yield reads.decode('utf-8').upper().strip()



@cython.boundscheck(False)
@cython.wraparound(False)
cdef char base_to_int(str base):
    if base == 'N':
        return 0
    elif base == 'A':
        return 1
    elif base == 'C':
        return 2
    elif base == 'G':
        return 3
    elif base == 'T':
        return 4
    elif base == '-':
        return 5



@cython.boundscheck(False)
@cython.wraparound(False)
cpdef char[::1] str_to_chars(str seq):
    cdef long long seqlen = len(seq)
    chars_buf = np.zeros(seqlen, dtype=np.uint8)
    cdef char[::1] chars = chars_buf
    for i in range(seqlen):
        chars[i] = <char>(seq[i] - char('\x00'))
    return chars



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef tuple calc_confusion_matrices(range_tuple):
    ''' Measure basecaller SUB/INDEL error profile. '''

    cdef str ctg, reads, c
    cdef long long start, end, pos
    cdef int[:,:,::1] np_info
    cdef char[::1] bases
    cdef int L, L_IDX, i, num_chunks, l, l_idx, indel, n, n_idx
    cdef char ref_base
    cdef int max_n = cfg.args.max_n
    cdef int max_l = cfg.args.max_l
    cdef int nbases = cfg.nbases

    # initialize results matrices
    subs_buf = np.zeros((nbases, nbases), dtype=np.longlong)
    nps_buf = np.zeros((max_n, max_l+1, max_l+1), dtype=np.longlong)
    inss_buf = np.zeros((max_l+1), dtype=np.longlong)
    dels_buf = np.zeros((max_l+1), dtype=np.longlong)
    num_chunks = count_chunks(cfg.args.regions)

    cdef long long[:,::1] subs = subs_buf
    cdef long long[:,:,::1] nps = nps_buf
    cdef long long[::1] inss = inss_buf
    cdef long long[::1] dels = dels_buf

    # get reference FASTA
    ctg, start, end = range_tuple
    pileups = get_pileups(cfg.args.bam, ctg, start, end)

    # calculate n-polymer info for region
    np_info = get_np_info(bases_to_int(cfg.args.refs[ctg][start:end+1]))
    L, L_IDX = 0, 1

    # calculate confusion matrices
    pos = 0
    bases = bases_to_int(cfg.args.refs[ctg][start:end])
    for reads in pileups:
        was_del = was_ins = True
        ref_base = bases[pos]

        i = 0
        while i < len(reads):
            c = reads[i]

            if c == '^': # ignore start char and mapping quality
                i += 2

            elif c == '$' or c == '*': # ignore end char and deletion
                i += 1

            elif c in ['N', 'A', 'C', 'G', 'T']:  # substitution
                subs[ ref_base, base_to_int(c) ] += 1
                i += 1

                # record absence of insertions/deletions
                if not was_ins:
                    inss[0] += 1
                if not was_del:
                    dels[0] += 1
                if not was_ins and not was_del:
                    for n in range(1, max_n+1):
                        n_idx = n - 1
                        l = np_info[pos+1, L, n_idx]
                        l_idx = np_info[pos+1, L_IDX, n_idx]
                        if l != 0 and l_idx == 0:
                            nps[n_idx, l, l] += 1
                was_ins = was_del = False

            elif c == '-': # deletion
                was_del = True
                
                # get deletion length
                indel = 0
                i += 1
                c = reads[i]
                while c in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                    indel += int(c)
                    indel *= 10
                    i += 1
                    c = reads[i]
                indel //= 10

                # determine if n-polymer deletion
                cnv = False
                for n in range(1, max_n+1):
                    n_idx = n-1
                    l = np_info[pos+1, L, n_idx]
                    l_idx = np_info[pos+1, L_IDX, n_idx]
                    if l != 0 and l_idx == 0 and indel % n == 0 and indel <= l*n:
                        cnv = True
                        nps[n_idx, l, l - indel//n] += 1

                    elif l != 0 and l_idx == 0:
                        nps[n_idx, l, l] += 1

                # if not, count as general deletion
                if not cnv:
                    dels[min(max_l,indel)] += 1
                i += indel

            elif c == '+': # insertion
                was_ins = True

                # get insertion length
                indel = 0
                i += 1
                c = reads[i]
                while c in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                    indel += int(c)
                    indel *= 10
                    i += 1
                    c = reads[i]
                indel //= 10

                # determine if n-polymer insertion
                cnv = False
                for n in range(1, max_n+1):
                    n_idx = n-1
                    l = np_info[pos+1, L, n_idx]
                    l_idx = np_info[pos+1, L_IDX, n_idx]
                    if l != 0 and l_idx == 0 and indel % n == 0 and \
                            cfg.args.refs[ctg][start+pos+1:start+pos+n+1] * \
                            (indel//n) == reads[i:i+indel]:
                        cnv = True
                        nps[n_idx, l, min(max_l, l + indel//n)] += 1

                    elif l != 0 and l_idx == 0:
                        nps[n_idx, l, l] += 1

                # if not, count as general insertion
                if not cnv:
                    inss[min(max_l,indel)] += 1
                i += indel

            else:
                print(f"ERROR: unexpected character '{c}'.")
                print(f'{ctg}:{pos+start} [{ref_base}] {reads}')
                break

        # record absence of insertions/deletions for last read at pos
        if not was_ins:
            inss[0] += 1
        if not was_del:
            dels[0] += 1
        if not was_ins and not was_del:
            for n in range(1, max_n+1):
                n_idx = n - 1
                l = np_info[pos+1, L, n_idx]
                l_idx = np_info[pos+1, L_IDX, n_idx]
                if l != 0 and l_idx == 0:
                    nps[n_idx, l, l] += 1

        pos += 1

    with cfg.counter.get_lock():
        cfg.counter.value += 1
        print(f"\r    {cfg.counter.value} of {num_chunks} chunks processed.", 
                end='', flush=True)

    return subs_buf, nps_buf, inss_buf, dels_buf



def hap_to_bam(hap_data, bam_fn_pre=''):
    bam_fn = f'{bam_fn_pre}haps.bam'

    # generate header
    ln, sn = 0, 0
    for ctg, hap, seq, ref, cig in hap_data:
        ln = len(ref)
        sn = ctg
        break
    header = { 'HD': {
                   'VN': '1.6', 
                   'SO': 'coordinate'
               },
               'SQ': [{'LN': ln, 'SN': sn}],
               'PG': [{
                   'PN': 'realigner',
                   'ID': 'realigner',
                   'VN': cfg.__version__,
                   'CL': ' '.join(sys.argv)
               }]
             }

    # write aligned haplotype
    with pysam.Samfile(bam_fn, 'wb', header=header) as fh:
        for ctg, hap, seq, ref, cig in hap_data:
            new_alignment = pysam.AlignedSegment()
            new_alignment.query_name      = f"hap{hap}"
            new_alignment.query_sequence  = seq
            new_alignment.flag            = 0
            new_alignment.reference_start = 0
            new_alignment.mapping_quality = 60
            new_alignment.query_qualities = pysam.qualitystring_to_array('X'*len(seq))
            new_alignment.tags            = [('HP', hap)]
            new_alignment.reference_id    = 0
            new_alignment.cigarstring     = collapse_cigar(cig)
            fh.write(new_alignment)
    print(f"    wrote '{bam_fn}'")

