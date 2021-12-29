import multiprocessing as mp
import os, re, itertools, sys, gc, psutil
from bisect import bisect_left, bisect_right
from time import perf_counter

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
    if cfg.args.indel_cigar:
        cigar = cigar.replace('X', 'DI').replace('=','M')
    else: 
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
        cfg.counter.value += 1
        print(f"\r    {cfg.counter.value} reads processed.", end='', flush=True)

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
    if cfg.args.indel_cigar:
        cigar = cigar.replace('X', 'DI').replace('=','M')
    else: 
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
        contig_stops = [ min(stop, st + cfg.args.chunk_width) for st in starts ]
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
        print(f"0 of {num_chunks} chunks processed", end='', flush=True)
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
            exit(1)

        return subs, nps, inss, dels



def plot_confusion_matrices(subs, nps, inss, dels, max_np_len = 20):
    ''' Generate confusion matrix plots for each n-polymer, substitutions, 
        and INDELs.
    '''

    # plot homopolymer confusion matrices
    for n in range(cfg.args.max_np):
        fig, ax = plt.subplots(figsize=(max_np_len,max_np_len))
        ax.matshow(nps[n,:max_np_len,:max_np_len] / \
                (1 + np.sum(nps[n,:max_np_len,:max_np_len],axis=1)[:, np.newaxis]), 
                cmap=plt.cm.Blues, alpha=0.5)

        # add labels
        for i in range(max_np_len):
            total = np.sum(nps[n, i, :max_np_len])
            for j in range(max_np_len):
                count = int(nps[n, i, j])
                frac = (count + 0.1 + int(i == j)*10) / (total + 10 + max_np_len*0.1)
                ax.text(x=j, y=i, 
                        s=f'{count}\n{frac*100:.1f}%\n{-np.log(frac):.2f}', 
                        va='center', ha='center')

        # formatting
        plt.ylabel('Actual')
        plt.xlabel('Predicted')
        plt.title(f'{n+1}-Polymer Confusion Matrix')
        ax.set_xticks(range(max_np_len))
        ax.set_yticks(range(max_np_len))
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
            frac = (count + 0.1 + int(i==j)*10) / (total + 10 + max_np_len*0.1)
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


    fig, ax = plt.subplots(2, 1, figsize=(max_np_len,5))
    ax[0].matshow(inss[np.newaxis,:max_np_len], cmap=plt.cm.Greens, alpha=0.5)
    ax[1].matshow(dels[np.newaxis,:max_np_len], cmap=plt.cm.Reds, alpha=0.5)

    # add labels
    total = np.sum(inss)
    for i in range(max_np_len):
        count = int(inss[i])
        frac = (count + 0.1 + int(i == j)*10) / (total + 10 + max_np_len*0.1)
        ax[0].text(x=i, y=0, 
                s=f'{count}\n{frac*100:.1f}%\n{-np.log(frac):.2f}', 
                va='center', ha='center')
    total = np.sum(dels)
    for i in range(max_np_len):
        count = int(dels[i])
        frac = (count + 0.1 + int(i == j)*10) / (total + 10 + max_np_len*0.1)
        ax[1].text(x=i, y=0, 
                s=f'{count}\n{frac*100:.1f}%\n{-np.log(frac):.2f}', 
                va='center', ha='center')

    # formatting
    ax[0].set_ylabel('INSs')
    ax[1].set_ylabel('DELs')
    ax[0].set_xticks(range(max_np_len))
    ax[1].set_xticks(range(max_np_len))
    ax[0].set_yticks([])
    ax[1].set_yticks([])
    plt.suptitle(f'INDEL Confusion Matrices')
    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/indels_cm.png', dpi=300)
    plt.close()



def calc_confusion_matrices(range_tuple):
    ''' Measure basecaller SUB/INDEL error profile. '''

    # initialize results matrices
    subs = np.zeros((cfg.nbases, cfg.nbases))
    nps = np.zeros((cfg.args.max_np, cfg.args.max_np_len, cfg.args.max_np_len))
    inss = np.zeros((cfg.args.max_np_len))
    dels = np.zeros((cfg.args.max_np_len))

    # check that BAM exists, initialize
    try:
        bam = pysam.AlignmentFile(cfg.args.bam, 'rb')
    except FileNotFoundError:
        print(f"\nERROR: BAM file '{cfg.args.bam}' not found.")
        exit(1)

    # iterate over all reference positions
    num_chunks = count_chunks(cfg.args.regions)
    contig, window_start, window_end = range_tuple
    for read in bam.fetch(contig, window_start, window_end):

        # get reference, precalculate n-polymer stats
        try:
            # ref = reference[read_start:window_end]
            ref = read.get_reference_sequence().upper() \
                    [max(0, window_start-read.reference_start) : 
                            window_end-read.reference_start]
        except TypeError: 
            # throwing error due to NoneType, checking for None fails...
            continue
        int_ref = np.zeros(len(ref), dtype=np.uint8)
        for i in range(len(ref)): 
            int_ref[i] = cfg.base_dict[ref[i]]

        np_info = get_np_info(int_ref)
        RPTS = 0
        RPT = 1
        N = 2
        IDX = 3

        # find read substring overlapping region
        cigar_types = [ c[0] for c in read.cigartuples ]
        cigar_counts = [ c[1] for c in read.cigartuples ]
        read_idx, ref_idx = 0, read.reference_start
        read_start = max(read.reference_start, window_start)
        prev_cigar = None
        prev_count = 0

        # walk along reference, keeping stats
        while ref_idx < window_end and ref_idx < read.reference_end:

            read_move, ref_move = None, None
            cigar = cigar_types[0]
            count = cigar_counts[0]

            # determine whether to move on read/ref
            if cigar == Cigar.S:    # soft-clipped
                read_move = True
                ref_move = False
            elif cigar == Cigar.H:    # hard-clipped
                read_move = False
                ref_move = False
            elif cigar == Cigar.X:    # substitution
                read_move = True
                ref_move = True
            elif cigar == Cigar.I:    # insertion
                read_move = True
                ref_move = False
            elif cigar == Cigar.D:    # deletion
                read_move = False
                ref_move = True
            elif cigar == Cigar.E:    # match
                read_move = True
                ref_move = True
            elif cigar == Cigar.M:    # match/sub
                read_move = True
                ref_move = True
            else:
                print(f"\nERROR: unexpected CIGAR type for {read.query_name}")
                exit(1)

            if ref_idx > read_start and cigar != Cigar.S and cigar != Cigar.H:

                # update SUB matrix
                if ref_move and read_move and ref[ref_idx-read_start] != 'N':
                    subs[ cfg.base_dict[ref[ref_idx-read_start]], 
                          cfg.base_dict[read.query_sequence[read_idx]]] += 1

                # update INS matrix, start of ins
                if cigar == Cigar.I and cigar != prev_cigar:
                    inss[min(count, cfg.args.max_np_len-1)] += 1
                else:
                    inss[0] += 1

                # update DEL matrix, start of del
                if cigar == Cigar.D and cigar != prev_cigar:
                    dels[min(count, cfg.args.max_np_len-1)] += 1
                else:
                    dels[0] += 1

                n = np_info[N, ref_idx-read_start]
                np_len = np_info[RPTS, ref_idx-read_start]
                idx = np_info[IDX, ref_idx-read_start]
                rpt = np_info[RPT, ref_idx-read_start]
                prev_n = np_info[N, ref_idx-read_start-1]
                prev_np_len = np_info[RPTS, ref_idx-read_start-1]

                # update N-POLYMER matrix, start of np
                if n > 0 and rpt == 0 and idx == 0:

                    if prev_cigar == Cigar.I and prev_count % n == 0:
                        if read.query_sequence[read_idx-prev_count:read_idx] == \
                                ref[ref_idx-read_start:ref_idx-read_start+n] * \
                                int(prev_count/n):
                            indel = int(prev_count / n)
                    elif cigar == Cigar.D and count % n == 0:
                        indel = - int(min(np_len, count / n))
                    else:
                        indel = 0

                    if np_len < cfg.args.max_np_len and np_len+indel<cfg.args.max_np_len:
                        nps[n-1, np_len, np_len+indel] += 1

            # store previous action (to detect indels directly prior to HP)
            if cigar != prev_cigar:
                prev_cigar = cigar
                prev_count = count

            # shift reference index by one base or deleted section
            if ref_move:
                if read_move:
                    ref_idx += 1
                else:
                    ref_idx += count

            # shift read index
            if read_move:
                cigar_counts[0] -= 1
                read_idx += 1
            else:
                cigar_counts[0] = 0

            # move to next CIGAR section of interest
            if cigar_counts[0] == 0:
                del cigar_counts[0]
                del cigar_types[0]

    with cfg.counter.get_lock():
        cfg.counter.value += 1
        print(f"\r{cfg.counter.value} of {num_chunks} chunks processed.", 
                end='', flush=True)

    return subs, nps, inss, dels



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

