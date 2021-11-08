import multiprocessing as mp
import os, re, itertools, sys
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

    # count reads
    try:
        bam = pysam.AlignmentFile(bam_fn, 'rb')
    except FileNotFoundError:
        print(f"\nERROR: BAM file '{bam_fn}' not found.")
        exit(1)

    reads = None
    if cfg.args.contig:
        reads = bam.fetch(
                    cfg.args.contig, 
                    cfg.args.contig_beg, 
                    cfg.args.contig_end)
    else:
        reads = bam.fetch()
    nreads = sum(1 for _ in reads)

    # get all reads in region of interest
    if cfg.args.contig:
        reads = bam.fetch(
                    cfg.args.contig, 
                    cfg.args.contig_beg, 
                    cfg.args.contig_end)
    else:
        reads = bam.fetch()

    # convert BAM to workable mp format: [(id, ctg, pos, cigar, ref, seq)...]
    read_data = []
    rds = 0
    kept = 0
    print(f'\r    0 of {nreads} reads processed.', end='', flush=True)
    for read in reads:
        if not read.is_secondary and not read.is_supplementary and not read.is_unmapped:
            read_data.append((
                read.query_name,
                read.reference_name,
                read.reference_start,
                read.reference_start + read.reference_length,
                expand_cigar(read.cigarstring).replace('S','').replace('H',''),
                read.get_reference_sequence().upper(),
                read.query_alignment_sequence.upper(),
                0 if not read.has_tag('HP') else int(read.get_tag('HP'))
            ))
            kept += 1
        rds += 1
        print(f'\r    {rds} of {nreads} reads processed, {kept} primary reads kept.', end='', flush=True)
        if cfg.args.max_reads and kept >= cfg.args.max_reads:
            break

    return read_data



def realign_read(read_data):
    '''
    Re-align reads using better estimates of SNP frequencies in new alignment
    algorithm, taking n-polymer indels into account.
    '''

    # unpack
    read_id, ref_name, start, stop, cigar, ref, seq, hap = read_data

    # convert strings to np character arrays for efficiency
    int_ref = np.zeros(len(ref), dtype=np.uint8)
    for i in range(len(ref)): 
        int_ref[i] = cfg.base_dict[ref[i]]
    int_seq = np.zeros(len(seq), dtype=np.uint8)
    for i in range(len(seq)): 
        int_seq[i] = cfg.base_dict[seq[i]]

    # align
    new_cigar = align(int_ref, int_seq, cigar, cfg.args.sub_scores, cfg.args.np_scores)

    with cfg.counter.get_lock():
        cfg.counter.value += 1
        print(f"\r    {cfg.counter.value} reads realigned.", end='', flush=True)

    return (read_id, ref_name, start, stop, new_cigar, ref, seq, hap)



def realign_read2(read_data):
    '''
    Re-align reads using SNP frequency information from other read alignments.
    '''

    read_id, ref_name, start, stop, cigar, ref, seq, hap = read_data
    padded_cigar = get_padded_cigar(cigar, start, stop)
    ref_start = 0 if start == 0 else cfg.args.pileup_positions[start-1]+1
    pad_ref = cfg.args.padded_ref[
            ref_start : 
            cfg.args.pileup_positions[stop-1]+1
    ]

    # convert strings to np character arrays for efficiency
    int_ref = np.zeros(len(pad_ref), dtype=np.uint8)
    for i in range(len(pad_ref)): 
        int_ref[i] = cfg.base_dict[pad_ref[i]]
    int_seq = np.zeros(len(seq), dtype=np.uint8)
    for i in range(len(seq)): 
        int_seq[i] = cfg.base_dict[seq[i]]

    # align
    # NOTE: considering unphased (hap0) reads to be hap1
    new_cigar = align2(int_ref, int_seq, ref_start, 
            padded_cigar, max(0, hap-1), cfg.args.pileup_scores)

    with cfg.counter.get_lock():
        cfg.counter.value += 1
        print(f"\r    {cfg.counter.value} reads realigned.", end='', flush=True)

    return (read_id, ref_name, start, stop, new_cigar, ref, seq, hap)



def get_padded_cigar(cigar, ref_start, ref_stop):
    '''
    After padding reference with '-' bases, we need to insert the appropriate
    deletions 'D' into each read's CIGAR string.
    '''

    padded_cigar = ""
    ref_idx = ref_start
    for cig in cigar:
        if cig == 'I':
            padded_cigar += cig
        elif cig in 'M=XD':
            if ref_idx == 0:
                padded_cigar += 'D' * cfg.args.pileup_positions[ref_idx]
            else:
                padded_cigar += 'D' * (cfg.args.pileup_positions[ref_idx] - \
                        cfg.args.pileup_positions[ref_idx-1] - 1)
            padded_cigar += cig
            ref_idx += 1
        else:
            print("\nERROR: unexpected CIGAR in 'get_padded_cigar()'")

    assert(ref_idx == ref_stop)
    return padded_cigar


    
def write_results(read_data, outfile):
    '''
    Write a `.bam` file for a set of alignments.
    '''
    print("> creating BAM index")
    start = perf_counter()
    bam = pysam.AlignmentFile(cfg.args.bam, 'rb')
    bam_index = pysam.IndexedReads(bam)
    bam_index.build()
    print(f'    runtime: {perf_counter()-start:.2f}s')

    print("> writing results")
    start = perf_counter()
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
    with pysam.Samfile(outfile, 'wb', header=header) as fh:

        for read_id, refname, read_start, read_end, cigar, ref, seq, hap in read_data:

            # find corresponding read in original BAM
            try:
                aln_itr = bam_index.find(read_id)
                old_alignment = next(aln_itr)
                while   old_alignment.is_secondary or \
                        old_alignment.is_supplementary or \
                        old_alignment.is_unmapped:
                    old_alignment = next(aln_itr)

            except (KeyError, StopIteration) as e:
                print(f"\nERROR: could not find primary read {read_id} in BAM file '{cfg.args.bam}'.")
                exit(1)

            # overwrite CIGAR string
            new_alignment = pysam.AlignedSegment()
            new_alignment.query_name      = old_alignment.query_name
            new_alignment.query_sequence  = old_alignment.query_alignment_sequence
            new_alignment.flag            = old_alignment.flag
            new_alignment.reference_start = old_alignment.reference_start
            new_alignment.mapping_quality = old_alignment.mapping_quality
            new_alignment.query_qualities = old_alignment.query_alignment_qualities
            new_alignment.tags            = old_alignment.tags
            if new_alignment.has_tag('MD'):
                new_alignment.set_tag('MD', None)
            new_alignment.reference_id    = bam.references.index(cfg.args.contig)
            new_alignment.cigarstring     = collapse_cigar(cigar)
            fh.write(new_alignment)

            # print progress
            with cfg.counter.get_lock():
                cfg.counter.value += 1
                print(f"\r    {cfg.counter.value} of {len(read_data)} alignments written.", end='', flush=True)

    pysam.index(outfile)
    print(f'\n    runtime: {perf_counter()-start:.2f}s')
    return



def get_ranges(start, stop):
    ''' Split (start, stop) into `n` even chunks. '''
    starts = list(range(start, stop, cfg.args.chunk_width))
    stops = [ min(stop, st + cfg.args.chunk_width) for st in starts ]
    return list(zip(starts, stops))



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
        print(f"0 of {(cfg.args.contig_end-cfg.args.contig_beg+cfg.args.chunk_width-1)//cfg.args.chunk_width} chunks processed"
            f" ({cfg.args.contig}:{cfg.args.contig_beg}-{cfg.args.contig_end}).", end='', flush=True)
        ranges = get_ranges(cfg.args.contig_beg, cfg.args.contig_end)
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

        return subs, nps, inss, dels



def get_pileup_info(verbose=False):

    cfg.args.reference = get_fasta(cfg.args.ref, cfg.args.contig)
    max_inss = get_max_inss()
    cfg.args.pileup_positions = calc_positions(max_inss)
    cfg.args.padded_ref = get_padded_reference()
    cfg.args.pileup_scores = get_pileup_scores()

    # print debug pileup counts
    if verbose:
        print("\n\nFULL PILEUP SCORE MATRIX")
        for hap in range(2):
            print("\n\t\tN\tA\tC\tG\tT\t-")
            ct = 0
            for idx, ref_base in enumerate(cfg.args.padded_ref):
                if ref_base in cfg.bases:
                    print(f"{ct}\t{ref_base}", end="")
                    ct += 1
                else: # padding '-'
                    print(f"\t{ref_base}", end="")
                for i in range(6):
                    print(f"\t{cfg.args.pileup_scores[hap,i,idx]:.2f}", end="")
                print(" ")



def get_padded_reference():
    padded_len = cfg.args.pileup_positions[-1] + 1
    padded_ref = ["-"] * padded_len
    for idx, pos in enumerate(cfg.args.pileup_positions):
        padded_ref[pos] = cfg.args.reference[idx]
    return ''.join(padded_ref)
    


def get_pileup_scores():
    ''' Load cached pileup info, or calculate from BAM. 
    '''

    if not cfg.args.recalc_pileups and os.path.isfile(f'{cfg.args.stats_dir}/pileup_scores.npy'):
        print("> loading pileup scores")
        return np.float32(np.load(f'{cfg.args.stats_dir}/pileup_scores.npy'))

    else:
        print("> calculating pileup scores")
        start = perf_counter()
        reflen = len(cfg.args.reference)
        with cfg.counter.get_lock():
            cfg.counter.value = 0
        print(f"    0 of {(reflen+cfg.args.chunk_width-1)//cfg.args.chunk_width}"
                f" chunks processed.", end='', flush=True)

        ranges = get_ranges(0, reflen)
        with mp.Pool() as pool: # multi-threaded
            results = list(pool.map(calc_pileup_scores, ranges))
        pileup_scores = np.concatenate(results, axis=2)
        print(f'\n    runtime: {perf_counter()-start:.2f}s')
        np.save(f'{cfg.args.stats_dir}/pileup_scores', pileup_scores)

    return pileup_scores



def calc_pileup_scores(range_tuple):
    '''
    Calculate negative log-likelihood of N/A/C/G/T/- per pileup column.
    '''
    BIAS = 0.25
    window_start, window_end = range_tuple
    # include preceding insertions, handling edge case at start
    if window_start == 0:
        pileup_start_idx = 0
    else:
        pileup_start_idx = cfg.args.pileup_positions[window_start-1]+1
    pileup_end_idx = cfg.args.pileup_positions[window_end-1]+1
    pileups = np.zeros((2, 6, pileup_end_idx-pileup_start_idx), 
            dtype=np.float32) + BIAS

    # for positions with zero coverage, consider ref to be a single read
    for hap in range(2):
        for pileup_idx in range(pileup_start_idx, pileup_end_idx):
            pileups[hap, 
                    cfg.base_dict[cfg.args.padded_ref[pileup_idx]], 
                    pileup_idx-pileup_start_idx] += 1

    # check that BAM exists, initialize
    try:
        bam = pysam.AlignmentFile(f'{cfg.args.out_prefix}{cfg.args.itr}.bam', 'rb')
    except FileNotFoundError:
        print(f"\nERROR: BAM file '{cfg.args.out_prefix}{cfg.args.itr}.bam' not found.")
        exit(1)

    # iterate over all reference positions
    for read in bam.fetch(cfg.args.contig, window_start, window_end):

        try:
            s = read.query_sequence[0]
        except TypeError: 
            # throwing error due to NoneType, checking for None fails...
            continue

        # ignore unphased reads for this
        hap = None
        if read.has_tag('HP'):
            hap = read.get_tag('HP')-1
        else:
            continue

        # find read substring overlapping region
        cigar_types = [ c[0] for c in read.cigartuples ]
        cigar_counts = [ c[1] for c in read.cigartuples ]
        read_idx, ref_idx = 0, read.reference_start
        read_start = max(read.reference_start, window_start)
        prev_cigar = None
        prev_count = 0
        pileup_idx = read_start - window_start

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

            if ref_idx >= read_start:
                if cigar == Cigar.M or cigar == Cigar.X or cigar == Cigar.E:
                    while pileup_idx + pileup_start_idx < \
                            cfg.args.pileup_positions[ref_idx]:
                        pileups[hap, cfg.base_dict['-'], pileup_idx] += 1
                        pileup_idx += 1
                    pileups[hap, cfg.base_dict[read.query_sequence[read_idx]], 
                            pileup_idx] += 1
                    pileup_idx += 1
                elif cigar == Cigar.I:
                    pileups[hap, cfg.base_dict[read.query_sequence[read_idx]], 
                            pileup_idx] += 1
                    pileup_idx += 1
                elif cigar == Cigar.D:
                    while pileup_idx + pileup_start_idx < \
                            cfg.args.pileup_positions[ref_idx]:
                        pileups[hap, cfg.base_dict['-'], pileup_idx] += 1
                        pileup_idx += 1
                    pileups[hap, cfg.base_dict['-'], pileup_idx] += 1
                    pileup_idx += 1

            # store previous action (to detect indels directly prior to HP)
            if cigar != prev_cigar:
                prev_cigar = cigar
                prev_count = count

            # shift
            if ref_move: ref_idx += 1
            if read_move: read_idx += 1
            cigar_counts[0] -= 1

            # move to next CIGAR section of interest
            if cigar_counts[0] == 0:
                del cigar_counts[0]
                del cigar_types[0]

    with cfg.counter.get_lock():
        cfg.counter.value += 1
        print(f"\r    {cfg.counter.value} of "
            f"{(len(cfg.args.reference)+cfg.args.chunk_width-1)//cfg.args.chunk_width}"
            f" chunks processed.", end='', flush=True)

    pileup_sum = np.sum(pileups, axis=1) + BIAS*6
    for h in range(pileups.shape[0]):
        for b in range(pileups.shape[1]):
            for i in range(pileups.shape[2]):
                pileups[h,b,i] = -np.log(pileups[h,b,i] / pileup_sum[h,i])
    return pileups



def get_max_inss():
    ''' Load cached INSs info, or calculate from BAM. 
    '''

    if not cfg.args.recalc_pileups and os.path.isfile(f'{cfg.args.stats_dir}/max_inss.npy'):
        print("> loading max insertions")
        return np.load(f'{cfg.args.stats_dir}/max_inss.npy')

    else:
        print("> calculating max insertions")
        reflen = len(cfg.args.reference)
        with cfg.counter.get_lock():
            cfg.counter.value = 0
        print(f"    0 of {(reflen+cfg.args.chunk_width-1)//cfg.args.chunk_width}"
                f" chunks processed.", end='', flush=True)

        ranges = get_ranges(0, reflen)
        start = perf_counter()
        with mp.Pool() as pool: # multi-threaded
            results = list(pool.map(calc_max_inss, ranges))
        max_inss = np.concatenate(results, axis=None)
        print(f'\n    runtime: {perf_counter()-start:.2f}s')
        np.save(f'{cfg.args.stats_dir}/max_inss', max_inss)

    return max_inss



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
    window_start, window_end = range_tuple
    for read in bam.fetch(cfg.args.contig, window_start, window_end):

        # get reference, precalculate n-polymer stats
        try:
            # ref = reference[read_start:window_end]
            ref = read.get_reference_sequence().upper() \
                    [max(0, window_start-read.reference_start) : window_end-read.reference_start]
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
        print(f"\r{cfg.counter.value} of "
            f"{(cfg.args.contig_end-cfg.args.contig_beg+cfg.args.chunk_width-1)//cfg.args.chunk_width} chunks processed"
            f" ({cfg.args.contig}:{cfg.args.contig_beg}-{cfg.args.contig_end}).", end='', flush=True)

    return subs, nps, inss, dels



def calc_max_inss(range_tuple):
    ''' Measure read 'error' profile, per pileup position. '''

    # check that BAM exists, initialize
    try:
        bam = pysam.AlignmentFile(f'{cfg.args.out_prefix}{cfg.args.itr}.bam', 'rb')
    except FileNotFoundError:
        print(f"\nERROR: BAM file '{cfg.args.out_prefix}{cfg.args.itr}.bam' not found.")
        exit(1)

    # iterate over all reference positions
    window_start, window_end = range_tuple
    max_inss = np.zeros(window_end-window_start, dtype=int)
    for read in bam.fetch(cfg.args.contig, window_start, window_end):

        # find read substring overlapping region
        cigar_types = [ c[0] for c in read.cigartuples ]
        cigar_counts = [ c[1] for c in read.cigartuples ]
        read_idx, ref_idx = 0, read.reference_start
        read_start = max(read.reference_start, window_start)
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

            if ref_idx >= read_start and cigar == Cigar.I:
                max_inss[ref_idx-window_start] = \
                        max(max_inss[ref_idx-window_start], count)

            # move on to next section
            if ref_move: ref_idx += count
            if read_move: read_idx += count
            del cigar_counts[0]
            del cigar_types[0]

    with cfg.counter.get_lock():
        cfg.counter.value += 1
        print(f"\r    {cfg.counter.value} of "
            f"{(len(cfg.args.reference)+cfg.args.chunk_width-1)//cfg.args.chunk_width}"
            f" chunks processed.", end='', flush=True)

    return max_inss



def calc_positions(inss):
    """
    Calculate the position of a reference base in the new padded reference matrix.
    """
    return np.cumsum(inss+1)-1
