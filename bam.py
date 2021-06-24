import multiprocessing as mp
from collections import defaultdict
import numpy as np
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt
import os, re, itertools, sys

import pysam

import cfg
from aln import *
from cig import *
from vcf import *

def realign_bam():
    '''
    Wrapper for multi-threaded read re-alignment.
    '''

    # convert BAM to workable mp format: [(id, ctg, pos, cigar, ref, seq)...]
    print('    > extracting read data from BAM')
    read_data = get_read_data()

    # update 'ref' based on SUB variant calls (optional)
    if cfg.args.splice_subs:
        print('\n    > parsing VCF')
        cfg.args.subs = get_vcf_data()
        print('\n    > splicing subs into reference')
        with mp.Pool() as pool:
            read_data = pool.map(splice_subs_into_ref(read_data))

    # align 'seq' to 'ref', update 'cigar'
    print('\n    > computing read realignments')
    with cfg.read_count.get_lock(): cfg.read_count.value = 0
    with mp.Pool() as pool:
        read_data = pool.map(realign_read, read_data)

    # standardize CIGAR format, replace subs with indels
    if cfg.args.indels_only:
        with cfg.read_count.get_lock(): cfg.read_count.value = 0
        print('\n    > converting to standard INDEL format')
        with mp.Pool() as pool:
            read_data = pool.map(standardize_cigar, read_data)

    return read_data
   


def get_read_data():

    # count reads
    bam = pysam.AlignmentFile(cfg.args.bam, 'rb')
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
    print(f'\r        0 of {nreads} reads processed.', end='', flush=True)
    for read in reads:
        read_data.append((
            read.query_name,
            read.reference_name,
            read.reference_start,
            read.cigarstring,
            read.get_reference_sequence().upper(),
            read.query_alignment_sequence.upper(),
        ))
        rds += 1
        print(f'\r        {rds} of {nreads} reads processed.', end='', flush=True)

    return read_data



def realign_read(read_data):

    read_id, ref_name, start, cigar, ref, seq = read_data
    new_cigar = align(ref, seq, ref, cigar, cfg.args.sub_scores, cfg.args.hp_scores)

    with cfg.read_count.get_lock():
        cfg.read_count.value += 1
        print(f"\r        {cfg.read_count.value} reads realigned.", end='', flush=True)

    return (read_id, ref_name, start, new_cigar, ref, seq)


    
def splice_realignments(read_data):
    '''
    Given ('read_id', [(pos1, cigar1), (pos2, cigar2)...]), return
    the full alignment with a modified CIGAR string.
    '''

    # extract useful read data
    read_id, subread_data = read_data
    read_start = cfg.starts[read_id]
    cigartuples = cfg.cigars[read_id]
    cigar_types = [ c[0] for c in cigartuples ]
    cigar_counts = [ c[1] for c in cigartuples ]
    ref_idx = read_start
    cigar = ''

    for pos, subcigar in subread_data:
        ref_start = pos - cfg.args.window
        ref_end = pos + cfg.args.window

        while ref_idx < ref_end and cigar_types:

            # store cigar prefix
            if ref_idx < ref_start:
                cigar += cfg.cigar[cigar_types[0]]

            # determine whether to move on read/ref
            if cigar_types[0] == Cigar.S:    # soft-clipped
                pass
            elif cigar_types[0] == Cigar.H:    # hard-clipped
                pass
            elif cigar_types[0] == Cigar.X:    # substitution
                ref_idx += 1
            elif cigar_types[0] == Cigar.I:    # insertion
                pass
            elif cigar_types[0] == Cigar.D:    # deletion
                ref_idx += 1
            elif cigar_types[0] == Cigar.E:    # match
                ref_idx += 1
            elif cigar_types[0] == Cigar.M:    # match/sub
                ref_idx += 1
            else:
                print(f"ERROR: unexpected CIGAR type for {read.query_name}")
                exit(1)

            # move to next CIGAR section of interest
            cigar_counts[0] -= 1
            if cigar_counts[0] == 0:
                del cigar_types[0]
                del cigar_counts[0]

        # add cigar from this subsection of read
        cigar += subcigar

    cigar += extend_pysam_cigar(cigar_types, cigar_counts)
    cigar = collapse_cigar(cigar)
    with cfg.read_count.get_lock():
        cfg.read_count.value += 1
        print(f"\r    {cfg.read_count.value} reads processed.", end='', flush=True)

    return (read_id, cigar)



def realign_pos(data):
    '''
    Re-align all reads covering a single reference column (within window).
    '''

    pos, alleles, gt = data
    alignments = []
    bam = pysam.AlignmentFile(cfg.args.bam, 'rb')

    ref_start = pos-cfg.args.window
    ref_end = pos+cfg.args.window

    for read in bam.fetch(cfg.args.contig, pos, pos+1):

        # skip read if it doesn't fully cover the window
        ref = read.get_reference_sequence().upper() \
                [ref_start-read.reference_start : ref_end-read.reference_start]
        if len(ref) < ref_end - ref_start: continue

        # find read substring overlapping region
        cigar_types = [ c[0] for c in read.cigartuples ]
        cigar_counts = [ c[1] for c in read.cigartuples ]
        read_idx, ref_idx = 0, read.reference_start
        read_start, read_end, first = None, None, True
        pre_cigar = ''

        while ref_idx < ref_end:
            cigar = cigar_types[0]

            # first read position overlapping region
            if first and ref_idx >= ref_start:
                first = False
                read_start = read_idx

            # store cigar prefix
            if ref_idx < ref_start:
                pre_cigar += cfg.cigar[cigar]

            # determine whether to move on read/ref
            if cigar == Cigar.S:    # soft-clipped
                read_idx += 1
            elif cigar == Cigar.H:    # hard-clipped
                pass
            elif cigar == Cigar.X:    # substitution
                read_idx += 1
                ref_idx += 1
            elif cigar == Cigar.I:    # insertion
                read_idx += 1
            elif cigar == Cigar.D:    # deletion
                ref_idx += 1
            elif cigar == Cigar.E:    # match
                read_idx += 1
                ref_idx += 1
            elif cigar == Cigar.M:    # match/sub
                read_idx += 1
                ref_idx += 1
            else:
                print(f"ERROR: unexpected CIGAR type for {read.query_name}")
                exit(1)

            # move to next CIGAR section of interest
            cigar_counts[0] -= 1
            if cigar_counts[0] == 0:
                del cigar_counts[0]
                del cigar_types[0]

        # extract read section
        read_end = read_idx
        seq = read.query_sequence[read_start:read_end]

        # if there's a substitution variant, alter reference
        ref_base = ref[cfg.args.window]
        hap = int(read.get_tag('HP'))-1
        if hap < 0: # unphased
            if gt[0] == gt[1]: # haplotype variants identical
                if not None in gt and len(alleles[0]) == 1 and len(alleles[gt[0]]) == 1: # sub
                    real_base = alleles[gt[0]]
                else: 
                    real_base = ref_base
            else:
                real_base = ref_base
        else:
            # TODO: combo of sub and deletion will fail here (uncommon, ignoring)
            if not None in gt and len(alleles[0]) == 1 and len(alleles[gt[hap]]) == 1: # hap sub
                real_base = alleles[gt[hap]]
            else:
                real_base = ref_base
        # print(f"position {pos}, hap {hap}, gt {gt}, alleles {alleles}, called {real_base}")
        orig_ref = ref
        ref = ref[:cfg.args.window] + real_base + ref[cfg.args.window+1:]
        
        cigar = align(ref, seq, orig_ref, cfg.args.sub_scores, cfg.args.hp_scores)
        # dump(ref, seq, cigar)
        new_cigar = standardize_cigar(cigar, ref, seq)
        alignments.append((read.query_name, pos, new_cigar))

    with cfg.pos_count.get_lock():
        cfg.pos_count.value += 1
        print(f"\r    {cfg.pos_count.value} positions processed.", end='', flush=True)

    return alignments



def write_results(alignments, outfile):
    '''
    Write a `.bam` file for a set of alignments.
    '''
    print("    > creating BAM index")
    bam = pysam.AlignmentFile(cfg.args.bam, 'rb')
    bam_index = pysam.IndexedReads(bam)
    bam_index.build()

    print("    > writing results")
    header = { 'HD': {
                   'VN': '1.6', 
                   'SO': 'coordinate'
               },
               'SQ': [{'LN': l, 'SN': ctg} for l, ctg in \
                       zip(bam.lengths, bam.references) if \
                       re.match("^chr[0-9A-Za-z][0-9a-zA-Z]?$", ctg)],
               'PG': [{
                   'PN': 'realigner',
                   'ID': 'realigner',
                   'VN': cfg.__version__,
                   'CL': ' '.join(sys.argv)
               }]
             }
    with pysam.Samfile(outfile, 'wb', header=header) as fh:

        for read_id, refname, pos, cigar, ref, seq in alignments:

            # find corresponding read in original BAM
            try:
                aln_itr = bam_index.find(read_id)
                old_alignment = next(aln_itr)
            except (KeyError, StopIteration) as e:
                print(f"ERROR: could not find read {read_id} in BAM file '{cfg.args.bam}'.")
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
            new_alignment.reference_id    = list([ctg for ctg in bam.references \
                    if re.match("^chr[0-9A-Za-z][0-9a-zA-Z]?$", ctg)]).index(cfg.args.contig)
            new_alignment.cigarstring     = cigar
            fh.write(new_alignment)

            # print progress
            with cfg.results_count.get_lock():
                cfg.results_count.value += 1
                print(f"\r    {cfg.results_count.value} of {len(alignments)} alignments written.", end='', flush=True)

    return



def get_ranges(start, stop):
    ''' Split (start, stop) into `n` even chunks. '''
    starts = list(range(start, stop, cfg.args.chunk_width))
    stops = [ min(stop, st + cfg.args.chunk_width) for st in starts ]
    return list(zip(starts, stops))



def get_confusion_matrices():
    ''' Load cached CMs if they exist. '''
    if not cfg.args.recalc_cms and \
            os.path.isfile(f'{cfg.args.stats_dir}/subs_cm.npy') and \
            os.path.isfile(f'{cfg.args.stats_dir}/hps_cm.npy'):
        print("> loading confusion matrices\r", end='')
        return np.load(f'{cfg.args.stats_dir}/subs_cm.npy'), \
                np.load(f'{cfg.args.stats_dir}/hps_cm.npy')
    else:
        print("> calculating confusion matrices")
        print(f"0 of {(cfg.args.contig_end-cfg.args.contig_beg)//cfg.args.chunk_width} chunks processed"
            f" ({cfg.args.contig}:{cfg.args.contig_beg}:{cfg.args.contig_end}).", end='', flush=True)
        ranges = get_ranges(cfg.args.contig_beg, cfg.args.contig_end)
        with mp.Pool() as pool: # multi-threaded
            results = list(pool.map(calc_confusion_matrices, ranges))

        # sum results
        sub_cms, hp_cms = list(map(np.array, zip(*results)))
        subs = np.sum(sub_cms, axis=0)
        hps = np.sum(hp_cms, axis=0)

        # cache results
        np.save(f'{cfg.args.stats_dir}/subs_cm', subs)
        np.save(f'{cfg.args.stats_dir}/hps_cm', hps)

        return subs, hps



def plot_hp_len_dists(hps):
    hps = np.sum(hps, axis=0)
    for l in range(cfg.args.max_hp):
        fig, ax = plt.subplots(1, 1, figsize=(10,10))
        hp_lens = hps[l, :] / (1 + np.sum(hps[l, :]))
        ax.step(range(cfg.args.max_hp), hp_lens, alpha=0.5, color='r')
        ax.set_title(str(l))
        plt.tight_layout()
        plt.savefig(f'{cfg.args.stats_dir}/hp{l}_dist.png', dpi=200)
        plt.close()



def plot_confusion_matrices(subs, hps):

    # plot homopolymer confusion matrices
    fig, ax = plt.subplots(2, 2, figsize=(30,30))
    cmaps = [plt.cm.Reds, plt.cm.Blues, plt.cm.Oranges, plt.cm.Greens]
    for base, idx in cfg.bases.items():
        x_idx, y_idx = idx % 2, idx // 2
        frac_matrix = hps[idx] / (1 + hps[idx].sum(axis=1))[:,np.newaxis]
        ax[x_idx, y_idx].matshow(frac_matrix, cmap=cmaps[idx], alpha=0.5)
        for i in range(cfg.args.max_hp):
            total = np.sum(hps[idx,i])
            for j in range(cfg.args.max_hp):
                count = int(hps[idx,i,j])
                frac = (count + 0.1 + int(i == j)*10) / (total + 10 +  cfg.args.max_hp/10)
                ax[x_idx, y_idx].text(x=j, y=i, 
                        s=f'{count}\n{frac*100:.1f}%\n{-np.log(frac):.2f}', 
                        va='center', ha='center')
        ax[x_idx, y_idx].set_xlabel('Predicted')
        ax[x_idx, y_idx].set_ylabel('Actual')
        ax[x_idx, y_idx].set_title(f'{base}')
    plt.suptitle('Homopolymer Confusion Matrices')
    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/hps.png', dpi=300)
    plt.close()

    # plot substitution confusion matrix
    fig, ax = plt.subplots(figsize=(5,5))
    ax.matshow(subs, cmap=plt.cm.Greys, alpha=0.5)
    for i in range(len(cfg.bases)):
        total = np.sum(subs[i])
        for j in range(len(cfg.bases)):
            count = int(subs[i,j])
            frac = (count + 0.1 + int(i==j)*10) / (total + 10 + cfg.args.max_hp/10)
            ax.text(x=j, y=i, 
                    s=f'{count}\n{frac*100:.1f}%\n{-np.log(frac):.2f}', 
                    va='center', ha='center')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    ax.set_xticks(range(len(cfg.bases)))
    ax.set_xticklabels(list(cfg.bases.keys()))
    ax.set_yticks(range(len(cfg.bases)))
    ax.set_yticklabels(list(cfg.bases.keys()))
    plt.title(f'Substitutions CM')
    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/subs.png', dpi=300)
    plt.close()



def calc_confusion_matrices(range_tuple):
    ''' Measure basecaller SUB/INDEL error profile. '''

    # initialize results matrices
    subs = np.zeros((len(cfg.bases), len(cfg.bases)))
    hps = np.zeros((len(cfg.bases), cfg.args.max_hp, cfg.args.max_hp))

    # check that BAM exists, initialize
    try:
        bam = pysam.AlignmentFile(cfg.args.bam, 'rb')
    except FileNotFoundError:
        print(f"ERROR: BAM file '{cfg.args.bam}' not found.")
        exit(1)

    # get contig length (for pysam iteration)
    try:
        contig_idx = list(bam.references).index(cfg.args.contig)
        contig_len = bam.lengths[contig_idx]
    except ValueError:
        print(f"ERROR: contig '{cfg.args.contig}' not found in '{bam}'.")
        exit(1)

    # iterate over all reference positions
    window_start, window_end = range_tuple
    for read in bam.fetch(cfg.args.contig, window_start, window_end):

        ref = read.get_reference_sequence().upper() \
                [max(0, window_start-read.reference_start) : window_end-read.reference_start]

        # find read substring overlapping region
        cigar_types = [ c[0] for c in read.cigartuples ]
        cigar_counts = [ c[1] for c in read.cigartuples ]
        read_idx, ref_idx = 0, read.reference_start
        read_start = max(read.reference_start, window_start)
        prev_cigar = None
        prev_count = 0

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
                print(f"ERROR: unexpected CIGAR type for {read.query_name}")
                exit(1)

            if ref_idx >= read_start:
                # update homopolymer confusion matrix
                hp_ptr = ref_idx - read_start

                hp_off_end = (hp_ptr+1) >= len(ref) # read ended prematurely, not INDEL
                if not hp_off_end:
                    # read just started, or new base seen (starting new HP)
                    hp_start = (hp_ptr == 0) or ref[hp_ptr] != ref[hp_ptr-1]
                    if hp_start and cigar != Cigar.S and cigar != Cigar.H: 
                        while (not hp_off_end) and ref[hp_ptr] == ref[hp_ptr+1]: # get homopolymer length
                            hp_ptr += 1
                            hp_off_end = hp_ptr+1 >= len(ref)

                        if not hp_off_end: # full HP in read
                            hp_len = hp_ptr+1 - (ref_idx-read_start)

                            # calculate INDEL length (if present)
                            indel = 0
                            if prev_cigar == Cigar.I:
                                indel = prev_count
                            elif cigar == Cigar.D:
                                indel = -count

                            # only do stats on homopolymers which fit in confusion mat
                            if not (hp_len >= cfg.args.max_hp or hp_len+indel < 0 or \
                                    hp_len+indel >= cfg.args.max_hp):
                                hps[cfg.bases[ref[ref_idx-read_start]], hp_len, hp_len+indel] += 1

            # store previous action (to detect indels directly prior to HP)
            if cigar != prev_cigar:
                prev_cigar = cigar
                prev_count = count

            # shift reference index by one base or deleted section
            if ref_move:
                if read_move:
                    if ref_idx >= read_start:
                        if ref[ref_idx-read_start] != 'N':
                            subs[ cfg.bases[ref[ref_idx-read_start]], 
                                  cfg.bases[read.query_sequence[read_idx]]] += 1
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

    with cfg.pos_count.get_lock():
        cfg.pos_count.value += 1
        print(f"\r{cfg.pos_count.value} of "
            f"{(cfg.args.contig_end-cfg.args.contig_beg)//cfg.args.chunk_width} chunks processed"
            f" ({cfg.args.contig}:{cfg.args.contig_beg}:{cfg.args.contig_end}).", end='', flush=True)

    return subs, hps

