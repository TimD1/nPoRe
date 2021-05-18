import multiprocessing as mp
from collections import defaultdict
import numpy as np
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt
import os, re, itertools, sys

import pysam

import cfg
from aln import *

class Cigar():
    ''' Enum for pysam's cigartuples encoding.  '''
    M = 0
    I = 1
    D = 2
    N = 3
    S = 4
    H = 5
    P = 6
    E = 7
    X = 8
    B = 9



def get_cigars(bam):

    cigars = dict()
    starts = dict()
    bam = pysam.AlignmentFile(cfg.args.bam, 'rb')

    for read in bam.fetch():
        cigars[read.query_name] = read.cigartuples
        starts[read.query_name] = read.reference_start

    return cigars, starts



def realign_bam(positions):
    '''
    Wrapper for multi-threaded read re-alignment at each BAM position.
    '''

    print('    > getting all CIGAR strings')
    cfg.cigars, cfg.starts = get_cigars(cfg.args.bam)
    
    # calculate all realignments for specified positions
    print('    > computing subread realignments')
    cfg.pos_count.value = 0
    with mp.Pool() as pool:
        print(positions[:1])
        subread_alignments = pool.map(realign_pos, positions[:1])
    subread_alignments = list(itertools.chain(*subread_alignments))

    # store results in dict, grouped by read
    print('\n    > sorting subread realignments by read')
    read_alignments = defaultdict(list)
    for read, pos, cigar in subread_alignments:
        read_alignments[read].append((pos, cigar))
    print(f'    {len(read_alignments)} reads found.')

    # update per-read CIGAR strings in parallel
    print('    > splicing together realignments')
    with mp.Pool() as pool:
        read_alignments = pool.map(splice_realignments, 
                list(read_alignments.items()))

    return read_alignments


    
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



def realign_pos(pos):
    '''
    Re-align all reads covering a single reference column (within window).
    '''

    aligner = Aligner(cfg.args.sub_scores, cfg.args.hp_scores)
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

        alignment = aligner.align(ref, seq)
        alignment.dump()
        this_cigar = extend_cigar_str(alignment.extended_cigar_str)
        alignments.append((read.query_name, pos, this_cigar))

    with cfg.pos_count.get_lock():
        cfg.pos_count.value += 1
        print(f"\r    {cfg.pos_count.value} positions processed.", end='', flush=True)


    return alignments


def read_len(extended_cigar):
    length = 0
    for op in extended_cigar:
        if op in 'SXI=M':
            length += 1
    return length

def ref_len(extended_cigar):
    length = 0
    for op in extended_cigar:
        if op in 'XD=M':
            length += 1
    return length



def extend_cigar_str(cig):
    groups = re.findall(r"(?P<len>\d+)(?P<op>\D+)", cig)
    return ''.join([int(count)*op for (count, op) in groups])


def extend_pysam_cigar(ops, counts):
    return ''.join([int(count)*cfg.cigar[op] for (count, op) in zip(counts, ops)])


def collapse_cigar(extended_cigar):
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

    out = ''
    for num, op in groups:
        out += '%s%s' % (num, op)
    return out



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

        for read_id, cigar in alignments:

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
            new_alignment.query_sequence  = old_alignment.query_sequence
            new_alignment.flag            = old_alignment.flag
            new_alignment.reference_start = old_alignment.reference_start
            new_alignment.mapping_quality = old_alignment.mapping_quality
            new_alignment.query_qualities = old_alignment.query_qualities
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
    if not cfg.args.force and \
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



def fit_curve(hps):

    # merge complement bps and convert cm to percentages
    hps = hps.sum(axis=0)
    pct = (hps+0.00001) / (1 + np.sum(hps, axis=1))[:, np.newaxis]

    # 0: PERCENT HP CORRECT
    # plot percentage correctness (by hp length)
    fig, ax = plt.subplots(3, 1, figsize=(15,25))
    hp_lens = range(1, 81) # at most max_hp-1
    all_hp_lens = range(cfg.args.max_hp)
    pct_correct = [pct[i,i] for i in hp_lens]
    ax[0].plot(hp_lens, pct_correct, color='green', linestyle='', marker='o', alpha=0.5)
    ax[0].set_xlabel('Actual Homopolymer Length')
    ax[0].set_ylabel('Log Percent Correct')

    # show polynomial fit
    pct_cor_coeff = np.polyfit(hp_lens, pct_correct, 3)
    pct_cor_fit = np.poly1d(pct_cor_coeff)
    ax[0].plot(all_hp_lens, pct_cor_fit(all_hp_lens), color='k', linestyle=':', alpha=0.5)

    # 1: HP LENGTH MEAN
    # overwrite correct hp length (by interpolating), for separate distribution
    for l in hp_lens:
        pct[l,l] = (pct[l,l-1] + pct[l,l+1]) / 2

    # estimate gaussian distributions for hp lengths
    means = [np.dot(pct[l], np.array(range(cfg.args.max_hp))) for l in hp_lens]
    ax[1].plot(hp_lens, means, color='b', linestyle='', marker='o', alpha=0.5)
    mean_coeff = np.polyfit(hp_lens, means, 2)
    mean_fit = np.poly1d(mean_coeff)
    ax[1].plot(all_hp_lens, mean_fit(all_hp_lens), color='k', linestyle=':', alpha=0.5)
    ax[1].plot(all_hp_lens, all_hp_lens, color='r', linestyle='--', alpha=0.5)
    ax[1].set_xlabel('Actual Homopolymer Length')
    ax[1].set_ylabel('Mean Length Called')

    # 2: HP LENGTH STDDEV
    stds = []
    for l in hp_lens:
        vals = []
        for i in range(cfg.args.max_hp):
            vals.extend(int(1000*pct[l,i])*[i])
        stds.append(np.std(vals))
    ax[2].plot(hp_lens, stds, color='g', linestyle='', marker='o', alpha=0.5)
    std_coeff = np.polyfit(hp_lens, stds, 3)
    std_fit = np.poly1d(std_coeff)
    ax[2].plot(all_hp_lens, std_fit(all_hp_lens), color='k', linestyle=':', alpha=0.5)
    ax[2].set_xlabel('Actual Homopolymer Length')
    ax[2].set_ylabel('Standard Deviation of Length Called')

    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/best_fit.png', dpi=200)



def plot_dists(hps):
    hps = np.sum(hps, axis=0)
    for l in range(cfg.args.max_hp):
        fig, ax = plt.subplots(1, 1, figsize=(10,10))
        hp_lens = hps[l, :] / (1 + np.sum(hps[l, :]))
        ax.step(range(cfg.args.max_hp), hp_lens, alpha=0.5, color='r')
        ax.set_title(str(l))
        plt.tight_layout()
        plt.savefig(f'{cfg.args.stats_dir}/hp{l}_dist.png', dpi=200)
        plt.close()



def show_scores(hp_scores):
    max_hp = 20
    plt.figure(figsize=(15,15))
    plt.matshow(hp_scores[:max_hp,:max_hp], cmap=plt.cm.Reds, alpha=0.5)
    for i in range(max_hp):
        for j in range(max_hp):
            plt.text(x=j, y=i, s=f'{hp_scores[j,i]:.1f}', fontsize=7, va='center', ha='center')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title('Homopolymer Score Matrix')
    plt.savefig(f'{cfg.args.stats_dir}/hp_scores.png', dpi=300)
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

