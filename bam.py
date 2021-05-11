import multiprocessing as mp
from collections import defaultdict
import numpy as np
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt
import os

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



def realign_bam(positions):
    '''
    Wrapper for multi-threaded read re-alignment at each BAM position.
    '''
    
    # with mp.Pool() as pool:
    #     read_alignments = pool.map(realign_pos, positions)

    read_alignments = realign_pos(positions[0])
    


def realign_pos(pos):
    '''
    Re-align all reads covering a single reference column (within window).
    '''

    alignments = defaultdict(list)
    bam = pysam.AlignmentFile(cfg.args.bam, 'rb')

    ref_start = pos-cfg.args.window
    ref_end = pos+cfg.args.window
    print(f'pos {pos}: {ref_start}-{ref_end}')

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

        while ref_idx < ref_end:

            # first read position overlapping region
            if first and ref_idx >= ref_start:
                first = False
                read_start = read_idx

            read_move, ref_move = None, None
            cigar = cigar_types[0]

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

            # shift reference index by one base or deleted section
            if ref_move:
                if read_move:
                    ref_idx += 1
                else:
                    ref_idx += cigar_counts[0]

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

        # extract read section
        read_end = read_idx
        seq = read.query_sequence[read_start:read_end]
        print(f'\nref:\t{ref}')
        print(f'seq:\t{seq}')

        sw = LocalAlignment(IdentityScoringMatrix())
        sw.align(seq, ref).dump()

    return alignments



def write_bam(fname, alignments, header, bam=True):
    '''
    Write a `.bam` file for a set of alignments.
    '''
    with pysam.AlignmentFile(fname, 'wb', header=header) as fh:
        for ref_id, subreads in enumerate(alignments):
            for aln in sorted(subreads, key=lambda x: x.rstart):
                a = pysam.AlignedSegment()
                a.reference_id = ref_id
                a.query_name = aln.qname
                a.query_sequence = aln.seq
                a.reference_start = aln.rstart
                a.cigarstring = aln.cigar
                a.flag = aln.flag
                a.mapping_quality = 60
                fh.write(a)
    pysam.index(fname) 



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
    at_hps = hps[cfg.bases['A']] + hps[cfg.bases['T']]
    gc_hps = hps[cfg.bases['G']] + hps[cfg.bases['C']]
    at_pct = (at_hps+0.00001) / (1 + np.sum(at_hps, axis=1))[:, np.newaxis]
    gc_pct = (gc_hps+0.00001) / (1 + np.sum(gc_hps, axis=1))[:, np.newaxis]

    # 0: PERCENT HP CORRECT
    # plot percentage correctness (by hp length)
    fig, ax = plt.subplots(3, 1, figsize=(15,25))
    hp_lens = range(1, 15) # at most max_hp-1
    at_pct_correct = [at_pct[i,i] for i in hp_lens]
    gc_pct_correct = [gc_pct[i,i] for i in hp_lens]
    ax[0].plot(hp_lens, at_pct_correct, color='red', linestyle='', marker='o', alpha=0.5)
    ax[0].plot(hp_lens, gc_pct_correct, color='blue', linestyle='', marker='o', alpha=0.5)
    ax[0].set_xlabel('Actual Homopolymer Length')
    ax[0].set_ylabel('Percent Correct')
    ax[0].legend(['AT', 'GC'])

    # show polynomial fit
    at_b_cor, at_m_cor = polyfit(hp_lens, np.log(at_pct_correct), 1)
    gc_b_cor, gc_m_cor = polyfit(hp_lens, np.log(gc_pct_correct), 1)
    # ax[0].plot(hp_lens, np.exp(at_m_cor*hp_lens + at_b_cor), color='k', linestyle=':', alpha=0.5)
    # ax[0].plot(hp_lens, np.exp(gc_m_cor*hp_lens + gc_b_cor), color='k', linestyle=':', alpha=0.5)

    # 1: HP LENGTH MEAN
    # overwrite correct hp length (by interpolating), for separate distribution
    for l in hp_lens:
        at_pct[l,l] = (at_pct[l,l-1] + at_pct[l,l+1]) / 2
        gc_pct[l,l] = (gc_pct[l,l-1] + gc_pct[l,l+1]) / 2

    # estimate gaussian distributions for hp lengths
    at_means = [np.dot(at_pct[l], np.array(range(cfg.args.max_hp))) for l in hp_lens]
    gc_means = [np.dot(gc_pct[l], np.array(range(cfg.args.max_hp))) for l in hp_lens]
    ax[1].plot(hp_lens, at_means, color='red', linestyle='', marker='o', alpha=0.5)
    ax[1].plot(hp_lens, gc_means, color='blue', linestyle='', marker='o', alpha=0.5)
    at_b_mean, at_m_mean = polyfit(hp_lens, at_means, 1)
    gc_b_mean, gc_m_mean = polyfit(hp_lens, gc_means, 1)
    ax[1].plot(hp_lens, at_m_mean*hp_lens + at_b_mean, color='k', linestyle=':', alpha=0.5)
    ax[1].plot(hp_lens, gc_m_mean*hp_lens + gc_b_mean, color='k', linestyle=':', alpha=0.5)
    ax[1].set_xlabel('Actual Homopolymer Length')
    ax[1].set_ylabel('Mean Length Called')

    # 2: HP LENGTH STDDEV
    at_stds = []
    for l in hp_lens:
        vals = []
        for i in range(cfg.args.max_hp):
            vals.extend(int(1000*at_pct[l,i])*[i])
        at_stds.append(np.std(vals))
    gc_stds = []
    for l in hp_lens:
        vals = []
        for i in range(cfg.args.max_hp):
            vals.extend(int(1000*gc_pct[l,i])*[i])
        gc_stds.append(np.std(vals))
    ax[2].plot(hp_lens, at_stds, color='red', linestyle='', marker='o', alpha=0.5)
    ax[2].plot(hp_lens, gc_stds, color='blue', linestyle='', marker='o', alpha=0.5)
    at_b_std, at_m_std = polyfit(hp_lens, at_stds, 1)
    gc_b_std, gc_m_std = polyfit(hp_lens, gc_stds, 1)
    ax[2].plot(hp_lens, at_m_std*hp_lens + at_b_std, color='k', linestyle=':', alpha=0.5)
    ax[2].plot(hp_lens, gc_m_std*hp_lens + gc_b_std, color='k', linestyle=':', alpha=0.5)
    ax[2].set_xlabel('Actual Homopolymer Length')
    ax[2].set_ylabel('Standard Deviation of Length Called')

    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/best_fit.png', dpi=200)



def plot_dists(hps):
    colors = ['red', 'blue', 'orange', 'green']
    for l in range(cfg.args.max_hp):
        fig, ax = plt.subplots(1, 1, figsize=(10,10))
        for base, base_idx in cfg.bases.items():
            hp_lens = hps[base_idx, l, :] / (1 + np.sum(hps[base_idx, l, :]))
            ax.step(range(cfg.args.max_hp), hp_lens, alpha=0.5, color=colors[base_idx])
        ax.set_title(str(l))
        ax.legend(list(cfg.bases.keys()))
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

