import multiprocessing as mp
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import os

import pysam

import cfg

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


# global enum for bases ACGT
bases = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


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


        # alignments[read.alignment.query_name].append((ref_start, ref_end, cigar))

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


def plot_confusion_matrices(subs, hps):

    # plot homopolymer confusion matrices
    fig, ax = plt.subplots(2, 2, figsize=(30,30))
    cmaps = [plt.cm.Reds, plt.cm.Blues, plt.cm.Oranges, plt.cm.Greens]
    for base, idx in bases.items():
        x_idx, y_idx = idx % 2, idx // 2
        frac_matrix = hps[idx] / (1 + hps[idx].sum(axis=1))[:,np.newaxis]
        ax[x_idx, y_idx].matshow(frac_matrix, cmap=cmaps[idx], alpha=0.3)
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

    # plot substitution confusion matrix
    fig, ax = plt.subplots(figsize=(5,5))
    ax.matshow(subs, cmap=plt.cm.Greys, alpha=0.3)
    for i in range(len(bases)):
        total = np.sum(subs[i])
        for j in range(len(bases)):
            count = int(subs[i,j])
            frac = (count + 0.1 + int(i==j)*10) / (total + 10 + cfg.args.max_hp/10)
            ax.text(x=j, y=i, 
                    s=f'{count}\n{frac*100:.1f}%\n{-np.log(frac):.2f}', 
                    va='center', ha='center')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    ax.set_xticks(range(len(bases)))
    ax.set_xticklabels(list(bases.keys()))
    ax.set_yticks(range(len(bases)))
    ax.set_yticklabels(list(bases.keys()))
    plt.title(f'Substitutions CM')
    plt.tight_layout()
    plt.savefig(f'{cfg.args.stats_dir}/subs.png', dpi=300)



def calc_confusion_matrices(range_tuple):
    ''' Measure basecaller SUB/INDEL error profile. '''

    # initialize results matrices
    subs = np.zeros((len(bases), len(bases)))
    hps = np.zeros((len(bases), cfg.args.max_hp, cfg.args.max_hp))

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
                                hps[bases[ref[ref_idx-read_start]], hp_len, hp_len+indel] += 1

            # store previous action (to detect indels directly prior to HP)
            if cigar != prev_cigar:
                prev_cigar = cigar
                prev_count = count

            # shift reference index by one base or deleted section
            if ref_move:
                if read_move:
                    if ref_idx >= read_start:
                        subs[ bases[ref[ref_idx-read_start]], 
                              bases[read.query_sequence[read_idx]]] += 1
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

