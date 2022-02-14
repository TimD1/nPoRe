import argparse, os, subprocess
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
import multiprocessing as mp
from collections import defaultdict

import pysam

from aln import *
from bam import *
import cfg as cfg


def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--ref", required=True,
            help="Input reference FASTA.")
    parser.add_argument("--bed", required=True,
            help="Input BED containing regions for which to compute reference "
            "n-polymer information. If '--contig*' is used, '--bed' is still "
            "required to determine complement non n-polymer regions.")

    # other options for region (must also specify BED for complement)
    parser.add_argument("--contig", type=str,
            help="Allows specifying a single contig to for which to compute "
            "n-polymer information; it can be used in combination with "
            "'--contig_beg' and '--contig_end'.")
    parser.add_argument("--contig_beg", type=int,
            help='Start of region for which n-polymers are calculated.')
    parser.add_argument("--contig_end", type=int,
            help='End of region for which n-polymers are calculated.')
    parser.add_argument("--contigs", type=str,
            help="Allows specifying multiple contigs.")

    parser.add_argument("-chunk_width", type=int, default=1000000,
            help="Reference is considered in chunks of size '--chunk_width' when "
            "calculating n-polymer information.")

    parser.add_argument("--max_n", type=int, default=6,
            help="Maximum n-polymer length (period of repeating sequence) "
            "considered.")
    parser.add_argument("--max_l", type=int, default=100,
            help="Maximum length (number of times a repeated sequence occurs)"
            " considered.")

    parser.add_argument("--out_prefix", required=True,
            help="Output BED file prefix.")

    return parser



def get_np_regions(region):
    ''' 
    For each n, return list of valid n-polymer regions (ctg, start, stop).
    Naive, will contain overlaps that are merged later.
    '''
    ctg, start, stop = region
    np_info = get_np_info(bases_to_int(cfg.args.refs[ctg][start:stop]))
    np_regions = [[] for _ in range(cfg.args.max_n)]
    L, L_IDX = 0, 1

    for pos in range(start, stop):
        idx = pos - start
        for n in range(1, cfg.args.max_n+1):
            n_idx = n - 1
            if np_info[idx, L, n_idx] and not np_info[idx, L_IDX, n_idx]: # start
                np_regions[n_idx].append((ctg, pos, pos + n * np_info[idx, L, n_idx]))

    # for idx, n in enumerate(np_regions):
    #     print(f'{ctg}:{start}-{stop}, n = {idx+1}: {len(n)}')

    return np_regions



def save_np_region_beds(np_regions, out_prefix, slop=1):

    # sort and merge BEDs by n-polymer
    print(f'> saving n-polymer BEDs, n = 1-{cfg.args.max_n}')
    start_timer = perf_counter()
    for n in range(1, cfg.args.max_n+1):

        # print all data to file
        n_idx = n-1
        with open(f'{out_prefix}_{n}.bed', 'w') as bedfile:
            for ctg_data in np_regions:
                for ctg, start, stop in ctg_data[n_idx]:
                    print(f'{ctg}\t{max(0,start-slop)}\t{stop+slop}', file=bedfile)

        # sort and merge overlaps
        merge = subprocess.Popen(["bedtools", "merge", "-i", 
                f'{out_prefix}_{n}.bed'],
                stdout = subprocess.PIPE)
        sed1 = subprocess.Popen(["sed", "s/^chr//"], stdin = merge.stdout,
                stdout = subprocess.PIPE)
        sort1 = subprocess.Popen(["sort", "-k1,1n", "-k2,2n", "-k3,3n"], 
                stdin = sed1.stdout, stdout = subprocess.PIPE)
        with open(f'{out_prefix}_{n}tmp.bed', 'w') as tmpbedfile:
            sed2 = subprocess.run(["sed", "s/^[0-9]/chr&/"], 
                    stdin = sort1.stdout, stdout = tmpbedfile)
        subprocess.run(["mv", f'{out_prefix}_{n}tmp.bed', f'{out_prefix}_{n}.bed'])
    print(f'    runtime: {perf_counter()-start_timer:.2f}s')

    # merge into single n-polymer BED
    print(f'> merging n-polymer BEDs')
    start_timer = perf_counter()
    beds = [f'{out_prefix}_{n}.bed' for n in range(1, cfg.args.max_n+1)]
    cat = subprocess.Popen(["cat"] + beds, stdout = subprocess.PIPE)
    sed3 = subprocess.Popen(["sed", "s/^chr//"], stdin = cat.stdout,
            stdout = subprocess.PIPE)
    sort2 = subprocess.Popen(["sort", "-k1,1n", "-k2,2n", "-k3,3n"], 
            stdin = sed3.stdout, stdout = subprocess.PIPE)
    sed4 = subprocess.Popen(["sed", "s/^[0-9]/chr&/"], 
            stdin = sort2.stdout, stdout = subprocess.PIPE)
    with open(f'{out_prefix}_all.bed', 'w') as allbedfile:
        subprocess.run(["bedtools", "merge"],
                stdin = sed4.stdout, stdout = allbedfile)
    print(f'    runtime: {perf_counter()-start_timer:.2f}s')

    print(f'> converting supplied .BED to .GENOME file')
    if not cfg.args.bed:
        print("ERROR: 'cfg.args.bed' must be supplied.")
        exit(1)
    elif cfg.args.bed[-4:] == ".bed":
        cfg.args.genome = f'{cfg.args.bed[:-4]}.genome'
        with open(cfg.args.genome, 'w') as genomefile:
            subprocess.run(["cut", "-f1,3", cfg.args.bed], stdout = genomefile)
    else:
        print("ERROR: 'cfg.args.bed' is not BED file.")
        exit(1)

    print(f'> finding complement')
    start_timer = perf_counter()
    with open(f'{out_prefix}_0.bed', 'w') as blandbedfile:
        subprocess.run(["bedtools", "complement", "-L", 
            "-i", f'{out_prefix}_all.bed', "-g", cfg.args.genome], 
            stdout = blandbedfile)
    print(f'    runtime: {perf_counter()-start_timer:.2f}s')




def main():

    print(f'> extracting reference contigs')
    start = perf_counter()
    get_bam_regions()
    cfg.args.refs = {}
    for ctg, _, _ in cfg.args.regions:
        cfg.args.refs[ctg] = get_fasta(cfg.args.ref, ctg)

    print(f'> subdividing into chunks')
    ranges = get_ranges(cfg.args.regions)

    print(f'> computing repeat BEDs, n = 1-{cfg.args.max_n}')
    start = perf_counter()
    with mp.Pool() as pool:
        np_regions = pool.map(get_np_regions, ranges)
    print(f'    runtime: {perf_counter()-start:.2f}s')

    save_np_region_beds(np_regions, cfg.args.out_prefix, slop=1)



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()

    try:
        main()
    except KeyboardInterrupt:
        print("\nERROR: Program terminated.")
        exit(1)
