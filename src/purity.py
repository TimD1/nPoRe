import argparse, os, subprocess
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
import multiprocessing as mp
from collections import defaultdict

import pysam


def compute_purity(pileup_str):

    pileup_str = pileup_str.decode('utf-8')
    bases = defaultdict(int)
    inss = defaultdict(int)

    i = 0
    while i < len(pileup_str):
        c = pileup_str[i]

        if c == '^': # ignore start char and mapping quality
            i += 2

        elif c == '$': # ignore end char
            i += 1

        elif c in 'ACGT*': # add to dict
            bases[c] += 1
            i += 1

        elif c == '-': # deletion after, ignore
            skip = 0
            i += 1
            c = pileup_str[i]
            while c in '0123456789':
                skip += int(c)
                skip *= 10
                i += 1
                c = pileup_str[i]
            skip /= 10
            i += int(skip) # ignore 'N' bases

        elif c == '+': # insertion after, add to dict
            skip = 0
            i += 1
            c = pileup_str[i]
            while c in '0123456789':
                skip += int(c)
                skip *= 10
                i += 1
                c = pileup_str[i]
            skip /= 10
            inss[ pileup_str[i:i+int(skip)] ] += 1
            i += int(skip) # skip over insertion

        else:
            print(f"ERROR: unexpected character '{c}'.")
            print(pileup_str)
            break

    # calculate bases score
    n = sum(bases.values())
    if not n: return None # avoid zero-coverage
    bases_score = 0
    for b in 'ACGT*':
        bases_score += (bases[b]/n)**2

    # calculate insertion score
    not_inss = n - sum(inss.values())
    inss_score = (not_inss/n)**2
    for v in inss.values():
        inss_score += (v/n)**2

    # # debug print
    # print(pileup_str)
    # print("Bases:", bases_score)
    # for k,v in bases.items():
    #     print(k,v)
    # print("INSs:", inss_score)
    # for k,v in inss.items():
    #     print(k,v)
    # print(" ")

    return bases_score, inss_score



def plot_purity(bam_scores, out):

    plt.rcParams.update({'font.size': 22})
    fig, ax = plt.subplots(2,2,figsize=(20,8))
    labels = ['clair3-hap 1', 'clair3-hap 2', 'clair3-npore-hap 1', 'clair3-npore-hap 2']
    colors = 'rygb'

    prev_base_scores, prev_ins_scores = None, None
    prev_base_counts, prev_ins_counts = None, None
    for bam_idx, scores in enumerate(bam_scores):
        base_scores, ins_scores = list(zip(*scores))

        if bam_idx % 2:
            base_counts = [0]*100
            for x in base_scores:
                base_counts[int(x*100-0.00001)] += 1
            for x in prev_base_scores:
                base_counts[int(x*100-0.00001)] += 1

            ins_counts = [0]*100
            for x in ins_scores:
                ins_counts[int(x*100-0.00001)] += 1
            for x in prev_ins_scores:
                ins_counts[int(x*100-0.00001)] += 1

            if bam_idx == 3:

                ax[1][0].bar(np.linspace(0-0.005,1-0.005,100), [0 if not x or not y else x/y for x,y in zip(base_counts,prev_base_counts)], width=0.01)
                ax[1][0].axhline(1, color='k', linestyle=':')
                ax[1][0].set_xlim(0,1)
                ax[1][0].set_ylabel('Ratio')
                ax[1][0].set_title('Ratio: $\\frac{\\textrm{clair3-npore-hap}}{\\textrm{clair3-hap}}$')

                ax[1][1].bar(np.linspace(0-0.005,1-0.005,100), [0 if not x or not y else x/y for x,y in zip(ins_counts,prev_ins_counts)], width=0.01)
                ax[1][1].set_xlim(0,1)
                ax[1][1].axhline(1, color='k', linestyle=':')
                ax[1][1].set_title('Ratio: $\\frac{\\textrm{clair3-npore-hap}}{\\textrm{clair3-hap}}$')


        ax[0][0].hist(base_scores, bins=np.linspace(0,1,100), linewidth=3,
                histtype='step', color=colors[bam_idx], alpha=0.8)
        ax[0][0].set_xticks(np.linspace(0,1,11))
        ax[0][0].set_xlim(0,1)
        ax[0][0].set_yscale('log')
        ax[0][0].set_title('Pileup Gini Purity Histogram')
        ax[0][0].set_ylabel('Counts')

        ax[0][1].hist(ins_scores, bins=np.linspace(0,1,100), linewidth=3,
                histtype='step', color=colors[bam_idx], alpha=0.8)
        ax[0][1].set_xticks(np.linspace(0,1,11))
        ax[0][1].set_xlim(0,1)
        ax[0][1].set_yscale('log')
        ax[0][1].set_title('Insertion Gini Purity Histogram')

        prev_base_scores = base_scores
        prev_ins_scores = ins_scores

        if bam_idx % 2:
            prev_base_counts = base_counts
            prev_ins_counts = ins_counts

    ax[0][1].legend(labels)
    plt.tight_layout()
    plt.savefig(f'{out}.png', dpi=300)
    plt.close()
    


def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bams", nargs=4, required=True,
            help="Input BAMs from which to calculate Gini purity of pileups. "
            " Four BAMS are expected, in this order: baseline hap1, baseline "
            "hap2, realigned hap1, realigned hap2")
    parser.add_argument("--region", type=str,
            help="Region for which to compute pileup Gini purity.")
    parser.add_argument("--out", default="out",
            help="Output prefix for where to save calculations and plots.")
    parser.add_argument("--plot_only", action='store_true',
            help="Load cached calculations, and just re-plot.")
    return parser



def main():
    bam_scores = []

    if not args.plot_only:
        for bam_idx, bam in enumerate(args.bams):
            print(f'> computing samtools mpileup for {bam}')
            start = perf_counter()
            mpileup_args = ["-r", args.region] if args.region else []
            ps = subprocess.Popen(["samtools", "mpileup"] + mpileup_args + [bam], 
                    stdout=subprocess.PIPE)
            pileups = subprocess.check_output(["cut", "-f5"], stdin=ps.stdout).upper().split()
            ps.wait()
            print(f'    runtime: {perf_counter()-start:.2f}s')

            print(f'> computing purity for {bam}')
            start = perf_counter()
            with mp.Pool() as pool:
                bam_scores.append(list(filter(None, pool.map(compute_purity, pileups, chunksize=100))))
                pool.close()
                pool.join()
            print(f'    runtime: {perf_counter()-start:.2f}s')

        print(f'> caching results')
        start = perf_counter()
        for bam_idx in range(len(args.bams)):
            np.save(f'{args.out}{bam_idx}', bam_scores[bam_idx])
        print(f'    runtime: {perf_counter()-start:.2f}s')

    else:
        print(f'> loading results')
        start = perf_counter()
        for bam_idx in range(len(args.bams)):
            bam_scores.append(np.load(f'{args.out}{bam_idx}.npy'))
        print(f'    runtime: {perf_counter()-start:.2f}s')

    print('> plotting purity')
    start = perf_counter()
    plot_purity(bam_scores, args.out)
    print(f'    runtime: {perf_counter()-start:.2f}s')


if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()

    try:
        main()
    except KeyboardInterrupt:
        print("\nERROR: Program terminated.")
        exit(1)
