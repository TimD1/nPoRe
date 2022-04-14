import os, sys, argparse, pysam
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(f'{parentdir}/src')

import cfg as cfg
from util import get_fasta
from aln import get_np_info, print_np_info
from cig import bases_to_int
from sankey import sankey

# variant types
SUB = 0
INS = 1
DEL = 2
CPX = 3
variants = {
        "substitution": SUB,
        "insertion": INS,
        "deletion": DEL,
        "complex": CPX,
        }
types = { 
        "t": "substitution",
        "i": "insertion",
        "d": "deletion",
        "c": "complex",
        }

# call types
TP = 0
FN = 1
FP = 2
calls = {"TP": TP, "FN": FN, "FP": FP}

# copy number variants
FALSE = 0
TRUE = 1
L = 0
L_IDX = 1

class VcfCounts:
    ''' Store VCF aggregate variant counts by type. '''

    def __init__(self):
        self.types = np.zeros((4,3), dtype=int)
        self.cnvs = np.zeros(2, dtype=int)

    def __str__(self):
        return (
            f'Overview\n'
            f'SUBs:     {self.types[SUB][TP]:7} TP\t{self.types[SUB][FN]:5} FN\t{self.types[SUB][FP]:5} FP\n'
            f'INSs:     {self.types[INS][TP]:7} TP\t{self.types[INS][FN]:5} FN\t{self.types[INS][FP]:5} FP\n'
            f'DELs:     {self.types[DEL][TP]:7} TP\t{self.types[DEL][FN]:5} FN\t{self.types[DEL][FP]:5} FP\n'
            f'COMPLEXs: {self.types[CPX][TP]:7} TP\t{self.types[CPX][FN]:5} FN\t{self.types[CPX][FP]:5} FP\n'
            f'CNVs:     {self.cnvs[TRUE]} INDELs are, {self.cnvs[FALSE]} INDELs are not.\n'
        )

    def add(self, variant, call):
        if call and call != ".":
            self.types[variants[variant], calls[call]] += 1



def count(vcf_fn):

    vcf = pysam.VariantFile(vcf_fn, 'r')
    data = VcfCounts()

    for record in vcf.fetch():

        # count call types
        ref_call = record.samples['TRUTH']['BD']
        query_call = record.samples['QUERY']['BD']
        ref_gt = record.samples['TRUTH']['GT']
        query_gt = record.samples['QUERY']['GT']
        ref_type = record.samples['TRUTH']['BI']
        query_type = record.samples['QUERY']['BI']


        if len(record.alleles) > 2 or type(ref_type) == tuple:
            if ref_type != ".": 
                if type(ref_type) == tuple or \
                        len(set([x for x in ref_gt if x])) > 1:
                    data.add("complex", ref_call)
                else:
                    data.add(types[ref_type[0]], ref_call)
            if query_type != "." and query_call != "TP":
                if type(query_type) == tuple or \
                        len(set([x for x in query_gt if x])) > 1:
                    data.add("complex", query_call)
                else:
                    data.add(types[query_type[0]], query_call)
        else:
            if ref_type != ".": 
                data.add(types[ref_type[0]], ref_call)
            if query_type != "." and query_call != "TP":
                data.add(types[query_type[0]], query_call)

        # get all actual (TP and FN) INDELS, test if CNV
        if ref_type != ".":

            # ignore complex variants
            if (len(record.alleles) <= 2 and type(ref_type) != tuple) or ( \
                    type(ref_type) == tuple and \
                    len(set([x for x in ref_gt if x])) == 1):

                if types[ref_type[0]] == "insertion":
                    ref = record.alleles[0]
                    alt = record.alleles[[x for x in ref_gt if x != 0][0]] # first non-ref alleles
                    pos = record.pos-1 + len(ref)
                    ins = alt[len(ref):]
                    refseq = cfg.args.refs[record.contig][pos:pos+20]
                    np_info = get_np_info(bases_to_int(refseq))
                    np_info_seq = get_np_info(bases_to_int(ins+refseq))
                    cnv = False
                    for n in range(1, cfg.args.max_n+1):
                        n_idx = n-1
                        if np_info[0, L, n_idx] and np_info_seq[0, L, n_idx]:
                            if len(ins) % n == 0 and ins[:n] == refseq[:n]:
                                data.cnvs[TRUE] += 1
                                cnv = True
                                break
                    if not cnv:
                        data.cnvs[FALSE] += 1
                    # print(record.alleles, refseq, ref_call, query_call, 
                    #         ref_gt, query_gt, ref_type, query_type, cnv)

                elif types[ref_type[0]] == "deletion":
                    ref = record.alleles[0]
                    alt = record.alleles[[x for x in ref_gt if x != 0][0]] # first non-ref allele
                    pos = record.pos-1 + len(alt)
                    dell = ref[len(alt):]
                    refseq = cfg.args.refs[record.contig][pos:pos+20]
                    np_info = get_np_info(bases_to_int(refseq))
                    cnv = False
                    for n in range(1, cfg.args.max_n+1):
                        n_idx = n-1
                        if np_info[0, L, n_idx] and len(dell) % n == 0: # np start
                            data.cnvs[TRUE] += 1
                            cnv = True
                            break
                    if not cnv:
                        data.cnvs[FALSE] += 1
                    # print(record.alleles, refseq, ref_call, query_call, 
                    #         ref_gt, query_gt, ref_type, query_type, cnv)
    return data



def disc_pie(data, suffix=""):
    fig, ax = plt.subplots()
    plt.pie(data.types[:,TP] + data.types[:,FN], 
            labels=variants.keys(), autopct='%1.1f%%', startangle=90)
    plt.suptitle(suffix)
    plt.tight_layout()
    plt.savefig(f"img/disc_pie{'_' if suffix else ''}{suffix}.png", dpi=300)
    plt.close()



def error_pie(data, suffix=""):
    fig, ax = plt.subplots(2,2)
    for x in range(2):
        for y in range(2):
            i = x*2+y
            ax[x,y].pie(data.types[i,:], labels=calls.keys(), colors = 'gry',
                    autopct='%1.1f%%', startangle=90)
            ax[x,y].set_title(list(variants.keys())[i])
    plt.suptitle(suffix)
    plt.tight_layout()
    plt.savefig(f"img/call_pie{'_' if suffix else ''}{suffix}.png", dpi=300)
    plt.close()



def plot_sankey(all_data, np_data=None, np_sizes=None):

    # sankey 0
    label0 = ["Other"] + [f"{i}-Polymer" for i in range(1, cfg.args.max_n+1)]
    total_size = sum([np_sizes[f"np_{i}"] for i in range(cfg.args.max_n+1)])
    leftweight0 = [np_sizes[f"np_{i}"]/total_size for i in range(cfg.args.max_n+1)]
    total_errors = sum([np.sum(np_data[i].types[:,1:]) \
            for i in range(cfg.args.max_n+1)])
    rightweight0 = [np.sum(np_data[i].types[:,1:])/total_errors \
            for i in range(cfg.args.max_n+1)]

    # sankey 1
    left1 = (["Other"] + [f"{i}-Polymer" for i in range(1, cfg.args.max_n+1)])* 4
    right1 = ["Substitutions"]*7 + ["Insertions"]*7 + ["Deletions"]*7 + ["Complex"]*7
    total = np.sum([x.types[:,1:] for x in np_data])
    weight1 = []
    for error_type in range(4):
        for np_idx in range(cfg.args.max_n+1):
            weight1.append(sum(np_data[np_idx].types[error_type,1:])/total)

    bottom_labels = [
            "Relative\nRegion Sizes",
            "Errors\nby Region",
            "Variant Call\nError Types",
            ]
    sankey(
            lefts = [label0, left1],
            rights = [label0, right1],
            colors = cfg.args.colors,
            leftWeights = [leftweight0, weight1],
            rightWeights = [rightweight0, weight1],
            rightColors = [False, True],
            gaps = [False, False],
            bottoms = bottom_labels,
            fontsize = 12,
            figureName = "img/sankey1",
    )


    # # sankey 2
    # left2 = ["True Positive"]*(cfg.args.max_n+1) + \
    #         ["False Negative"]*(cfg.args.max_n+1) + \
    #         ["False Positive"]*(cfg.args.max_n+1)
    # right2 = (["Other"] + [f"{i}-Polymer" for i in range(1, cfg.args.max_n+1)])* 3
    # total = np.sum([x.types[:,:] for x in np_data])
    # leftweight2 = [np.sum(np_data[i].types[:,TP]) / total \
    #                 for i in range(cfg.args.max_n+1)] + \
    #         [np.sum(np_data[i].types[:,FN]) / total \
    #                 for i in range(cfg.args.max_n+1)] + \
    #         [np.sum(np_data[i].types[:,FP]) / total \
    #                 for i in range(cfg.args.max_n+1)]

    # total = np.sum([x.types[:,1:] for x in np_data])
    # rightweight2 = [0]*(cfg.args.max_n+1) + [np.sum(np_data[i].types[:,FN]) / total \
    #                 for i in range(cfg.args.max_n+1)] + \
    #         [np.sum(np_data[i].types[:,FP]) / total \
    #                 for i in range(cfg.args.max_n+1)]

    # # sankey 3
    # label3 = ["Other"] + [f"{i}-Polymer" for i in range(1, cfg.args.max_n+1)]
    # total_errors = sum([np.sum(np_data[i].types[:,1:]) \
    #         for i in range(cfg.args.max_n+1)])
    # leftweight3 = [np.sum(np_data[i].types[:,1:])/total_errors \
    #         for i in range(cfg.args.max_n+1)]
    # total_size = sum([np_sizes[f"np_{i}"] for i in range(cfg.args.max_n+1)])
    # rightweight3 = [np_sizes[f"np_{i}"]/total_size for i in range(cfg.args.max_n+1)]


################################################################################

    # sankey 4
    label4 = ["Other"] + [f"{i}-Polymer" for i in range(1, cfg.args.max_n+1)]
    total_size = sum([np_sizes[f"np_{i}"] for i in range(cfg.args.max_n+1)])
    leftweight4 = [np_sizes[f"np_{i}"]/total_size for i in range(cfg.args.max_n+1)]
    total = sum([np.sum(np_data[i].types[INS,:2]) + np.sum(np_data[i].types[DEL,:2]) \
            for i in range(cfg.args.max_n+1)])
    rightweight4 = [( np.sum(np_data[i].types[INS,:2]) + \
            np.sum(np_data[i].types[DEL,:2]) ) / total \
            for i in range(cfg.args.max_n+1)]

    # sankey 5
    left5 = ["Other"]*2 + [f"{n}-Polymer" for i in 
            range(1, cfg.args.max_n+1) for n in [i,i]]
    right5 = ["General INDEL", "Copy Number INDEL"]*(cfg.args.max_n+1)
    total = np.sum( [x.cnvs[:] for x in np_data] )
    weight5 = [w / total for i in range(cfg.args.max_n+1) for w in np_data[i].cnvs]

    bottom_labels = [
            "Relative\nRegion Sizes",
            "True INDELs\nby Region",
            "True INDEL\nVariant Types"
            ]
    sankey(
            lefts = [label4, left5],
            rights = [label4, right5],
            colors = cfg.args.colors,
            leftWeights = [leftweight4, weight5],
            rightWeights = [rightweight4, weight5],
            rightColors = [False, True],
            gaps = [False, False],
            bottoms = bottom_labels,
            fontsize = 12,
            figureName = "img/sankey2",
    )



def get_region_sizes(beds):

    # initialize dict
    sizes = {}
    sizes['all'] = 0
    sizes['np_all'] = 0
    for i in range(cfg.args.max_n+1):
        sizes[f'np_{i}'] = 0

    CONTIG, START, END = 0, 1, 2
    for name in sizes.keys():
        bedfilename = beds.replace("$", name)
        with open(bedfilename, 'r') as bedfile:
            for line in bedfile:
                region = line.strip().split()
                sizes[name] += int(region[END]) - int(region[START])

    return sizes



def main():

    print(f'> extracting reference contigs')
    cfg.args.refs = {}
    ref = pysam.FastaFile(cfg.args.ref)
    # for ctg in ref.references:
    for ctg in ["chr20", "chr21", "chr22"]:
    # for ctg in ["chr20"]:
        cfg.args.refs[ctg] = get_fasta(cfg.args.ref, ctg)

    print("> calculating 'all' stats")
    all_data = count(cfg.args.vcfs.replace("$", "all"))

    print("> plotting 'all'")
    disc_pie(all_data)
    error_pie(all_data)

    print("ALL")
    print(all_data)

    # print("> calculating BED sizes")
    # sizes = get_region_sizes(cfg.args.beds)
    # for name, size in sizes.items():
    #     print(f'{name}: {size}')

    sizes = {
            "all": 2875001522,
            "np_all": 1050014093,
            "np_0": 1824987429,
            "np_1": 960760575,
            "np_2": 83289186,
            "np_3": 16369108,
            "np_4": 9977189,
            "np_5": 3393943,
            "np_6": 1162724,
    }

    with mp.Pool() as pool:
        np_data = pool.map(count, [cfg.args.vcfs.replace("$", f"np_{i}") 
                for i in range(cfg.args.max_n+1)])

    for i in range(cfg.args.max_n+1):
        disc_pie(np_data[i], suffix=f'np{i}')
        error_pie(np_data[i], suffix=f'np{i}')
        print(f"NP {i}")
        print(np_data[i])

    plot_sankey(all_data, np_data, sizes)



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("--ref", 
            default= "/x/gm24385/reference/"
            "GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta")
    parser.add_argument("--vcfs", 
            default="/home/timdunn/NanoporeScripts/results/data/"
            "par-happy/g5-$-truth41std-eval41.vcf")
    parser.add_argument("--beds",
            default="/x/gm24385/reference/$.bed")
    parser.add_argument("--max_n", type=int, default=6)
    parser.add_argument("--max_l", type=int, default=100)
    return parser



def set_colors():
    cfg.args.colors = {
        'Substitutions':'#f78c1b',
        'General INDEL':'#1b7ef7',
        'Copy Number INDEL':'#9912c9',
        'Insertions':'#9bd937',
        'Deletions':'#f71b1b',
        'Complex':'#ffff66',
        'False Negative':'#f71b1b',
        'True Positive':'#12e23f',
        'False Positive':'#f78c1b',
    }

    # grayscale for n-polymers
    for n in range(cfg.args.max_n+1):
        chars = '0123456789ABCDEF'
        cfg.args.colors[f'{n}-Polymer' if n else "Other"] = f'#{chars[12-2*n]*6}'



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    set_colors()
    main()
