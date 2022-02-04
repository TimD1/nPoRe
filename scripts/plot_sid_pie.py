import os, sys, argparse, pysam
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import cfg as cfg
from util import get_fasta
from aln import get_np_info, print_np_info
from cig import bases_to_int
from mysankey import sankey

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
        # print(variant, call)
        if call and call != ".":
            self.types[variants[variant], calls[call]] += 1



def count(vcf_fn):

    vcf = pysam.VariantFile(vcf_fn, 'r')
    data = VcfCounts()

    for record in vcf.fetch():

        # happy truth VCF contains two samples
        happy = 'TRUTH' in list(record.samples) and \
                'QUERY' in list(record.samples)

        # count call types
        if happy:
            ref_call = record.samples['TRUTH']['BD']
            query_call = record.samples['QUERY']['BD']
        else:
            ref_call = query_call = 'TP'

        # count variant types
        ref_gt = record.samples['TRUTH']['GT']
        query_gt = record.samples['QUERY']['GT']
        ref_type = record.samples['TRUTH']['BI']
        query_type = record.samples['QUERY']['BI']
        # print(record.alleles, ref_gt, query_gt, ref_type, query_type)

        if len(record.alleles) > 2 or type(ref_type) == tuple:

            if happy:
                if ref_type != "." and ref_call != "TP": 
                    if type(ref_type) == tuple or \
                            len(set([x for x in ref_gt if x])) > 1:
                        data.add("complex", ref_call)
                    else:
                        data.add(types[ref_type[0]], ref_call)
                if query_type != ".":
                    if type(query_type) == tuple or \
                            len(set([x for x in query_gt if x])) > 1:
                        data.add("complex", query_call)
                    else:
                        data.add(types[query_type[0]], query_call)

            else:
                if ref_type != "." and ref_call != "TP": 
                    data.add("complex", ref_call)
                if query_type != ".":
                    data.add("complex", query_call)
        else:
            if ref_type != "." and ref_call != "TP": 
                data.add(types[ref_type[0]], ref_call)
            if query_type != ".":
                data.add(types[query_type[0]], query_call)

        # get all actual (TP and FN) INDELS, test if CNV
        if ref_type != ".":

            # ignore complex variants
            if type(ref_type) != tuple and \
                    len(set([x for x in ref_gt if x])) == 1 and \
                    types[ref_type[0]] == "insertion":
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

            if type(ref_type) != tuple and \
                    len(set([x for x in ref_gt if x])) == 1 and \
                    types[ref_type[0]] == "deletion":
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

    left = ["Substitutions"]*3 + ["Insertions"]*3 + ["Deletions"]*3 + ["Complex"]*3
    right = ["True Positive", "False Negative", "False Positive"] * 4
    sankey(
            left = left,
            right = right,
            leftWeight = all_data.types.flatten(),
            rightWeight = all_data.types.flatten(),
            colorDict = cfg.args.colors,
            rightColor = True,
            fontsize = 12
    )
    fig = plt.gcf()
    fig.set_size_inches(4, 6)
    plt.savefig("img/sankey0.png", bbox_inches="tight", dpi=300)

    left = ["True Positive", "False Negative", "False Positive"]
    sankey(
            left = left,
            right = left,
            leftWeight = np.sum(all_data.types, axis=0) / np.sum(all_data.types),
            rightWeight = [ 0, 
                    np.sum(all_data.types[:,FN])/np.sum(all_data.types[:,1:]),
                    np.sum(all_data.types[:,FP])/np.sum(all_data.types[:,1:])
                    ],
            colorDict = cfg.args.colors,
            fontsize = 12
    )
    fig = plt.gcf()
    fig.set_size_inches(4, 6)
    plt.savefig("img/sankey1.png", bbox_inches="tight", dpi=300)

    left = ["False Negative"]*(cfg.args.max_n+1) + \
            ["False Positive"]*(cfg.args.max_n+1)
    right = (["Other"] + [f"{i}-Polymer" for i in range(1, cfg.args.max_n+1)])* 2
    weights = [np.sum(np_data[i].types[:,FN]) for i in range(cfg.args.max_n+1)] + \
            [np.sum(np_data[i].types[:,FP]) for i in range(cfg.args.max_n+1)]
    sankey(
            left = left,
            right = right,
            leftWeight = weights,
            rightWeight = weights,
            aspect = 20,
            colorDict = cfg.args.colors,
            fontsize = 12
    )
    fig = plt.gcf()
    plt.savefig("img/sankey2.png", bbox_inches="tight", dpi=300)

    labels = ["Other"] + [f"{i}-Polymer" for i in range(1, cfg.args.max_n+1)]
    total_errors = sum([np.sum(np_data[i].types[:,1:]) \
            for i in range(cfg.args.max_n+1)])
    left_weights = [np.sum(np_data[i].types[:,1:])/total_errors \
            for i in range(cfg.args.max_n+1)]
    total_size = sum([np_sizes[f"np_{i}"] for i in range(cfg.args.max_n+1)])
    right_weights = [np_sizes[f"np_{i}"] for i in range(cfg.args.max_n+1)] / \
            total_size
    sankey(
            left = labels,
            right = labels,
            leftWeight = left_weights,
            rightWeight = right_weights,
            aspect = 20,
            colorDict = cfg.args.colors,
            fontsize = 12
    )
    fig = plt.gcf()
    plt.savefig("img/sankey3.png", bbox_inches="tight", dpi=300)

    labels = ["Other"] + [f"{i}-Polymer" for i in range(1, cfg.args.max_n+1)]
    total_size = sum([np_sizes[f"np_{i}"] for i in range(cfg.args.max_n+1)])
    left_weights = [np_sizes[f"np_{i}"] for i in range(cfg.args.max_n+1)] / total_size
    total = sum([np.sum(np_data[i].types[INS,:2]) + np.sum(np_data[i].types[DEL,:2]) \
            for i in range(cfg.args.max_n+1)])
    right_weights = [( np.sum(np_data[i].types[INS,:2]) + \
            np.sum(np_data[i].types[DEL,:2]) ) / total \
            for i in range(cfg.args.max_n+1)]
    sankey(
            left = labels,
            right = labels,
            leftWeight = left_weights,
            rightWeight = right_weights,
            aspect = 20,
            colorDict = cfg.args.colors,
            fontsize = 12
    )
    fig = plt.gcf()
    plt.savefig("img/sankey4.png", bbox_inches="tight", dpi=300)

    left_labels = ["Other"]*2 + [f"{n}-Polymer" for i in 
            range(1, cfg.args.max_n+1) for n in [i,i]]
    right_labels = ["General INDEL", "Copy Number Variant"]*(cfg.args.max_n+1)
    weights = [w for i in range(cfg.args.max_n+1) for w in np_data[i].cnvs]
    sankey(
            left = left_labels,
            right = right_labels,
            leftWeight = weights,
            rightWeight = weights,
            aspect = 20,
            colorDict = cfg.args.colors,
            fontsize = 12
    )
    fig = plt.gcf()
    plt.savefig("img/sankey5.png", bbox_inches="tight", dpi=300)



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
    for ctg in ref.references:
        cfg.args.refs[ctg] = get_fasta(cfg.args.ref, ctg)

    print("> calculating 'all' stats")
    all_data = count(cfg.args.vcfs.replace("$", "all"))

    print("> plotting 'all'")
    disc_pie(all_data)
    error_pie(all_data)

    print("ALL")
    print(all_data)

    print("> calculating BED sizes")
    sizes = get_region_sizes(cfg.args.beds)
    for name, size in sizes.items():
        print(f'{name}: {size}')

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
        'Substitutions':'#1b7ef7',
        'General Indel':'#1b7ef7',
        'Copy Number Variant':'#9BD937',
        'Insertions':'#9BD937',
        'Deletions':'#f71b1b',
        'Complex':'#9912c9',
        'False Negative':'#f71b1b',
        'True Positive':'#12e23f',
        'False Positive':'#f78c1b',
    }

    # grayscale for n-polymers
    for n in range(cfg.args.max_n+1):
        chars = '0123456789ABCDEF'
        cfg.args.colors[f'{n}-Polymer' if n else "Other"] = f'#{chars[2*n]*6}'



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    set_colors()
    main()
