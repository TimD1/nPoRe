import os, sys, argparse, pysam
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

import cfg as cfg

# variant types
SUB = 0
INS = 1
DEL = 2
CPX = 3
variants = {
        "substitution": SUB,
        "insertion": INS,
        "deletion": DEL,
        "complex": CPX}
types = { 
        "t": "substitution",
        "i": "insertion",
        "d": "deletion"}

# call types
TP = 0
FN = 1
FP = 2
calls = {"TP": TP, "FN": FN, "FP": FP}

class VcfCounts:
    ''' Store VCF aggregate variant counts by type. '''

    def __init__(self):
        self.matrix = np.zeros((4,3), dtype=int)

    def __str__(self):
        return (
            f'Overview\n'
            f'SUBs:     {self.matrix[SUB][TP]:7} TP\t{self.matrix[SUB][FN]:5} FN\t{self.matrix[SUB][FP]:5} FP\n'
            f'INSs:     {self.matrix[INS][TP]:7} TP\t{self.matrix[INS][FN]:5} FN\t{self.matrix[INS][FP]:5} FP\n'
            f'DELs:     {self.matrix[DEL][TP]:7} TP\t{self.matrix[DEL][FN]:5} FN\t{self.matrix[DEL][FP]:5} FP\n'
            f'COMPLEXs: {self.matrix[CPX][TP]:7} TP\t{self.matrix[CPX][FN]:5} FN\t{self.matrix[CPX][FP]:5} FP\n'
        )

    def add(self, variant, call):
        # print(variant, call)
        if call and call != ".":
            self.matrix[variants[variant], calls[call]] += 1



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

        if len(record.alleles) > 2:

            if happy:
                if ref_type != "." and ref_call != "TP": 
                    if type(ref_type) == tuple:
                        data.add("complex", ref_call)
                        # print("complex", ref_call)
                    else:
                        data.add(types[ref_type[0]], ref_call)
                        # print(types[ref_type[0]], ref_call)
                if query_type != ".":
                    if type(query_type) == tuple:
                        data.add("complex", query_call)
                        # print("complex", query_call)
                    else:
                        data.add(types[query_type[0]], query_call)
                        # print(types[query_type[0]], query_call)

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

    return data



def disc_pie(data):
    fig, ax = plt.subplots()
    plt.pie(data.matrix[:,TP] + data.matrix[:,FN], 
            labels=variants.keys(), autopct='%1.1f%%', startangle=90)
    plt.tight_layout()
    plt.savefig("img/disc_pie.png", dpi=300)
    plt.close()



def error_pie(data):
    fig, ax = plt.subplots(2,2)
    for x in range(2):
        for y in range(2):
            i = x*2+y
            ax[x,y].pie(data.matrix[i,:], labels=calls.keys(), colors = 'gry',
                    autopct='%1.1f%%', startangle=90)
            ax[x,y].set_title(list(variants.keys())[i])
    plt.tight_layout()
    plt.savefig("img/call_pie.png", dpi=300)
    plt.close()



def sankey(data):

    labels = [
            "Substitution", "Insertion", "Deletion", "Complex",   # call type
            "True Positive", "False Negative", "False Positive", # error category
            " ", "FN", "FP"
    ]

    source = [
            0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
            4, 5, 6
    ]

    target = [
            4, 5, 6, 4, 5, 6, 4, 5, 6, 4, 5, 6,
            7, 8, 9
    ]

    t = np.sum(data.matrix)
    value = [ 
            data.matrix[SUB,TP]/t, data.matrix[SUB,FN]/t, data.matrix[SUB,FP]/t, 
            data.matrix[INS,TP]/t, data.matrix[INS,FN]/t, data.matrix[INS,FP]/t, 
            data.matrix[DEL,TP]/t, data.matrix[DEL,FN]/t, data.matrix[DEL,FP]/t, 
            data.matrix[CPX,TP]/t, data.matrix[CPX,FN]/t, data.matrix[CPX,FP]/t,
            0.1, np.sum(data.matrix[:,FN])/np.sum(data.matrix[:,1:]),
            np.sum(data.matrix[:,FP])/np.sum(data.matrix[:,1:]),
    ]

    link = dict(source=source, target=target, value=value)
    node = dict(label=labels, pad=0, thickness=5)
    data = go.Sankey(link=link, node=node)
    fig = go.Figure(data)
    fig.write_image("img/sankey.png")



def main():
    happy_data = count(cfg.args.happy_vcf)
    disc_pie(happy_data)
    error_pie(happy_data)
    sankey(happy_data)
    
    print("HAPPY")
    print(happy_data)
    # print("TRUTH")
    # print(truth_data)



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("--happy_vcf", 
            default="/home/timdunn/NanoporeScripts/results/data/"
            "par-happy/g5-all-truth41-eval41.vcf")
    parser.add_argument("--truth_vcf", 
            default="/x/gm24385/test/guppy_5_0_6/g5-ra-hap/ref/1_std.vcf.gz")
    parser.add_argument("--max_n", type=int, default=6)
    parser.add_argument("--max_l", type=int, default=100)
    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
