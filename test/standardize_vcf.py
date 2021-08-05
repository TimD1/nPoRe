import os, sys, subprocess
import pickle
import argparse
import numpy as np

currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

from vcf import *
from cig import *
from bam import hap_to_bam
import cfg as cfg

def fix_vcf(vcf):
    while True:
        try:
            file_open = open(vcf, 'a')
            if file_open: break
        except IOError:
            pass
    file_open.close()

    if vcf[-3:] == ".gz":
        subprocess.run([ "gunzip", "-f", vcf ])
        vcf = vcf[:-3]
    subprocess.run([ "sed", "-i", "-e", "s/END=0/\./g", vcf ])
    subprocess.run([ "bgzip", "-f", vcf ])
    subprocess.run([ "tabix", "-p", "vcf", vcf+".gz" ])


def main():

    print(f"> splitting vcf '{cfg.args.vcf}'")
    vcf1, vcf2 = split_vcf(cfg.args.vcf, "out/hap")

    print(f"> indexing '{vcf1}' and '{vcf2}'")
    subprocess.run(['tabix', '-p', 'vcf', vcf1])
    subprocess.run(['tabix', '-p', 'vcf', vcf2])

    print(f"> applying '{vcf1}'")
    seq1, cig1 = apply_vcf(vcf1, cfg.args.ref)
    print(f"> applying '{vcf2}'")
    seq2, cig2 = apply_vcf(vcf2, cfg.args.ref)

    print(f"> reading reference")
    ref = get_fasta(cfg.args.ref, cfg.args.contig)

    print(f"> standardizing cigars")
    cigar1_data = standardize_cigar(("1", "chr19", 0, cig1, ref, seq1))
    cigar2_data = standardize_cigar(("2", "chr19", 0, cig2, ref, seq2))

    print(f"\n> saving CIGAR data")
    pickle.dump(cigar1_data, open('out/cigar1_data.pkl', 'wb'))
    pickle.dump(cigar2_data, open('out/cigar2_data.pkl', 'wb'))

    print(f"> loading CIGAR data")
    cigar1_data = pickle.load(open('out/cigar1_data.pkl', 'rb'))
    cigar2_data = pickle.load(open('out/cigar2_data.pkl', 'rb'))

    print(f"> aligned haps to BAM")
    hap_to_bam(cigar1_data, "out/hap")
    hap_to_bam(cigar2_data, "out/hap")

    print('> generating VCF1')
    prefix = f"{args.vcf.split(os.extsep)[0]}_std"
    vcf1 = gen_vcf(cigar1_data, prefix)
    print('> generating VCF2')
    vcf2 = gen_vcf(cigar2_data, prefix)

    print(f"> indexing '{vcf1}' and '{vcf2}'")
    fix_vcf(vcf1)
    fix_vcf(vcf2)

    print(f"> merging vcfs: '{vcf1}' and '{vcf2}'")
    merge_vcfs(vcf1, vcf2)



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )
    parser.add_argument("--vcf", type=str, 
            default='/x/gm24385/chr19/ref/HG002_GRCh38_1_22_v4.1_draft_benchmark.vcf.gz')
    parser.add_argument("--ref", type=str, 
            default='/x/gm24385/chr19/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta')
    parser.add_argument("--contig", type=str, default="chr19")
    parser.add_argument("--contig_beg", type=int, default=1)
    parser.add_argument("--contig_end", type=int, default=40_000_000)
    parser.add_argument("--indels_only", default=True, action="store_true")
    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
