import os, sys, subprocess
import pickle
import argparse
import numpy as np

from vcf import *
from cig import *
from bam import *
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

    with open(vcf, 'r+') as vcf_file:
        lines = vcf_file.readlines()
        lines[2] = "##contig=<ID=chr19,length=58617616>\n"
        vcf_file.seek(0)
        vcf_file.writelines(lines)

    while True:
        try:
            file_open = open(vcf, 'a')
            if file_open: break
        except IOError:
            pass
    file_open.close()

    subprocess.run([ "sed", "-i", "-e", "s/END=0/\./g", vcf ])
    subprocess.run([ "bgzip", "-f", vcf ])
    subprocess.run([ "tabix", "-p", "vcf", vcf+".gz" ])


def main():

    # print(f"> splitting vcf")
    # vcf1, vcf2 = split_vcf(cfg.args.vcf, cfg.args.out)

    # print(f"> indexing vcfs")
    # subprocess.run(['tabix', '-p', 'vcf', vcf1])
    # subprocess.run(['tabix', '-p', 'vcf', vcf2])

    # print(f"> reading reference")
    # ref = get_fasta(cfg.args.ref, cfg.args.contig)

    # print(f"> converting vcfs and ref to sequences")
    # seq1, cig1 = apply_vcf(vcf1, ref)
    # seq2, cig2 = apply_vcf(vcf2, ref)

    # # package data
    # cigar1_data = ("1", "chr19", 0, 0, cig1, "="*len(ref), ref, ref, seq1, 1)
    # cigar2_data = ("2", "chr19", 0, 0, cig2, "="*len(ref), ref, ref, seq2, 2)

    # print("> calculating score matrices")
    # subs, hps = get_confusion_matrices()
    # cfg.args.sub_scores, cfg.args.hp_scores = calc_score_matrices(subs, hps)

    # print(f"> realigning hap sequences")
    # with mp.Pool() as pool:
    #     data = pool.map(realign_read, [cigar1_data, cigar2_data])

    # print(f"> saving CIGAR data")
    # cigar1_data = data[0]
    # cigar2_data = data[1]
    # pickle.dump(cigar1_data, open('cigar1_data.pkl', 'wb'))
    # pickle.dump(cigar2_data, open('cigar2_data.pkl', 'wb'))

    print(f"> loading CIGAR data")
    cigar1_data = pickle.load(open('cigar1_data.pkl', 'rb'))
    cigar2_data = pickle.load(open('cigar2_data.pkl', 'rb'))
    data = [cigar1_data, cigar2_data]

    if cfg.args.std_cigar:
        print(f"\n> chunking hap cigars")
        data = chunk_cigars(data)
        print(f"\n> standardizing hap cigars")
        with mp.Pool() as pool:
            data = pool.map(standardize_cigar, data)
        print(f"\n> stitching hap cigars")
        data = stitch_cigars(data)
    cigar1_data = data[0]
    cigar2_data = data[1]

    # print(f"> saving CIGAR data")
    # pickle.dump(cigar1_data, open('cigar1_data.pkl', 'wb'))
    # pickle.dump(cigar2_data, open('cigar2_data.pkl', 'wb'))

    # print(f"> loading CIGAR data")
    # cigar1_data = pickle.load(open('cigar1_data.pkl', 'rb'))
    # cigar2_data = pickle.load(open('cigar2_data.pkl', 'rb'))

    # print(f"> saving debug output to bam")
    # hap_to_bam(cigar1_data, cfg.args.out)
    # hap_to_bam(cigar2_data, cfg.args.out)
    # subprocess.run(['samtools', 'index', f'{cfg.args.out}1.bam'])
    # subprocess.run(['samtools', 'index', f'{cfg.args.out}2.bam'])

    print('\n> generating standardized vcfs')
    vcf1 = gen_vcf(cigar1_data, cfg.args.out+"_")
    vcf2 = gen_vcf(cigar2_data, cfg.args.out+"_")

    print(f"> indexing vcfs")
    fix_vcf(vcf1)
    fix_vcf(vcf2)

    print(f"> merging vcfs")
    out_fn = f"{cfg.args.out}.vcf.gz"
    merge_vcfs(vcf1, vcf2, out_fn)
    subprocess.run(['tabix', '-p', 'vcf', out_fn])



def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )

    parser.add_argument("vcf", type=str)
    parser.add_argument("ref", type=str)
    parser.add_argument("out", type=str)

    parser.add_argument("--contig", type=str)
    parser.add_argument("--contig_beg", type=int)
    parser.add_argument("--contig_end", type=int)

    parser.add_argument("--stats_dir", default="./stats")
    parser.add_argument("--recalc_cms", action="store_true")

    parser.add_argument("--max_hp", type=int, default=100)
    parser.add_argument("--chunk_width", type=int, default=10000)

    parser.add_argument("--std_cigar", action="store_true")
    parser.add_argument("--indels_only", action="store_true")
    return parser



if __name__ == "__main__":
    parser = argparser()
    cfg.args = parser.parse_args()
    main()
