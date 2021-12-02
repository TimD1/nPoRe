import os, sys, subprocess
import pickle
import argparse
import numpy as np

from vcf import *
from cig import *
from bam import *
import cfg as cfg

def main():

    print(f"> indexing input VCFs")
    for vcf in cfg.args.vcfs:
        if vcf[-4:] == ".vcf":
            subprocess.run(["bgzip", "-f", vcf])
            vcf = vcf + ".gz"
        elif vcf[-7:] == ".vcf.gz":
            pass
        else:
            print("\nERROR: input VCFs must be in VCF format.")
            exit(1)
        subprocess.run(['tabix', '-f', '-p', 'vcf', vcf])

    print(f"> fixing VCFs")
    for in_vcf, out_vcf in zip(cfg.args.vcfs, cfg.args.out_vcfs):
        print(f"    {in_vcf}")
        fix_phasing(in_vcf, out_vcf, cfg.args.plot)

    print(f"> indexing output VCFs")
    for vcf in cfg.args.out_vcfs:
        if vcf[-4:] == ".vcf":
            subprocess.run(["bgzip", "-f", vcf])
            vcf = vcf + ".gz"
        elif vcf[-7:] == ".vcf.gz":
            pass
        else:
            print("\nERROR: output VCFs must be in VCF format.")
            exit(1)
        subprocess.run(['tabix', '-f', '-p', 'vcf', vcf])


def argparser():
    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter,
            add_help = False
    )

    parser.add_argument("--bam", type=str)
    parser.add_argument("--vcfs", action='store', type=str, nargs='+')
    parser.add_argument("--out_vcfs", action='store', type=str, nargs='+')
    parser.add_argument("--plot", action='store_true')

    return parser



if __name__ == "__main__":

    parser = argparser()
    cfg.args = parser.parse_args()

    if len(cfg.args.vcfs) != len(cfg.args.out_vcfs):
        print("\nERROR: must have same number of input/output VCFs.")
        exit(1)
    if cfg.args.vcfs is None or cfg.args.vcfs == []:
        print("\nERROR: must supply at least one input VCF.")
        exit(1)

    main()
