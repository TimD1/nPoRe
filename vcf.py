import pysam
from collections import defaultdict

import cfg

def get_vcf_data():
    '''
    Parse VCF file into just sub information:
    { 'ctg1': [ [(posX, baseX)...],     < unknown haplotype (homozygous subs only)
                [(posY, baseY)...],     < hap1
                [(posZ, baseZ)...]]     < hap2
      'ctg2': ...
    }
    '''

    # count snps
    vcf_file = pysam.VariantFile(cfg.args.apply_vcf, 'r')
    snps = None
    if cfg.args.contig:
        snps = vcf_file.fetch(
                    cfg.args.contig, 
                    cfg.args.contig_beg, 
                    cfg.args.contig_end)
    else:
        snps = vcf_file.fetch()
    n = sum(1 for _ in snps)

    # get all snps in region of interest
    if cfg.args.contig:
        snps = vcf_file.fetch(
                    cfg.args.contig, 
                    cfg.args.contig_beg, 
                    cfg.args.contig_end)
    else:
        snps = vcf_file.fetch()

    vcf_dict = {}
    print(f'\r        0 of {n} SNPs processed.', end='', flush=True)
    for i, snp in enumerate(snps):
        if snp.qual > cfg.args.min_snp_qual:
            gt = None
            for sample in snp.samples:
                gts = snp.samples[sample]['GT']
                break

            if snp.contig not in vcf_dict:
                vcf_dict[snp.contig] = [ [], [], [] ]

            # homo vars
            if len(snp.alleles) == 2:
                if len(snp.alleles[0]) == 1 and len(snp.alleles[1]) == 1:
                    if gts[0]: # hap1
                        vcf_dict[snp.contig][1].append((snp.start, snp.alleles[1]))
                    if gts[1]: # hap2
                        vcf_dict[snp.contig][2].append((snp.start, snp.alleles[1]))
                    if gts[0] and gts[1]: # unknown hap
                        vcf_dict[snp.contig][0].append((snp.start, snp.alleles[1]))

            # hetero vars
            elif len(snp.alleles) == 3:
                if len(snp.alleles[1]) == len(snp.alleles[0]):
                    for hap_idx, gt in enumerate(gts):
                        if gt == 1:
                            vcf_dict[snp.contig][hap_idx].append((snp.start, snp.alleles[1][0]))
                if len(snp.alleles[2]) == len(snp.alleles[0]):
                    for hap_idx, gt in enumerate(gts):
                        if gt == 2:
                            vcf_dict[snp.contig][hap_idx].append((snp.start, snp.alleles[2][0]))
        print(f'\r        {i+1} of {n} SNPs processed.', end='', flush=True)

    # warn user if VCF is unphased. I've fallen for this twice now
    unphased = True
    for contig in vcf_dict:
        if vcf_dict[contig][0] != vcf_dict[contig][1]:
            unphased = False
    if unphased:
        print(f"WARNING: VCF file '{cfg.args.apply_vcf}' may be unphased.\n" + \
               "As a result, SNP splicing will not be correct."
                )

    return vcf_dict


def split_vcf(vcf_fn, vcf_out_pre=None):
    '''
    Splits phased VCF into hapVCFs.
    '''

    # open VCF
    vcf_file = pysam.VariantFile(vcf_fn, 'r')
    snps = None
    if cfg.args.contig:
        vcf = vcf_file.fetch(
                    cfg.args.contig, 
                    cfg.args.contig_beg, 
                    cfg.args.contig_end)
    else:
        vcf = vcf_file.fetch()

    # create output VCFs
    if vcf_out_pre is None:
        vcf_out_pre = vcf_fn + "_hap"
    vcf_out1_fn = vcf_out_pre + "1.vcf.gz"
    vcf_out2_fn = vcf_out_pre + "2.vcf.gz"
    vcf_out1 = pysam.VariantFile(vcf_out1_fn, 'w', header=vcf_file.header)
    vcf_out2 = pysam.VariantFile(vcf_out2_fn, 'w', header=vcf_file.header)

    # read diploid VCF, copying records into hap1 or hap2 VCF
    unphased = True
    for record in vcf:

        # only deal with 1-sample VCFs for now
        for sample in record.samples:
            gt = record.samples[sample]['GT']
            break
        assert(len(gt) == 2)

        if gt[0]: # hap1 variant
            vcf_out1.write(record.copy())
        if gt[1]: # hap2 variant
            vcf_out2.write(record.copy())
        if gt[0] and not gt[1]: # 1|0 means phased, otherwise all are 0|1
            unphased = False

    if unphased:
        print(f"WARNING: VCF file '{cfg.args.vcf}' may be unphased.")

    return vcf_out1_fn, vcf_out2_fn



def merge_vcf(vcf_fn1, vcf_fn2, out_fn=None):
    if out_fn is None:
        out_fn = vcf_fn1[:-5] # remove '_hapN'
    pass

def apply_vcf():
    pass

def gen_vcf():
    pass

