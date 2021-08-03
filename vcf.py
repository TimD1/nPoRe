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
