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
    vcf_file = pysam.VariantFile(cfg.args.vcf, 'r')
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
    return vcf_dict
    


def get_positions(vcf_filename: str, min_qual: int, window: int):
    '''
    Parse VCF file and extract candidate variant positions.
     - VCF filename must be zipped (supply <filename>.vcf.gz)
     - VCF index must be present (<filename>.vcf.gz.tbi)
    '''

    # get candidate variants from VCF
    vcf = pysam.VariantFile(vcf_filename, 'r')
    all_pos = []
    for record in vcf.fetch():
        if record.qual > min_qual:
            gt = None
            for sample in record.samples:
                gt = record.samples[sample]['GT']
                break
            all_pos.append((record.start, record.alleles, gt))

    # filter nearby variant candidates (within window)
    # TODO: merge windows instead
    last_pos = 0
    pos = []
    for data in all_pos:
        if data[0] > last_pos + 2*window:
            pos.append(data)
            last_pos = data[0]
    return pos
