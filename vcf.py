import pysam
from collections import defaultdict
import os
from Bio import SeqIO

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

        # TODO: with complex vars, can shorten alleles
        if len(record.alleles) == 3:  # different variants
            record1 = record.copy()
            for sample in record1.samples:
                record1.samples[sample]['GT'] = ()
            alleles1 = [record.alleles[0], record.alleles[gt[0]]]
            record1.alleles = alleles1
            vcf_out1.write(record1)

            record2 = record.copy()
            for sample in record2.samples:
                record2.samples[sample]['GT'] = ()
            alleles2 = [record.alleles[0], record.alleles[gt[1]]]
            record2.alleles = alleles2
            vcf_out2.write(record2)
            continue

        if gt[0]:                    # hap1 variant
            record1 = record.copy()
            for sample in record1.samples:
                record1.samples[sample]['GT'] = ()
            vcf_out1.write(record1)

        if gt[1]:                    # hap2 variant
            record2 = record.copy()
            for sample in record2.samples:
                record2.samples[sample]['GT'] = ()
            vcf_out2.write(record2)

        if gt[0] and not gt[1]: # 1|0 means phased, otherwise all are 0|1
            unphased = False

    if unphased:
        print(f"WARNING: VCF file '{cfg.args.vcf}' may be unphased.")

    return vcf_out1_fn, vcf_out2_fn



def merge_vcfs(vcf_fn1, vcf_fn2, out_fn=None):

    # create iterators for both input VCFs
    vcf1 = pysam.VariantFile(vcf_fn1, 'r')
    vcf2 = pysam.VariantFile(vcf_fn2, 'r')
    vcf1_itr = iter(vcf1.fetch(cfg.args.contig, 
                cfg.args.contig_beg, cfg.args.contig_end))
    vcf2_itr = iter(vcf2.fetch(cfg.args.contig, 
                cfg.args.contig_beg, cfg.args.contig_end))
    record1 = next(vcf1_itr, None)
    record2 = next(vcf2_itr, None)

    # create output VCF
    if out_fn is None:
        prefix1 = vcf_fn1.split(os.extsep)[0]
        prefix2 = vcf_fn2.split(os.extsep)[0]
        if prefix1[-1] == '1' and prefix2[-1] == '2' and \
                prefix1[:-1] == prefix2[:-1]:
            out_fn = prefix1[:-1] + '.vcf.gz'
        else:
            out_fn = 'out.vcf.gz'
    vcf_out = pysam.VariantFile(out_fn, 'w', header=vcf1.header)

    # step through both input VCFs
    while record1 or record2:
        pos1 = float('inf') if not record1 else record1.pos
        pos2 = float('inf') if not record2 else record2.pos
        pos = min(pos1, pos2)
        hap1 = pos1 == pos
        hap2 = pos2 == pos

        # merge variants, tracking genotypes
        if hap1 and hap2:
            if record1.alleles == record2.alleles: # same 1|1
                record = record1.copy()
                for sample in record.samples:
                    record.samples[sample]['GT'] = (1,1)
                vcf_out.write(record)

            else:                                        # diff 1|2
                record = record1.copy()
                reflen_diff = len(record1.alleles[0]) - len(record2.alleles[0])
                if reflen_diff > 0: # hap1 longer del
                    record.alleles = (
                            record1.alleles[0], 
                            record1.alleles[1], 
                            record2.alleles[1] + record1.alleles[0] \
                                [len(record2.alleles[0]):len(record1.alleles[0])])
                else:               # same, hap2 longer del
                    reflen_diff = abs(reflen_diff)
                    record.alleles = (
                            record1.alleles[0], 
                            record1.alleles[1] + record2.alleles[0]\
                                [len(record1.alleles[0]):len(record2.alleles[0])], 
                            record2.alleles[1])

                for sample in record.samples:
                    record.samples[sample]['GT'] = (1,2)
                vcf_out.write(record)

        elif hap1:                                       # 1|0
            record = record1.copy()
            for sample in record.samples:
                record.samples[sample]['GT'] = (1,0)
            vcf_out.write(record)

        elif hap2:                                       # 0|1
            record = record2.copy()
            for sample in record.samples:
                record.samples[sample]['GT'] = (0,1)
            vcf_out.write(record)

        # continue
        if hap1: record1 = next(vcf1_itr, None)
        if hap2: record2 = next(vcf2_itr, None)



def get_fasta(reference, contig):
    return str(SeqIO.to_dict(SeqIO.parse(reference, "fasta"))[contig].seq[:])



def apply_vcf(vcf_fn, ref_fn):

    ref = get_fasta(ref_fn, cfg.args.contig)
    vcf = pysam.VariantFile(vcf_fn, 'r')

    cig = ''
    seq = ''
    last_pos = 0
    for record in vcf.fetch(
            cfg.args.contig, cfg.args.contig_beg, cfg.args.contig_end):
        pos = record.pos - 1

        # all positions not in VCF are unchanged
        seq += ref[last_pos:pos]
        cig += '=' * (pos - last_pos)
        last_pos = pos

        # compare current position for sub/ins/del
        seq += record.alleles[1]
        if record.alleles[0][0] == record.alleles[1][0]:
            cig += '='
            last_pos += 1
        else:
            cig += 'X'
            last_pos += 1

        indel_len = len(record.alleles[1]) - len(record.alleles[0])
        if indel_len > 0:   # insertion
            cig += 'I' * indel_len

        elif indel_len < 0:   # deletion
            indel_len = abs(indel_len)
            cig += 'D' * indel_len
            last_pos += indel_len

    return seq, cig






def gen_vcf():
    pass

