import pysam
from collections import defaultdict
import os

from cig import *
from util import *
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
    print(f'\r    0 of {n} SNPs processed.', end='', flush=True)
    for i, snp in enumerate(snps):
        if snp.qual > cfg.args.min_qual:
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
        print(f'\r    {i+1} of {n} SNPs processed.', end='', flush=True)

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


def split_vcf(vcf_fn, vcf_out_pre=''):
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

        elif gt[0] and gt[1]:          # same hap1 and hap2 variant
            record1 = record.copy()
            for sample in record1.samples:
                record1.samples[sample]['GT'] = ()
            vcf_out1.write(record1)
            vcf_out2.write(record1)

        elif gt[0]:                    # hap1 variant only
            record1 = record.copy()
            for sample in record1.samples:
                record1.samples[sample]['GT'] = ()
            vcf_out1.write(record1)

        elif gt[1]:                    # hap2 variant only
            record2 = record.copy()
            for sample in record2.samples:
                record2.samples[sample]['GT'] = ()
            vcf_out2.write(record2)

        elif not gt[0] and not gt[1] and record.alleles[0] == record.alleles[1]:
            pass # ignore, same as ref

        else: # PySAM error, consider it homozygous variant
            # print(f"ERROR: no variant: {record.pos} {record.alleles} {gt}")
            # TODO: file bug in PySam repo?
            record1 = record.copy()
            for sample in record1.samples:
                record1.samples[sample]['GT'] = ()
            vcf_out1.write(record1)
            vcf_out2.write(record1)

        if gt[0] and not gt[1]: # 1|0 means phased, otherwise all are 0|1
            unphased = False

    if unphased:
        print(f"WARNING: VCF file may be unphased.")

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
        prefix1 = '.'.join(vcf_fn1.split(os.extsep)[:-2])
        prefix2 = '.'.join(vcf_fn2.split(os.extsep)[:-2])
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
                indel_len1 = len(record1.alleles[1]) - len(record1.alleles[0])
                indel_len2 = len(record2.alleles[1]) - len(record2.alleles[0])

                if indel_len1 >= 0 and indel_len2 >= 0: # both sub/ins
                    record.alleles = (
                            record1.alleles[0], 
                            record1.alleles[1], 
                            record2.alleles[1])

                if indel_len1 < 0 and indel_len2 < 0: # both del
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
                                record2.alleles[0], 
                                record1.alleles[1] + record2.alleles[0]\
                                    [len(record1.alleles[0]):len(record2.alleles[0])], 
                                record2.alleles[1])

                elif indel_len1 < 0 and indel_len2 >= 0:
                    record.alleles = (
                            record1.alleles[0], 
                            record1.alleles[1], 
                            record2.alleles[1] + record1.alleles[0][1:])

                elif indel_len2 < 0 and indel_len1 >= 0:
                    record.alleles = (
                            record2.alleles[0], 
                            record2.alleles[1], 
                            record1.alleles[1] + record2.alleles[0][1:])

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



def apply_vcf(vcf_fn, ref):

    len_ref = len(ref)

    cig = ''
    seq = ''
    ref_ptr = 0
    vcf = pysam.VariantFile(vcf_fn, 'r')
    for record in vcf.fetch(
            cfg.args.contig, cfg.args.contig_beg, cfg.args.contig_end):
        pos = record.pos - 1
        if record.qual < cfg.args.min_qual:
            continue

        if pos < ref_ptr: # overlapping indels, allow second if insertion
            indel_len = len(record.alleles[1]) - len(record.alleles[0])
            if indel_len > 0:
                seq += record.alleles[1][len(record.alleles[0]):]
                cig += 'I' * indel_len
            continue
        else:
            seq += ref[ref_ptr:pos]
            cig += '=' * (pos - ref_ptr)
            ref_ptr += pos-ref_ptr

        # compare current position for sub/ins/del
        seq += record.alleles[1]
        minlen = min(len(record.alleles[0]), len(record.alleles[1]))
        for i in range(minlen):
            if record.alleles[0][i] == record.alleles[1][i]:
                cig += '='
                ref_ptr += 1
            else:
                cig += 'X'
                ref_ptr += 1

        indel_len = len(record.alleles[1]) - len(record.alleles[0])
        if indel_len > 0:   # insertion
            cig += 'I' * indel_len

        elif indel_len < 0:   # deletion
            indel_len = abs(indel_len)
            cig += 'D' * indel_len
            ref_ptr += indel_len

        # print(f'pos: {pos}, len(seq): {len(seq)} seq_len(cig): {seq_len(cig)}, alleles: {record.alleles}')
        # if seq_len(cig) != len(seq) or ref_len(cig) != len_ref:
        #     input()

    # add remaining (all matches)
    cig += '=' * (len_ref - ref_ptr)
    seq += ref[ref_ptr:]

    return seq, cig



def gen_vcf(read_data, vcf_out_pre = ''):

    hap_id, ctg_name, start, stop, cigar, hap_cigar, ref, hap_ref, seq, hap = read_data

    # create VCF header
    vcf_header = pysam.VariantHeader()
    vcf_header.add_sample('SAMPLE')
    vcf_header.add_meta('contig', items=[('ID', ctg_name)])
    vcf_header.add_meta('FORMAT', 
            items=[('ID',"GT"), ('Number',1), ('Type','String'), 
                ('Description','Genotype')])
    vcf_header.add_meta('FORMAT', 
            items=[('ID',"GQ"), ('Number',1), ('Type','Integer'), 
                ('Description','Genotype quality score')])

    # create VCF file
    vcf_out_fn = f"{vcf_out_pre}{hap_id}.vcf.gz"
    vcf_out = pysam.VariantFile(vcf_out_fn, 'w', header=vcf_header)

    # convert CIGAR to VCF
    ref_ptr = start
    seq_ptr = 0
    cig_ptr = 0
    cig_len = len(cigar)
    while cig_ptr < cig_len:

        if cigar[cig_ptr] == '=':
            ref_ptr += 1
            seq_ptr += 1
            cig_ptr += 1

        elif cigar[cig_ptr] == 'X':
            record = vcf_out.header.new_record(
                    contig = ctg_name,
                    start = ref_ptr,
                    alleles = (ref[ref_ptr], seq[seq_ptr]),
                    qual = 10,
                    filter = 'PASS'
            )
            vcf_out.write(record)
            ref_ptr += 1
            seq_ptr += 1
            cig_ptr += 1

        elif cigar[cig_ptr] == 'M':
            # match, don't add to VCF
            if ref[ref_ptr] == seq[seq_ptr]:
                ref_ptr += 1
                seq_ptr += 1
                cig_ptr += 1
            else: # simple sub
                record = vcf_out.header.new_record(
                        contig = ctg_name,
                        start = ref_ptr,
                        alleles = (ref[ref_ptr], seq[seq_ptr]),
                        qual = 10,
                        filter = 'PASS'
                )
                vcf_out.write(record)
                ref_ptr += 1
                seq_ptr += 1
                cig_ptr += 1

        elif cigar[cig_ptr] == 'D':
            del_len = 0
            while cig_ptr < cig_len and cigar[cig_ptr] == 'D':
                del_len += 1
                cig_ptr += 1
            record = vcf_out.header.new_record(
                    contig = ctg_name,
                    start = ref_ptr-1,
                    alleles = (ref[ref_ptr-1:ref_ptr+del_len], ref[ref_ptr-1]),
                    qual = 10,
                    filter = 'PASS'
            )
            vcf_out.write(record)
            ref_ptr += del_len

        elif cigar[cig_ptr] == 'I':
            ins_len = 0
            while cig_ptr < cig_len and cigar[cig_ptr] == 'I':
                ins_len += 1
                cig_ptr += 1
            record = vcf_out.header.new_record(
                    contig = ctg_name,
                    start = ref_ptr-1,
                    alleles = (ref[ref_ptr-1], seq[seq_ptr-1:seq_ptr+ins_len]),
                    qual = 10,
                    filter = 'PASS'
            )
            vcf_out.write(record)
            seq_ptr += ins_len
        else:
            print(f"ERROR: unrecognized CIGAR operation '{cigar[cig_ptr]}'")
            exit(1)
    vcf_out.close()

    return vcf_out_fn
