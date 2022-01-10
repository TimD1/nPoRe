import pysam
from collections import defaultdict
import os, subprocess

from cig import *
from util import *
import cfg


def filter_overlaps(in_vcf_fn, out_vcf_fn):
    ''' 
    Filter overlapping variants in VCF

    3   ATTTTTTT   A        # kept
    5   T          C        # removed
    6   TTTT       T        # removed
    '''
    in_vcf = pysam.VariantFile(in_vcf_fn, 'r')
    out_vcf = pysam.VariantFile(out_vcf_fn, 'w', header=in_vcf.header)

    prev_contig = ""
    prev_stop = 0
    for record in in_vcf.fetch():
        if record.contig != prev_contig:
            prev_stop = 0
            prev_contig = record.contig
        if record.start < prev_stop: # exclusive
            continue
        else:
            out_vcf.write(record)
            prev_stop = record.stop
    out_vcf.close()



def split_vcf(vcf_fn, vcf_out_pre='', filter_unphased=False):
    '''
    Splits phased VCF into hapVCFs.
    '''

    # create VCFs
    vcf = pysam.VariantFile(vcf_fn, 'r')
    vcf_out1_fn = vcf_out_pre + "1.vcf.gz"
    vcf_out2_fn = vcf_out_pre + "2.vcf.gz"
    vcf_out1 = pysam.VariantFile(vcf_out1_fn, 'w', header=vcf.header)
    vcf_out2 = pysam.VariantFile(vcf_out2_fn, 'w', header=vcf.header)

    # read diploid VCF, copying records into hap1 or hap2 VCF
    unphased = True
    records = False
    for ctg, start, stop in cfg.args.regions:
        for record in vcf.fetch(ctg, start, stop):
            records = True

            # only deal with 1-sample VCFs for now
            for sample in record.samples:
                gt = record.samples[sample]['GT']

            # TODO: with complex vars, can shorten alleles
            if len(record.alleles) == 3:  # different variants

                if record.alleles[gt[0]] != '*': # skip spanning deletion
                    record1 = record.copy()
                    for sample in record1.samples:
                        record1.samples[sample]['GT'] = ()
                    alleles1 = [record.alleles[0], record.alleles[gt[0]]]
                    record1.alleles = alleles1
                    vcf_out1.write(record1)

                if record.alleles[gt[1]] != '*':
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
                if filter_unphased:
                    phased_snp = False
                    for sample in record1.samples:
                        if 'PS' in record1.samples[sample]:
                            phased_snp = True
                    if not phased_snp: continue

                for sample in record1.samples:
                    record1.samples[sample]['GT'] = ()
                vcf_out1.write(record1)

            elif gt[1]:                    # hap2 variant only
                record2 = record.copy()
                if filter_unphased:
                    phased_snp = False
                    for sample in record2.samples:
                        if 'PS' in record2.samples[sample]:
                            phased_snp = True
                    if not phased_snp: continue

                for sample in record2.samples:
                    record2.samples[sample]['GT'] = ()
                vcf_out2.write(record2)

            elif not gt[0] and not gt[1] and record.alleles[0] == record.alleles[1]:
                pass # ignore, same as ref

            else: # PySAM error, consider it homozygous variant
                # TODO: file bug in PySam repo?
                record1 = record.copy()
                for sample in record1.samples:
                    record1.samples[sample]['GT'] = ()
                vcf_out1.write(record1)
                vcf_out2.write(record1)

            if gt[0] and not gt[1]: # 1|0 means phased, otherwise all are 0|1
                unphased = False

    if not records:
        print(f"\nWARNING: VCF file has no variants in selected region.")
    elif unphased:
        print(f"\nWARNING: VCF file may be unphased.")

    # close and index haploid VCFs
    vcf_out1.close()
    vcf_out2.close()
    subprocess.run(['tabix', '-f', '-p', 'vcf', vcf_out1_fn])
    subprocess.run(['tabix', '-f', '-p', 'vcf', vcf_out2_fn])

    return vcf_out1_fn, vcf_out2_fn



def merge_vcfs(vcf_fn1, vcf_fn2, out_fn=None):

    # create iterators for both input VCFs
    vcf1 = pysam.VariantFile(vcf_fn1, 'r')
    vcf2 = pysam.VariantFile(vcf_fn2, 'r')

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

    for contig, start, end in cfg.args.regions:

        # get iterators
        vcf1_itr = iter(vcf1.fetch(contig, start, end))
        vcf2_itr = iter(vcf2.fetch(contig, start, end))
        record1 = next(vcf1_itr, None)
        record2 = next(vcf2_itr, None)

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



def apply_vcf(vcf_fn, hap):

    vcf = pysam.VariantFile(vcf_fn, 'r')
    data = []
    for contig, start, stop in cfg.args.regions:
        cig = ''
        seq = ''
        ref_ptr = 0
        ref = get_fasta(cfg.args.ref, contig)
        len_ref = len(ref)
        for record in vcf.fetch(contig, start, stop):
            pos = record.pos - 1
            if not record.qual or record.qual < cfg.args.min_qual:
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

        # add remaining (all matches)
        cig += '=' * (len_ref - ref_ptr)
        seq += ref[ref_ptr:]
        data.append( (contig, hap, seq, ref, cig) )

    return data



def gen_vcf(hap_data, hap, vcf_out_pre = ''):


    # create VCF header
    vcf_header = pysam.VariantHeader()
    vcf_header.add_sample('SAMPLE')
    for ctg, start, stop in cfg.args.regions:
        vcf_header.add_meta('contig', items=[('ID', ctg)])
    vcf_header.add_meta('FORMAT', 
            items=[('ID',"GT"), ('Number',1), ('Type','String'), 
                ('Description','Genotype')])
    vcf_header.add_meta('FORMAT', 
            items=[('ID',"GQ"), ('Number',1), ('Type','Integer'), 
                ('Description','Genotype quality score')])

    # create VCF file
    vcf_out_fn = f"{vcf_out_pre}{hap}.vcf.gz"
    vcf_out = pysam.VariantFile(vcf_out_fn, 'w', header=vcf_header)

    # for each contig, convert CIGAR to VCF
    contig_lengths = ""
    for contig, hap, seq, ref, cigar in hap_data:
        contig_lengths += f"##contig=<ID={contig},length={len(ref)}>\n"
        ref_ptr = 0
        seq_ptr = 0
        cig_ptr = 0
        cig_len = len(cigar)
        typ = 0
        while cig_ptr < cig_len:

            if cigar[cig_ptr] == '=':
                typ = 0
                ref_ptr += 1
                seq_ptr += 1
                cig_ptr += 1

            elif cigar[cig_ptr] == 'X':
                typ = 0
                record = vcf_out.header.new_record(
                        contig = contig,
                        start = ref_ptr,
                        alleles = (ref[ref_ptr], seq[seq_ptr]),
                        qual = 60,
                        filter = 'PASS'
                )
                vcf_out.write(record)
                ref_ptr += 1
                seq_ptr += 1
                cig_ptr += 1

            elif cigar[cig_ptr] == 'M':
                typ = 0
                # match, don't add to VCF
                if ref[ref_ptr] == seq[seq_ptr]:
                    ref_ptr += 1
                    seq_ptr += 1
                    cig_ptr += 1
                else: # simple sub
                    record = vcf_out.header.new_record(
                            contig = contig,
                            start = ref_ptr,
                            alleles = (ref[ref_ptr], seq[seq_ptr]),
                            qual = 60,
                            filter = 'PASS'
                    )
                    vcf_out.write(record)
                    ref_ptr += 1
                    seq_ptr += 1
                    cig_ptr += 1

            elif cigar[cig_ptr] == 'D':
                typ = 1
                del_len = 0
                while cig_ptr < cig_len and cigar[cig_ptr] == 'D':
                    del_len += 1
                    cig_ptr += 1
                if ref_ptr > 0:
                    record = vcf_out.header.new_record(
                            contig = contig,
                            start = ref_ptr-1,
                            alleles = (ref[ref_ptr-1:ref_ptr+del_len], ref[ref_ptr-1]),
                            qual = 60,
                            filter = 'PASS'
                    )
                    vcf_out.write(record)
                ref_ptr += del_len

            elif cigar[cig_ptr] == 'I':
                typ = 2
                ins_len = 0
                while cig_ptr < cig_len and cigar[cig_ptr] == 'I':
                    ins_len += 1
                    cig_ptr += 1
                if ref_ptr > 0 and seq_ptr > 0:
                    record = vcf_out.header.new_record(
                            contig = contig,
                            start = ref_ptr-1,
                            alleles = (ref[ref_ptr-1], seq[seq_ptr-1:seq_ptr+ins_len]),
                            qual = 60,
                            filter = 'PASS'
                    )
                    vcf_out.write(record)
                seq_ptr += ins_len

            else:
                print(f"\nERROR: unrecognized CIGAR operation '{cigar[cig_ptr]}'")
                exit(1)
    vcf_out.close()
    # except IndexError:
    #     print(f'INDEX ERROR: ctg: {contig}, hap: {hap}, ref_ptr: {ref_ptr}, len(ref): {len(ref)}, seq_ptr: {seq_ptr}, len(seq): {len(seq)}, cig_ptr: {cig_ptr}, cig_len: {cig_len}, ref_len: {ref_len(cigar)}, seq_len: {seq_len(cigar)}')
    #     exit(1)
    # except ValueError:
    #     print(f'VALUE ERROR: ctg: {contig}, hap: {hap}, ref_ptr: {ref_ptr}, len(ref): {len(ref)}, seq_ptr: {seq_ptr}, len(seq): {len(seq)}, cig_ptr: {cig_ptr}, cig_len: {cig_len}, ref_len: {ref_len(cigar)}, seq_len: {seq_len(cigar)}')
    #     if typ == 0:
    #         print(f'M: {ref[ref_ptr]}, {seq[seq_ptr]}')
    #     elif typ == 1:
    #         print(f'D: {ref[ref_ptr-1:ref_ptr+del_len]}, {ref[ref_ptr-1]}', del_len)
    #     elif typ == 2:
    #         print(f'I: {ref[ref_ptr-1]}, {seq[seq_ptr-1:seq_ptr+ins_len]}', ins_len)
    #     exit(1)

    # wait until file is closed
    while True:
        try:
            file_open = open(vcf_out_fn, 'a')
            if file_open: break
        except IOError:
            pass
    file_open.close()

    subprocess.run([ "gunzip", "-f", vcf_out_fn ])
    vcf_out_fn = vcf_out_fn[:-3]

    # fix header
    with open(vcf_out_fn, 'r+') as vcf_file:
        lines = [line for line in vcf_file.readlines() if line[:8] != "##contig"]
        lines.insert(2, contig_lengths)
        vcf_file.seek(0)
        vcf_file.writelines(lines)

    # wait until file is closed
    while True:
        try:
            file_open = open(vcf_out_fn, 'a')
            if file_open: break
        except IOError:
            pass
    file_open.close()

    subprocess.run([ "sed", "-i", "-e", "s/END=0/\./g", vcf_out_fn ])
    subprocess.run([ "bgzip", "-f", vcf_out_fn ])
    subprocess.run([ "tabix", "-f", "-p", "vcf", vcf_out_fn+".gz" ])

    return vcf_out_fn+".gz"
