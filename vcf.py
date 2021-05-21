import pysam

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
