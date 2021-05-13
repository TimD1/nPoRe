import pysam

def get_positions(vcf_filename: str, min_qual: int, window: int):
    '''
    Parse VCF file and extract candidate variant positions.
     - VCF filename must be zipped (supply <filename>.vcf.gz)
     - VCF index must be present (<filename>.vcf.gz.tbi)
    '''

    # get candidate variants from VCF
    vcf = pysam.VariantFile(vcf_filename, 'r')
    all_pos = [record.start for record in vcf.fetch() if record.qual > min_qual]

    # filter nearby variant candidates (within window)
    last_pos = 0
    pos = []
    for p in all_pos:
        if p > last_pos + 2*window:
            pos.append(p)
            last_pos = p
    return pos
