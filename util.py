from Bio import SeqIO
import pysam

import cfg


def get_fasta(reference, contig):
    return str(SeqIO.to_dict(SeqIO.parse(reference, "fasta"))[contig].seq[:])



def get_bam_regions():

    # verify supplied BAM/FASTA exist
    try:
        ref = pysam.FastaFile(cfg.args.ref)
    except (AttributeError, IOError, ValueError):
        print("\nERROR: could not open FASTA '{cfg.args.ref}'.")
        exit(1)
    try:
        bam = pysam.AlignmentFile(cfg.args.bam)
    except (AttributeError, IOError, ValueError):
        print("\nERROR: could not open BAM '{cfg.args.bam}'.")
        exit(1)

    # just align selected region
    if cfg.args.contig:
        if not cfg.args.contig_beg:
            cfg.args.contig_beg = 0
        if not cfg.args.contig_end:
            cfg.args.contig_end = \
                    ref.get_reference_length(cfg.args.contig)-1
        cfg.args.regions = [(cfg.args.contig, 
                cfg.args.contig_beg-1, cfg.args.contig_end-1)]

    # add all contigs
    else:

        # can't choose this option
        if cfg.args.contig_beg or cfg.args.contig_end:
            print("\nERROR: 'contig' not supplied, but start/endpoints set.")
            exit(1)

        # verify that this contig is included in ref FASTA
        cfg.args.regions = []
        for ctg, l in zip(bam.references, bam.lengths):
            if ctg not in ref.references:
                print(f"WARNING: contig '{ctg}' present in '{cfg.args.bam}', but"
                      f" not '{cfg.args.ref}', skipping...")

            # verify that this contig has reads
            else:
                if bam.count(ctg, 0, l-1):
                    cfg.args.regions.append((ctg, 0, l-1))



def get_vcf_regions():

    # verify supplied VCF/FASTA exist
    try:
        ref = pysam.FastaFile(cfg.args.ref)
    except (AttributeError, IOError, ValueError):
        print("\nERROR: could not open FASTA '{cfg.args.ref}'.")
        exit(1)
    try:
        bcf = pysam.VariantFile(cfg.args.vcf)
    except (AttributeError, IOError, ValueError):
        print("\nERROR: could not open VCF '{cfg.args.vcf}'.")
        exit(1)

    # just use selected region
    if cfg.args.contig:
        if not cfg.args.contig_beg:
            cfg.args.contig_beg = 0
        if not cfg.args.contig_end:
            cfg.args.contig_end = \
                    ref.get_reference_length(cfg.args.contig)-1
        cfg.args.regions = [(cfg.args.contig, 
                cfg.args.contig_beg-1, cfg.args.contig_end-1)]

    # use all contigs
    else:

        # invalid option
        if cfg.args.contig_beg or cfg.args.contig_end:
            print("\nERROR: 'contig' not supplied, but start/endpoints set.")
            exit(1)

        cfg.args.regions = []
        for ctg in vcf.header.contigs:

            # check contig is in reference FASTA
            if ctg not in ref.references:
                print(f"WARNING: contig '{ctg}' present in '{cfg.args.vcf}', but"
                      f" not '{cfg.args.ref}', skipping...")
            else:
                l = ref.get_reference_length(ctg)

                 # check this contig has mutations
                 if sum(1 for _ in vcf.fetch(ctg, 0, l-1):
                    cfg.args.regions.append((ctg, 0, l-1))
