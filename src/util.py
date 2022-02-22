from Bio import SeqIO
import pysam

import cfg


def get_fasta(reference, contig, start=None, end=None):
    return str(SeqIO.to_dict(SeqIO.parse(reference, "fasta"))[contig].seq[start:end])


def count_chunks(regions):
    return sum([(end-start+cfg.args.chunk_width-1)//cfg.args.chunk_width 
        for ctg,start,end in cfg.args.regions])


def get_bam_regions():

    # verify supplied BAM/FASTA exist
    try:
        ref = pysam.FastaFile(cfg.args.ref)
    except (AttributeError, IOError, ValueError):
        print(f"\nERROR: could not open --ref FASTA '{cfg.args.ref}'.")
        exit(1)
    try:
        bam = pysam.AlignmentFile(cfg.args.bam)
    except (AttributeError, IOError, ValueError):
        print(f"WARNING: could not open --bam BAM '{cfg.args.bam}'. ")

    # just align selected region
    if cfg.args.contig:
        if cfg.args.contig not in ref.references:
            print(f"ERROR: contig '{cfg.args.contig}' not present in '{cfg.args.ref}'. "
                    f"Valid contigs are: {ref.references}")
            exit(1)
        if cfg.args.contigs:
            print("\nERROR: can't set 'contig' and 'contigs'.")
            exit(1)
        if not cfg.args.contig_beg:
            cfg.args.contig_beg = 0
        if not cfg.args.contig_end:
            cfg.args.contig_end = \
                    ref.get_reference_length(cfg.args.contig)-1
        max_end = ref.get_reference_length(cfg.args.contig)-1
        cfg.args.regions = [(cfg.args.contig, 
                cfg.args.contig_beg, min(max_end, cfg.args.contig_end))]

    # add several contigs
    elif cfg.args.contigs:
        if cfg.args.contig_beg or cfg.args.contig_end:
            print("\nERROR: can't set start/endpoints with multiple contigs.")
            exit(1)
        cfg.args.regions = []
        for contig in cfg.args.contigs.split(','):
            if contig not in ref.references:
                print(f"ERROR: contig '{contig}' not present in '{cfg.args.ref}'. "
                    f"Valid contigs are: {ref.references}")
                exit(1)
            end = ref.get_reference_length(contig)-1
            cfg.args.regions.append((contig, 0, end))

    elif cfg.args.bed:
        try: 
            bed = open(cfg.args.bed, 'r')
            regions = [x.strip().split() for x in bed.readlines()]
            cfg.args.regions = [(ctg, int(start), int(stop)) for 
                    ctg, start, stop in regions]
        except FileNotFoundError:
            print(f"\nERROR: could not open 'cfg.args.bed' BED.")
            exit(1)

    # add all contigs
    else:

        # can't choose this option
        if cfg.args.contig_beg or cfg.args.contig_end:
            print("\nERROR: 'contig' not supplied, but start/endpoints set.")
            exit(1)

        # verify that this contig is included in ref FASTA
        if cfg.args.bam:
            cfg.args.regions = []
            for ctg, l in zip(bam.references, bam.lengths):
                if ctg not in ref.references:
                    print(f"WARNING: contig '{ctg}' present in '{cfg.args.bam}', but"
                          f" not '{cfg.args.ref}', skipping...")

                # verify that this contig has reads
                else:
                    if bam.count(ctg, 0, l-1):
                        cfg.args.regions.append((ctg, 0, l-1))
        else:
            for ctg, l in zip(ref.references, ref.lengths):
                cfg.args.regions.append((ctg, 0, l-1))



def get_vcf_regions():

    # verify supplied VCF/FASTA exist
    try:
        ref = pysam.FastaFile(cfg.args.ref)
    except (AttributeError, IOError, ValueError):
        print(f"\nERROR: could not open FASTA '{cfg.args.ref}'.")
        exit(1)
    try:
        vcf = pysam.VariantFile(cfg.args.vcf)
    except (AttributeError, IOError, ValueError):
        print(f"\nERROR: could not open VCF '{cfg.args.vcf}'.")
        exit(1)

    # just use selected region
    if cfg.args.contig:
        if cfg.args.contigs:
            print("\nERROR: can't set 'contig' and 'contigs'.")
            exit(1)
        if not cfg.args.contig_beg:
            cfg.args.contig_beg = 0
        if not cfg.args.contig_end:
            cfg.args.contig_end = \
                    ref.get_reference_length(cfg.args.contig)-1
        cfg.args.regions = [(cfg.args.contig, 
                cfg.args.contig_beg, cfg.args.contig_end)]

    # add several contigs
    elif cfg.args.contigs:
        if cfg.args.contig_beg or cfg.args.contig_end:
            print("\nERROR: can't set start/endpoints with multiple contigs.")
            exit(1)
        cfg.args.regions = []
        for contig in cfg.args.contigs.split(','):
            end = ref.get_reference_length(contig)-1
            cfg.args.regions.append((contig, 0, end))

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
                if sum(1 for _ in vcf.fetch(ctg, 0, l-1)):
                    cfg.args.regions.append((ctg, 0, l-1))
