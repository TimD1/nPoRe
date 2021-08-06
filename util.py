from Bio import SeqIO


def get_fasta(reference, contig):
    return str(SeqIO.to_dict(SeqIO.parse(reference, "fasta"))[contig].seq[:])
