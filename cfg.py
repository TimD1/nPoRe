import multiprocessing as mp

# globals (for multi-processing without passing args)
args = None

# counter (for in-progress printing)
pos_count = mp.Value('i', 0)
read_count = mp.Value('i', 0)
results_count = mp.Value('i', 0)

# global enum for bases ACGT
bases = {'N': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4}

# convert Cigar to string
cigar = 'MIDNSHP=XB'

__version__ = "0.0.1"
