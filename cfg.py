import multiprocessing as mp
from collections import defaultdict

# globals (for multi-processing without passing args)
args = None

# counter (for in-progress printing)
pos_count = mp.Value('i', 0)
read_count = mp.Value('i', 0)
results_count = mp.Value('i', 0)

# global enum for bases ACGT
bases = ['N', 'A', 'C', 'G', 'T']
nbases = len(bases)
base_dict = defaultdict(int)
base_dict['N'] = 0
base_dict['A'] = 1
base_dict['C'] = 2
base_dict['G'] = 3
base_dict['T'] = 4
base_dict['n'] = 0
base_dict['a'] = 1
base_dict['c'] = 2
base_dict['g'] = 3
base_dict['t'] = 4
base_dict['-'] = 5

# convert Cigar to string
cigars = 'MIDNSHP=XB'
cigar_dict = {
    'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4,
    'H': 5, 'P': 6, '=': 7, 'X': 8, 'B': 9
}

__version__ = "0.0.1"
