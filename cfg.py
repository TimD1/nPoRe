import multiprocessing as mp

# globals (for multi-processing without passing args)
args = None

# counter (for in-progress printing)
pos_count = mp.Value('i', 0)

# global enum for bases ACGT
bases = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

