import multiprocessing as mp

# globals (for multi-processing without passing args)
args = None

# counter (for in-progress printing)
pos_count = mp.Value('i', 0)
