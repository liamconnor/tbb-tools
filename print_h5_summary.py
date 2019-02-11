import sys

import numpy as np

import read_tbb_data

if __name__=='__main__':
    fn = sys.argv[1]
    T = read_tbb_data.TBBh5_Reader(fn)
    T.print_summary_data(fn)
