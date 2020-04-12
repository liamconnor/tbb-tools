import sys 

import h5py

fnh5 = sys.argv[1]

f = h5py.File(fnh5,'r')
for pair in f["/SUB_ARRAY_POINTING_000/"]['BEAM_000'].attrs.items():
	print(pair)
