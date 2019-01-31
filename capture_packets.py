import sys

import numpy as np 
import socket, struct, time

import read_tbb_data

try:
    import matplotlib.pylab as plt
except:
    print("Couldn't import matplotlib")

#server_socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
#server_socket.bind(("", 31664))

def compare_dat_raw(fndat, fnraw):
    T = read_tbb_data.TBB_Rawdata('new')

    header_full, data_full, crc32_full = T.read_all_data(fndat, nframe=10000000)
    data_dat = np.array(data_full)
    data_dat = data_dat.flatten()

    data_raw = np.fromfile(fnraw, np.int16)

    print(".raw file length %d" % len(data_raw))
    print(".dat file length %d" % len(data_dat))

    n = min(len(data_raw), len(data_dat))

    dat_diff = np.abs(data_raw[:n] - data_dat[:n])

    for ii, vv in enumerate(dat_diff):
        if vv != 0:
            print("First nonzero frame: %d" % ii)
            break

    return data_dat, data_raw, dat_diff, data_full, header_full

def listen(addresses, fn=None):
    T = read_tbb_data.TBB_Rawdata('new')

    server_sockets = []

    for ad in addresses:
        print("Binding to %d" % ad)
        server_socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        server_socket.bind(("", ad))
        server_sockets.append(server_socket)

    while True:
        if fn is not None:
            f = open(fn, 'a+')

        for server_socket in server_sockets:
            recvdata, address = server_socket.recvfrom(2400)
            h = T.parse_header(recvdata[:88])
            print(address)
            print(T.print_one_frame(h))
            if fn is not None:
                f.write(recvdata)

if __name__=='__main__':
    fn = sys.argv[1]
    print(fn)
    addresses = []

    for address in sys.argv[2:]:
        addresses.append(int(address))

    listen(addresses, fn=None)
