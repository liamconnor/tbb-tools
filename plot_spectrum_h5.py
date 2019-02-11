import sys

import numpy as np
import matplotlib.pylab as plt

import read_tbb_data

if __name__=='__main__':
    fn = sys.argv[1]

    T = read_tbb_data.TBBh5_Reader(fn)
    station_name = list(set(T.get_stations_groups()[1]))
    data, mapping, tt = T.station_data(station_name[0])
    print(len(tt))
    print(len(tt[0]))

    data_I = T.voltage_to_intensity(data)
    fig = plt.figure()
    plt.plot(data_I.mean(0).mean(-1))
    plt.semilogy()
    plt.xlabel('subband number')
    plt.ylabel('intensity')

    fig2 = plt.figure()
    plt.imshow(np.log10(data_I.mean(0)), aspect='auto', cmap='Greys')
    plt.ylabel('subband number')
    plt.xlabel('time')
    plt.show()
    
#    times = T.get_time_arr(fn)
#    plt.plot(t[-1]+4150*10*(freq**-2-freq[-1]**-2), linewidth=3.5, color='orange')
#    plt.xlabel('Freq [MHz]', fontsize=15)
#    plt.ylabel('Arrival time [unix]', fontsize=15);plt.legend(['Theory', 'TBB_Writer Metadata'])
#    plt.show()

