import sys

import numpy as np
import matplotlib.pylab as plt

import read_tbb_data

if __name__=='__main__':
    fn = sys.argv[1]

    T = read_tbb_data.TBBh5_Reader(fn)
    station_name = list(set(T.get_stations_groups()[1]))
    data, mapping, tt = T.station_data(station_name[0])
#    tt = np.array(tt)
    data_I = T.voltage_to_intensity(data)

    data_tf = data_I.mean(0)
#    t0 = tt[0, :, 0]
    nfreq = data_tf.shape[0]
    nframe = data_tf.shape[1]

    for ii in range(nfreq):
        plt.plot(2**ii*data_tf[ii])

    plt.semilogy()
    plt.show()

#    index_diff = 1+(t0.max() - t0.min())//5.12e-6
#    data_tf_full = np.zeros([nfreq, nframe+index_diff])

#    for ii in range(0):
#        continue
#        offset = int((t0[ii]-t0.min())/5.12e-6)
#        print(ii, offset, data_tf.shape)
#        data_tf_full[ii, offset:offset+nframe] = data_tf[ii]

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

