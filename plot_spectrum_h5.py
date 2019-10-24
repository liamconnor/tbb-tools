import sys

import numpy as np
import matplotlib.pylab as plt

import read_tbb_data

if __name__=='__main__':
    fn = sys.argv[1]

    ax_dict = {0:'ndipole', 1:'subband', 2:'time'}

    try:
        axis = int(sys.argv[2])
    except:
        axis = 0

    T = read_tbb_data.TBBh5_Reader(fn)
    station_name = list(set(T.get_stations_groups()[1]))
    data, mapping, tt = T.station_data(station_name[0], nsubband_tot=None)
    #data_I = T.voltage_to_intensity(data)

    print("Data array shape is (ndipole, nsubband, nsample)=%s" % str((data.shape)))
    print("Averaging over %s axis" % ax_dict[axis])

    data = data.mean(axis)

    if axis==0:
        ylab = 'SUBBAND'
        xlab = 'TIME [rebinned samples]'
    if axis==1:
        ylab = 'DIPOLE INDEX'
        xlab = 'TIME [rebinned samples]'        
    if axis==2:
        ylab = 'DIPOLE INDEX'
        xlab = 'SUBBAND'     

    fig2 = plt.figure()
    plt.imshow(np.log10(data), aspect='auto', cmap='Greys')
    plt.ylabel(ylab, fontsize=20)
    plt.xlabel(xlab, fontsize=20)
    plt.show()
    
#    times = T.get_time_arr(fn)
#    plt.plot(t[-1]+4150*10*(freq**-2-freq[-1]**-2), linewidth=3.5, color='orange')
#    plt.xlabel('Freq [MHz]', fontsize=15)
#    plt.ylabel('Arrival time [unix]', fontsize=15);plt.legend(['Theory', 'TBB_Writer Metadata'])
#    plt.show()

