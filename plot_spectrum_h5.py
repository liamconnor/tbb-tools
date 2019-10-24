import sys

import numpy as np
import matplotlib.pylab as plt

import read_tbb_data

if __name__=='__main__':
    fn = sys.argv[1]

    ax_dict = {0:'ndipole', 1:'subband', 2:'time'}

    T = read_tbb_data.TBBh5_Reader(fn)
    station_name = list(set(T.get_stations_groups()[1]))
    data, mapping, tt = T.station_data(station_name[0], nsubband_tot=10)
    #data_I = T.voltage_to_intensity(data)

    try:
        axis = sys.argv[2]
    except:
        axis = 0

    try:
        axis = int(axis)
    except:
        pass 

    print(axis)

    if type(axis)==int:
        data = data.mean(axis)
        print("Data array shape is (ndipole, nsubband, nsample)=%s" % str((data.shape)))
        print("Averaging over %s axis" % ax_dict[axis])

    if axis==0:
        ylab = 'SUBBAND'
        xlab = 'TIME [rebinned samples]'
    if axis==1:
        ylab = 'DIPOLE INDEX'
        xlab = 'TIME [rebinned samples]'        
    if axis==2:
        ylab = 'DIPOLE INDEX'
        xlab = 'SUBBAND'     

    if axis is not 'all':
        fig2 = plt.figure()
        plt.imshow(np.log10(data), aspect='auto', cmap='Greys')
        plt.ylabel(ylab, fontsize=20)
        plt.xlabel(xlab, fontsize=20)
        plt.show()
    else:
        fig2 = plt.figure()
        plt.add_subplot(131)
        plt.imshow(np.log10(data.mean(0)), aspect='auto', cmap='Greys')
        plt.ylabel(ylab, fontsize=20)
        plt.xlabel(xlab, fontsize=20)
        plt.add_subplot(132)
        plt.imshow(np.log10(data.mean(1)), aspect='auto', cmap='Greys')
        plt.ylabel(ylab, fontsize=20)
        plt.xlabel(xlab, fontsize=20)   
        plt.add_subplot(133)
        plt.imshow(np.log10(data.mean(2)), aspect='auto', cmap='Greys')
        plt.ylabel(ylab, fontsize=20)
        plt.xlabel(xlab, fontsize=20)
        plt.show()
#    times = T.get_time_arr(fn)
#    plt.plot(t[-1]+4150*10*(freq**-2-freq[-1]**-2), linewidth=3.5, color='orange')
#    plt.xlabel('Freq [MHz]', fontsize=15)
#    plt.ylabel('Arrival time [unix]', fontsize=15);plt.legend(['Theory', 'TBB_Writer Metadata'])
#    plt.show()

