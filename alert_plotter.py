import os

import numpy as np
import matplotlib.pylab as plt 
import argparse 
import h5py

startsec, endsec=0, 5
dm = 26.8 

def read_h5(fn, startsec, endsec, tint=1, fint=1):
    file = h5py.File(fn, 'r')
    timeres=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_0/'].attrs['INCREMENT']
    freqaxis=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_1/'].attrs['AXIS_VALUES_WORLD']
    startsample=(int(startsec/timeres)/tint)*tint
    endsample=(int(endsec/timeres)/tint)*tint
    data=file["/SUB_ARRAY_POINTING_000/BEAM_000/STOKES_0"][startsample:endsample,:]
    ntime,nfreq=data.shape
    nsamples=endsample-startsample
    data3=np.sum(data.reshape(nsamples/tint,tint,nfreq),axis=1)
    data3=np.sum(data3.reshape(nsamples/tint,nfreq/fint,fint),axis=2)
    timeres *= tint 
    freqaxis = freqaxis[fint//2::fint]
    time_arr = np.linspace(startsec,endsec,data3.shape[0])
    data3 = data3.T

    return data3, timeres, time_arr, freqaxis

def rebin_tf(data, tint=1, fint=1):
    """ Rebin in time and frequency accounting 
    for zeros
    """
    nfreq, ntime = data.shape

    if fint>1:
        # Rebin in frequency
        data_ = data[:nfreq//fint*fint].reshape(nfreq//fint, fint, ntime)
        weights = (data_.mean(-1)>0).sum(1)
        data_ = np.sum(data_, axis=1) / weights[:, None]
        data_[np.isnan(data_)] = 0.0
    else:
        data_ = data

    if tint>1:
        # Rebin in time
        data_ = data_[:, :ntime//tint*tint].reshape(nfreq//fint, ntime//tint, tint)
        weights = (data_.mean(0)>0).sum(1)
        data_ = np.sum(data_, axis=-1) / weights[None]
        data_[np.isnan(data_)] = 0.0

    return data_

def plot_im(data, freq=(109863281.25, 187976074.21875), time_arr=None,
            taxis=1, vmax=3, vmin=-2):
    data_ = data.copy()
    if taxis==0:
        print("Transposing for plotter")
        data_ = data_.T

    nfreq, ntime = data.shape
    print("Plotter: Assuming nfreq: %d and ntime: %d" % (nfreq, ntime))
    data_ -= np.mean(data_,axis=-1,keepdims=True)
    data_ /= np.std(data_, axis=-1, keepdims=True)

    if time_arr is None:
        extent = [0, ntime, freq[0], freq[-1]]
        xlab_ = 'Time (sample)'
    else:
        extent = [time_arr[0],time_arr[1],freq[0],freq[-1]]
        xlab_ = 'Time (sec)'

    data_[data_==0] = np.inf
    fig = plt.figure()
    plt.imshow(data_,origin='lower',aspect='auto',
               vmin=vmin,vmax=vmax,
               extent=extent, cmap='hot')
    plt.xlabel(xlab_, fontsize=16)
    plt.ylabel('Freq', fontsize=16)
    plt.show()

def plot_dedisp(data_dd, time_arr=None, dm=0):
    """ Visualize dedispersed data
    """
    nfreq,ntime = data_dd.shape
    if time_arr is None:
        time_arr = np.linspace(0, ntime, ntime)
        xlab_ = 'Time (sample)'
    else:
        xlab_ = 'Time (sec)'

    dd_ts = data_dd.mean(0)
    fig = plt.figure()
    print(time_arr.shape, dd_ts.shape)
    plt.plot(time_arr, dd_ts)
    plt.xlabel(xlab_, fontsize=16)
    plt.legend(['Freq-avg time series DM=%0.2f'%dm])
    plt.show()

def dedisperse(data, dm, timeres=4.9152e-4, 
               freq=(109863281.25, 187976074.21875), 
               freq_ref=None):
    data = data.copy()
    
    nfreq, ntime = data.shape[0], data.shape[1]

    if len(freq)==2:
        freq_arr = np.linspace(freq[0], freq[-1], nfreq)
    else:
        freq_arr = freq

    # Convert from Hz to MHz
    freq_arr *= 1e-6

    if freq_ref is None:
        freq_ref = np.median(freq_arr)

    tdelay = 4.148e3*dm*(freq_arr**-2 - freq_ref**-2)
    ntime = len(data[0])

    maxind_arr = []

    for ii, f in enumerate(freq_arr):
        data[ii] = np.roll(data[ii], -np.int(tdelay[ii]/timeres))

    return data

def dumb_clean(data, taxis=1, plot_clean=False):
    """ Do a simple RFI cleaning based on 
    variance in time and frequency direction.
    """
    nfreq, ntime = data.shape
    print("RFI Cleaner: Assuming nfreq: %d and ntime: %d" % (nfreq, ntime))
    stat = np.std(data,axis=taxis)
    mask = np.where(stat>2*np.median(stat))[0]
    ind_all = range(0, nfreq)
    ind_use = np.delete(ind_all, mask)

    if plot_clean:
        fig = plt.figure()
        plt.plot(ind_use, data[ind_use].mean(1), '.')
        plt.plot(mask, data[mask].mean(1), '.', color='r')
        plt.legend(['Unmasked', 'Masked %0.1f percent' % (100*len(mask)/float(nfreq))])
        plt.xlabel('Freq channel', fontsize=16)
        plt.semilogy()
        plt.show()

    data[mask] = 0

    return data, ind_use, mask

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="RFI clean, dedisperse, plot beamformed LOFAR data",
                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fn', 
                        help='filename (assumes h5 for now)', 
                        type=str, required=True)
    parser.add_argument('-rfi', '--rfi', 
                        help='Clean RFI', action='store_true')
    parser.add_argument('-tint', '--tint', help='downsample in time', 
                        default=1, type=int)
    parser.add_argument('-fint', '--fint', help='downsample in frequency', 
                        default=1, type=int)
    parser.add_argument('-dm', '--dm', help='dispersion measure', 
                        default=0, type=float)
    parser.add_argument('-T', '--times', help='start and end time in seconds', nargs=2,
                        type=float, default=(0,5))
    parser.add_argument('-p', '--plot_all', 
                        help='Make plots along the way', action='store_true')

    inputs = parser.parse_args()
    startsec, endsec = inputs.times[0], inputs.times[1]

    data, timeres, time_arr, freqaxis = read_h5(inputs.fn, startsec, endsec, tint=1, fint=1)

    if inputs.rfi:
        data, ind_use, mask = dumb_clean(data, plot_clean=inputs.plot_all)
    if inputs.dm>0:
        data = dedisperse(data, dm, freq=freqaxis, timeres=timeres)
    if inputs.fint>1 or inputs.tint>1:
        data = rebin_tf(data, tint=inputs.tint, fint=inputs.fint)
        time_arr = time_arr[inputs.tint//2::inputs.tint]
        freqaxis = freqaxis[inputs.fint//2::inputs.fint]
    if inputs.plot_all:
        print(time_arr.shape, freqaxis.shape, data.shape)
        plot_im(data, time_arr, vmax=3, vmin=-2)
        plot_dedisp(data, time_arr, dm=inputs.dm)








