import os

import numpy as np
import matplotlib.pylab as plt 
import argparse 
import h5py
from datetime import datetime, timedelta
import time 

def read_h5(fn, time_range=(0,5)):
    """ Read in hdf5 beamformed data at native resolution.
    Transpose data to (nfreq, ntime)

    Parameters:
    ----------
    fn : str 
        path to .h5 file 
    startsec : float 
        number of second into file of start sample
    endsec : float 
        number of seconds into file of end sample

    Returns:
    -------
    data : ndarray
        (nfreq, ntime) intensity array
    timeres : float 
        time resolution of data 
    time_arr : ndarray 
        ntime length array of time samples w.r.t. start of file
    freqaxis : ndarray 
        frequency array in Hz
    """
    file = h5py.File(fn, 'r')
    start_time_file=datetime.strptime(file.attrs[u'OBSERVATION_START_UTC'][0:19],'%Y-%m-%dT%H:%M:%S')
    end_time_file=datetime.strptime(file.attrs[u'OBSERVATION_END_UTC'][0:19],'%Y-%m-%dT%H:%M:%S')
    start_time_file_unix=time.mktime(start_time_file.timetuple())
    end_time_file_unix=time.mktime(end_time_file.timetuple())
    print("Unix times in file: %d-%d" % (start_time_file_unix, end_time_file_unix))

    if len(time_range)==1:
        time_start = time_range[0]-0.25
        range_err = time_start>=start_time_file_unix and time_range[0]<=end_time_file_unix
        assert range_err, "That unix time is not in the file"
        print("Looking for data at unix time %s" % time_start)
        startsec=time_start-start_time_file_unix
        endsec=startsec+8
    elif len(time_range)==2:
        startsec, endsec = time_range

    if endsec<=startsec:
        print("Start time is larger than end time")
        exit()
    elif endsec-startsec>30.0:
        print("Are you sure you want %0.2f sec of data?" % (endsec-startsec))

    try:
        timeres=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_0/'].attrs['INCREMENT']
        freqaxis=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_1/'].attrs['AXIS_VALUES_WORLD']
        beamno=0
    except:
        timeres=file['SUB_ARRAY_POINTING_000/BEAM_001/COORDINATES/COORDINATE_0/'].attrs['INCREMENT']
        freqaxis=file['SUB_ARRAY_POINTING_000/BEAM_001/COORDINATES/COORDINATE_1/'].attrs['AXIS_VALUES_WORLD']
        beamno=1

    startsample=int(startsec/timeres)
    endsample=int(endsec/timeres)
    data=file["/SUB_ARRAY_POINTING_000/BEAM_00%d/STOKES_0" % beamno][startsample:endsample,:]
    if len(data)==0:
        print("No data in specified range")
        exit()
    ntime,nfreq=data.shape
    nsamples=endsample-startsample
    time_arr = np.linspace(startsec,endsec,ntime)

    return data.T, timeres, time_arr, freqaxis, start_time_file_unix

def read_cs_h5(fn): 
    """ Combine multiple files
    """
    pass    

def volt2intensity(data):
    assert data.shape[-1]%4==0, "Expecting 4*ntime (xr,xi,yr,yi)"
    I = data[:, ::4]**2 + data[:, 1::4]**2 \
        + data[:, 2::4]**2 + data[:, 3::4]**2
    return I

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

def read_npy(fn):
    data = np.load(fn)

    return data

def plot_im(data, freq=(109863281.25, 187976074.21875), time_arr=None,
            taxis=1, vmax=3, vmin=-2, figname=None):
    """ Plot the 2D time/freq data. Freq is in Hz. vmax and 
    vmin are in sigmas. 
    """
    data_ = data.copy()
    if taxis==0:
        print("Transposing for plotter")
        data_ = data_.T

    nfreq, ntime = data.shape
    print("Plotter: Assuming nfreq: %d and ntime: %d" % (nfreq, ntime))
    data_ -= np.mean(data_,axis=-1)[..., None]
    data_ /= np.std(data_, axis=-1)[..., None]

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
    plt.colorbar()
    if figname is not None:
        plt.savefig(figname)
    plt.show()

def plot_dedisp(data_dd, time_arr=None, dm=0, figname=None):
    """ Plot dedispersed data. data_dd should be a 
    (nfreq, ntime) array.
    """
    nfreq,ntime = data_dd.shape
    if time_arr is None:
        time_arr = np.linspace(0, ntime, ntime)
        xlab_ = 'Time (sample)'
    else:
        xlab_ = 'Time (sec)'

    dd_ts = data_dd.mean(0)
    fig = plt.figure()
    plt.plot(time_arr, dd_ts)
    plt.xlabel(xlab_, fontsize=16)
    plt.legend(['Freq-avg time series DM=%0.2f'%dm])
    if figname is not None:
        plt.savefig(figname)
    plt.show()

def dedisperse(data, dm, timeres=4.9152e-4, 
               freq=(109863281.25, 187976074.21875), 
               freq_ref=None):
    """ Dedisperse data to some dm

    Parameters:
    ----------
    data : ndarray
        (nfreq, ntime) array 
    dm : float 
        dispersion measure in pc cm**-3
    timeres : float 
        time resolution of data in seconds 
    freq : ndarray or tuple 
        should be the frequency array in Hz, 
        either length 2 or length nfreq 
    freq_ref : 
        if None, assumes center of band

    Returns: 
    -------
    Dedispersed data 
    """
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

def dumb_clean(data, taxis=1, plot_clean=False, figname=None):
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
        if figname is not None:
            plt.savefig(figname)
        plt.show()

    data[mask] = 0

    return data, ind_use, mask

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                     description="RFI clean, dedisperse, plot LOFAR data",
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
    parser.add_argument('-s', '--save_data', help='save data to .npy',
                        action='store_true')
    parser.add_argument('-T', '--times', 
                        help='start end times in sec or unix time', nargs='+',
                        type=float, default=(0,5.0))
    parser.add_argument('-p', '--plot_all', 
                        help='Make plots along the way', action='store_true')

    inputs = parser.parse_args()

    if len(inputs.times)==2:
        startsec, endsec = inputs.times[0], inputs.times[1]
    elif len(inputs.times)==1:
        startsec = inputs.times[0]
        endsec = startsec + 5.0

    if inputs.fn[-2:]=='h5':
        res = read_h5(inputs.fn, 
                      inputs.times)
        data, timeres, time_arr, freqaxis, start_time_file_unix = res
        ftype='.h5'
    elif inputs.fn[-3:]=='npy':
        data = read_npy(inputs.fn)
        ftype='.npy'
        if inputs.rfi:
            data, ind_use, mask = dumb_clean(data, 
                                        plot_clean=inputs.plot_all)
        if inputs.fint>1 or inputs.tint>1:
            data = rebin_tf(data, tint=inputs.tint, fint=inputs.fint)
        if inputs.plot_all:
            plot_im(data, vmax=1, vmin=-1)
            plot_dedisp(data, dm=inputs.dm)
        exit()

    # RFI clean data by zapping bad channels
    if inputs.rfi:
        data, ind_use, mask = dumb_clean(data, plot_clean=inputs.plot_all, 
                                         figname=inputs.fn.strip(ftype)+'_rfi.pdf')
    # Dedisperse data if given DM > 0
    if inputs.dm>0:
        data = dedisperse(data, inputs.dm, freq=freqaxis, 
                          timeres=timeres)
    # Downsample / channelize data
    if inputs.fint>1 or inputs.tint>1:
        data = rebin_tf(data, tint=inputs.tint, fint=inputs.fint)
        time_arr = np.linspace(time_arr[0], time_arr[-1], data.shape[1])
        freqaxis = np.linspace(freqaxis[0], freqaxis[-1], data.shape[0])
        timeres *= inputs.tint
    # Make plots
    if inputs.plot_all:
        time_arr += start_time_file_unix
        plot_im(data, time_arr, vmax=3, vmin=-2, 
                figname=inputs.fn.strip(ftype)+'_waterfall.pdf')
        plot_dedisp(data, time_arr, dm=inputs.dm,
                    figname=inputs.fn.strip(ftype)+'_dedisp_ts.pdf')
    # Save data to numpy arrays
    if inputs.save_data:
        np.save(inputs.fn.strip(ftype)+'_DM%0.2f' % inputs.dm, data)
        np.save(inputs.fn.strip(ftype)+'timeseries_DM%0.2f' % inputs.dm, data.mean(0))

    print("\nSaved all plots to /home/connor/*pdf\n")

