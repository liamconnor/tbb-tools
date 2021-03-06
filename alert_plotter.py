import os

import numpy as np
import argparse 
import h5py
from datetime import datetime, timedelta
import time 

def get_timefreq(file):
    try:
        timeres=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_0/'].attrs['INCREMENT']
        freqaxis=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_1/'].attrs['AXIS_VALUES_WORLD']
        beamno=0
    except:
        timeres=file['SUB_ARRAY_POINTING_000/BEAM_001/COORDINATES/COORDINATE_0/'].attrs['INCREMENT']
        freqaxis=file['SUB_ARRAY_POINTING_000/BEAM_001/COORDINATES/COORDINATE_1/'].attrs['AXIS_VALUES_WORLD']
        beamno=1

    return timeres, freqaxis, beamno

def read_h5(fn, time_range=(0,5), tint=1, fint=1, freqindex='all'):
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
    print("Unix times in file: %d-%d" % (start_time_file_unix, 
                                         end_time_file_unix))

    if len(time_range)==1:
        time_start = time_range[0]-0.25
        range_err = time_start>=start_time_file_unix and time_range[0]<=end_time_file_unix
        assert range_err, "That unix time is not in the file"
        print("Looking for data at unix time %s" % time_start)
        startsec=time_start-start_time_file_unix
        endsec=startsec+8
    elif len(time_range)==2:
        if time_range[0]>1400000000:
            time_start = time_range[0]-0.25
            range_err = time_start>=start_time_file_unix and time_range[0]<=end_time_file_unix
            assert range_err, "That unix time is not in the file"
            print("Looking for data at unix time %s" % time_start)
            startsec=time_start-start_time_file_unix
            endsec=startsec+time_range[1]
        else:
            startsec, endsec = time_range

    if endsec<=startsec:
        print("Start time is larger than end time")
        exit()
    elif endsec-startsec>30.0:
        print("Are you sure you want %0.2f sec of data?" % (endsec-startsec))

    if freqindex=='all':
        freqmin,freqmax=0,None
    elif len(freqindex)==2:
        freqmin,freqmax=freqindex[0],freqindex[1]
    else:
        print("Expecting either str or len(2) tuple")
        exit()

    timeres, freqaxis, beamno = get_timefreq(file)
    freqaxis = freqaxis[freqmin:freqmax]

    print("\nAssuming time resolution: %0.6f\n" % timeres)

    startsample=int(startsec/timeres)
    endsample=int(endsec/timeres)
    ntime=endsample-startsample
    nfreq=len(freqaxis)
    chunksize=10000.
    if ntime>chunksize:
        data = np.empty([nfreq,ntime])
        print("Reading the file in chunks of %d" % chunksize)
        nchunk=int(np.ceil(ntime/chunksize))
        for jj in range(nchunk):
            startjj,endjj = int(jj*chunksize),int((jj+1)*chunksize)
            data_=file["/SUB_ARRAY_POINTING_000/BEAM_00%d/STOKES_0" % beamno][startsample+startjj:startsample+endjj,freqmin:freqmax]

            if jj<nchunk-1:
                data_=file["/SUB_ARRAY_POINTING_000/BEAM_00%d/STOKES_0" % beamno][startsample+startjj:startsample+endjj,freqmin:freqmax]
                data[:, startjj:endjj] = data_.T
            else:
                data_=file["/SUB_ARRAY_POINTING_000/BEAM_00%d/STOKES_0" % beamno][startsample+startjj:endsample,freqmin:freqmax]  
                data[:, startjj:] = data_.T
    else:
        data=file["/SUB_ARRAY_POINTING_000/BEAM_00%d/STOKES_0" % beamno][startsample:endsample,freqmin:freqmax]
        data=data.T

    if len(data)==0:
        print("No data in specified range")
        exit()

    time_arr = start_time_file_unix + np.linspace(startsec,endsec,ntime)

    # if fint>1 or tint>1:
    #     data = rebin_tf(data, tint=tint, fint=fint)
    #     time_arr = np.linspace(time_arr[0], time_arr[-1], data.shape[1])
    #     freqaxis = np.linspace(freqaxis[0], freqaxis[-1], data.shape[0])

    return data, timeres, time_arr, freqaxis, start_time_file_unix

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

def plot_im(data, time_arr=None, freq=(109863281.25, 187976074.21875), 
            taxis=1, vmax=3, vmin=-2, figname=None):
    """ Plot the 2D time/freq data. Freq is in Hz. vmax and 
    vmin are in sigmas. 
    """
    data_ = data.copy()
    if taxis==0:
        print("Transposing for plotter")
        data_ = data_.T

    if len(freq)>2:
        freq = (freq[0], freq[-1])

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

    ivar = (np.var(data_dd, axis=1)[:, None])**-1
    ivar[ivar==np.inf] = 0
    ivar[ivar!=ivar] = 0

    dd_ts = (ivar*data_dd).mean(0)
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
    parser.add_argument('-o', '--outdir', help='output directory for data and plots',
                        default='/data/projects/COM_ALERT/', type=str)
    parser.add_argument('-s', '--save_data', help='save data to .npy',
                        action='store_true')
    parser.add_argument('-T', '--times', 
                        help='start end times in sec or unix time', nargs='+',
                        type=float, default=(0,5.0))
    parser.add_argument('-p', '--plot_all', 
                        help='Make plots along the way', action='store_true')
    parser.add_argument('--nfchunks', 
                        help='number of freq chunks if arrays are too big', 
                        default=1, type=int)

    inputs = parser.parse_args()
    
    if inputs.plot_all:
        try:
            import matplotlib.pylab as plt
            fig=plt.figure() 
            plt.plot()
            plt.close()
        except:
            print("Cannot plot. Will save figures down.")
            import matplotlib as mpl
            mpl.use('Agg')
            import matplotlib.pyplot as plt

    if inputs.fn[-3:]=='npy':
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

    data_full=[]
    time_arr_full=[]
    freqaxis_full=[]
    file=h5py.File(inputs.fn,'r')
    freqaxis_tot=get_timefreq(file)[1]
    nchantot=len(freqaxis_tot)
    file.close()
    fchunksize=int(np.ceil(nchantot/inputs.nfchunks))
    freq_ref_MHz=1e-6*np.median(freqaxis_tot)
    for chunk in range(inputs.nfchunks):
        if inputs.fn[-2:]=='h5':
            res = read_h5(inputs.fn, inputs.times, tint=inputs.tint,
                          fint=inputs.fint, freqindex=(chunk*fchunksize,(chunk+1)*fchunksize))
            data, timeres, time_arr, freqaxis, start_time_file_unix = res
            ftype='.h5'
        # RFI clean data by zapping bad channels
        if inputs.rfi:
            fignamerfi=inputs.outdir+'/plots/'+inputs.fn.strip(ftype)+'_rfi.pdf'
            data, ind_use, mask = dumb_clean(data, plot_clean=inputs.plot_all, 
                                             figname=fignamerfi)
        # Dedisperse data if given DM > 0
        if inputs.dm>0:
            data = dedisperse(data, inputs.dm, freq=freqaxis,  
                              timeres=timeres, freq_ref=freq_ref_MHz)

        # Downsample / channelize data
        if inputs.fint>1 or inputs.tint>1:
            data = rebin_tf(data, tint=inputs.tint, fint=inputs.fint)
            time_arr = np.linspace(time_arr[0], time_arr[-1], data.shape[1])
            freqaxis = np.linspace(freqaxis[0], freqaxis[-1], data.shape[0])
            timeres *= inputs.tint

        data_full.append(data)
        freqaxis_full.append(freqaxis)

    data = np.concatenate(data_full, axis=0)
    freqaxis = np.concatenate(freqaxis_full)

    # Make plots
    if inputs.plot_all:
        fignameim=inputs.outdir+'/plots/'+inputs.fn.strip(ftype)+'_waterfall.pdf'
        fignamedd=inputs.outdir+'/plots/'+inputs.fn.strip(ftype)+'_dedisp_ts.pdf'
        plot_im(data, time_arr, freq=1e6*freqaxis, vmax=3, vmin=-2, figname=fignameim)
        plot_dedisp(data, time_arr, dm=inputs.dm, figname=fignamedd)

    # Save data to numpy arrays
    if inputs.save_data:
        np.save(inputs.outdir+'/dedispbf/'+inputs.fn.strip(ftype)+'_DM%0.2f' % inputs.dm, data)
        np.save(inputs.outdir+'/dedispbf/'+inputs.fn.strip(ftype)+'timeseries_DM%0.2f' % inputs.dm, data.mean(0))
        print("\nSaved all data to %s\n" % (inputs.outdir+'/dedispbf/'))

    if inputs.plot_all:
        print("\nSaved all plots to %s\n" % (inputs.outdir+'/plots/'))








