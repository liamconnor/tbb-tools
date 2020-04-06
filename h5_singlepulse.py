import numpy as np
import argparse
import h5py

import filterbank

filhdr = {'telescope_id': 4,
      'az_start': 0.0,
      'nbits': 8,
      'source_name': 'B0329+54',
      'data_type': 1,
      'nchans': 0,
      'machine_id': 15,
      'tsamp': 0.001,
      'foff': -0.0244,
      'src_raj': 0,
      'src_dej': 0,
      'tstart': 58523.3437492,
      'nbeams': 1,
      'fch1' : 800.,
      'za_start': 0.0,
      'rawdatafile': '',
      'nifs': 1,
      #'nsamples': 12500
      }

def read_fil_data(fn, start=0, stop=1):
     print("Reading filterbank file %s \n" % fn)
     fil_obj = filterbank.FilterbankFile(fn)
     header = fil_obj.header
     delta_t = fil_obj.header['tsamp'] # delta_t in seconds
     fch1 = header['fch1']
     nchans = header['nchans']
     foff = header['foff']
     fch_f = fch1 + nchans*foff
     freq = np.linspace(fch1, fch_f, nchans)

     try:
          data = fil_obj.get_spectra(start, stop)
     except(ValueError):
          data = 0
     # turn array into time-major, for preprocess
#    data = data.transpose() 

     return data, freq, delta_t, header

def create_new_filterbank(fnh5, fn_fil_out, telescope='LOFAR'):
   file = h5py.File(fnh5,'r')
   dt=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_0/'].attrs['INCREMENT']
   freqaxis=file['SUB_ARRAY_POINTING_000/BEAM_000/COORDINATES/COORDINATE_1/'].attrs['AXIS_VALUES_WORLD']
   freqaxis *= 1e-6
   nchans=len(freqaxis)
   foff=np.diff(freqaxis)[0]

   filhdr['tsamp'] = dt
   filhdr['nchans'] = nchans
   filhdr['fch1'] = freqaxis[-1]
   filhdr['foff'] = foff

   try:
      import sigproc
      filhdr['rawdatafile'] = fn_fil_out

      newhdr = ""
      newhdr += sigproc.addto_hdr("HEADER_START", None)
      for k,v in filhdr.items():
          newhdr += sigproc.addto_hdr(k, v)
      newhdr += sigproc.addto_hdr("HEADER_END", None)
      print("Writing new header to '%s'" % fn_fil_out)
      outfile = open(fn_fil_out, 'wb')
      outfile.write(newhdr)
      spectrum = np.zeros([filhdr['nchans']], dtype=np.uint8)
      outfile.write(spectrum)
      outfile.close()
   except:
      print("Either could not load sigproc or create filterbank")

def write_to_fil(data, header, fn):
     filterbank.create_filterbank_file(
          fn, header, spectra=data, mode='readwrite')
     print("Writing to %s" % fn)

def read_fil_data(fn, start=0, stop=1e7):
     print("Reading filterbank file %s \n" % fn)
     fil_obj = filterbank.FilterbankFile(fn)
     header = fil_obj.header
     delta_t = fil_obj.header['tsamp'] # delta_t in seconds
     fch1 = header['fch1']
     nchans = header['nchans']
     foff = header['foff']
     fch_f = fch1 + nchans*foff
     freq = np.linspace(fch1, fch_f, nchans)

     try:
          data = fil_obj.get_spectra(start, stop)
     except(ValueError):
          data = 0
     # turn array into time-major, for preprocess
#    data = data.transpose() 

     return data, freq, delta_t, header

def h5_to_fil(fnh5, fn_fil_out, nchunk='all'):
   file = h5py.File(fnh5,'r')
   chunksize = int(1e4)

   if nchunk=='all':
     nchunk=int(1e6)

   for ii in range(nchunk):
     print(ii, nchunk)
     startsample, endsample = ii*chunksize, (ii+1)*chunksize
     data = file["/SUB_ARRAY_POINTING_000/BEAM_000/STOKES_0"][startsample:endsample,:]
     data /= 1e12
     data *= 50
     data = data.astype(np.int8)

     if len(data)==0:
          print("Empty data")
          continue

     if ii==0:
          header = read_fil_data(fn_fil_out, start=0, stop=1)[-1]
          fn_rfi_clean = write_to_fil(data, header, fn_fil_out)
     elif ii>0:
          print("Writing chunk %d" % ii)
          fil_obj = filterbank.FilterbankFile(fn_fil_out, mode='readwrite')
          fil_obj.append_spectra(data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                     description="h5 to filterbank",
                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-fnh5', '--fnh5', 
                        help='filename (assumes h5 for now)', 
                        type=str, required=True)
    parser.add_argument('-fnfil', '--fnfil', 
                        help='filename (assumes h5 for now)', 
                        type=str, required=True)
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
    create_new_filterbank(inputs.fnh5, inputs.fnfil, telescope='LOFAR')
    h5_to_fil(inputs.fnh5, inputs.fnfil, nchunk='all')
    arg_str = (inputs.fnh5.strip('.h5'), inputs.dm, inputs.fnfil)
    os.system('prepdata -o %s -dm %f %s' % arg_str)



