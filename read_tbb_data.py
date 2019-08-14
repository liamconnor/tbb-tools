"""
Code to read and manipulate TBB packets. 
The data are organized in 2040 byte frames, each 
of which starts with an 88 byte header of format 
'BBBBIIIHH64BHH'. That is followed by a 1948 byte 
data frame. 

Example usage:
==============

fn = './data/sb_20171004_121157_10131.dat'

TBB = TBB_Rawdata()

Suppose you want to read in all the data after parsing headers:

headers, data, crc32 = TBB.read_all_data(fn, print_headers=True)

Suppose you want to print the first 10 frames, including 
their data payloads, of some file:

TBB.print_frames(fn, print_data=True)

Suppose you want to alter a data file, forcing it to have 
multiple subbands (in this case 100 and 300), 
but not changing its dipoles:

headers, data, crc32 = TBB.write_for_new_image(fn, 
                                  fnout='./outfile.dat', nframe=10000,
                                  RCU_ID=None, subbands=[100, 300],
                                  RSP_ID=None, crc16=None)

"""


import numpy as np
import struct
import h5py


class TBB_Rawdata():
    """ Class to deal with the raw data from the 
    transient buffer boards
    # 013001010 is station 13, RSP 001, RCU 010 
    """
    def __init__(self, version='new'):
        self.filename = None

        self.fmt = 'BBBBIIIHH64BHH'
        self.bytes_per_frame = 2040
        self.header_len = 88 # number of bytes in header
        self.data_len = 1948 # number of bytes in data packet
        self.filename = None

        if version is 'new':
          self.bytes_per_frame = 2012
          self.data_len = 1920

    def header_dict(self):
        # dictionary whose keys are the frame parameters 
        # key values are index in header array
        header_dict = { 'station_ID'    : 0,
                        'RSP_ID'        : 1,
                        'RCU_ID'        : 2,
                        'sample_freq'   : 3,
                        'seq_num'       : 4,
                        'utime'         : 5,
                        'samp_per_freq' : 7,
                        'crc16'         : -1,
                        'band_num'      : 999,
                        'slice_num'     : 999,
                       }

        return header_dict

    def f_open(fname):
        return open(fname, 'r')

    def parse_header(self, header_str):
        header = struct.unpack(self.fmt, header_str)

        return header 

    def read_data(self, f):
        """ Read in one frame of data, length nfram bytes
        and return data, header, and final four bytes (whose 
            purpose I don't yet know)
        """
        data_bi = f.read(self.bytes_per_frame)

        if len(data_bi)==0:
            return

        header = self.parse_header(data_bi[:self.header_len])
        data = np.fromstring(data_bi[self.header_len:\
                            self.header_len+self.data_len],\
                            dtype=np.int16)
        crc32 = data_bi[-4:]

        return header, data, crc32

    def read_all_data(self, fn, nframe=10000, print_headers=False):
        """ Take filename fn, read in nframe frames 
        and return 
        a list of headers (header_full), 
        a list of payloads (data_full), 
        and a list of crc32s (crc32_full).
        """
        f = open(fn,'r')

        data_full = []
        header_full = []
        crc32_full = []

        for ii in range(nframe):
            try: 
                header, data, crc32 = self.read_data(f)
          
                data_full.append(data)
                header_full.append(header)
                crc32_full.append(crc32)

                if print_headers is True:
                    self.print_one_frame(header)
            except:
                break

        return np.array(header_full), np.array(data_full), np.array(crc32_full)

    def construct_rcu_arr(self, data_full, header_full):
        """ Take data and header arrays and reshape into 
        an np.arr of shape (ncru, nframe, nsample_per_frame)

        returning the rcu list and data array
        """
        rcus_full = header_full[:, 2]
        rcus = list(set(header_full[:,2]))
        nrcu = len(rcus)
        samples_per_frame = data_full.shape[-1]
        nframe = data_full.shape[0]//nrcu
        
        data_rcu = np.zeros([nrcu, nframe, samples_per_frame], dtype=data_full.dtype)
        
        for ii, rcu in enumerate(rcus):
            data_rcu[ii] = data_full[rcus_full==rcu]

        print(data_rcu.shape)

        return data_rcu, rcus

    def alter_header(self, header, 
                band=None, station_ID=None, RSP_ID=None, 
                sample_freq=None, RCU_ID=None, crc16=None):
        """  Set header 
        """
        h = np.array(header)
        header_dict = self.header_dict()

        if band is not None:
            h = self.change_band(header, bands=band)
            h = np.array(h)
        if station_ID is not None:
            h[header_dict['station_ID']] = station_ID

        if RSP_ID is not None:
            h[header_dict['RSP_ID']] = RSP_ID

        if RCU_ID is not None:
            h[header_dict['RCU_ID']] = RCU_ID

        if sample_freq is not None:
            h[header_dict['sample_freq']] = sample_freq

        if crc16 is not None:
            h[header_dict['crc16']] = crc16

        return tuple(h)

    def get_bands_present(self, header):
      """ Find bands present given header. Based 
      on code in 
      https://svn.astron.nl/LOFAR/branches/TBB-DataWriter-11012/
                RTCP/Cobalt/OutputProc/scripts/tbb-printframes.cc

      As far as I can tell, the mapping is as follows:
      bytes 9:73 in the header contain the bandSel 
      
      bandSel = [b1, b2, ..., b64]
      b1 encodes bands [0, 1, ..., 7]
      b64 encodes bands [504, ..., 511]

      if b1 == 1, then log2(b1) is the bandnumber. In other 
      words, each element of bandSel can only be [1, 2, 4, 16, 32, 64, 128] 
      this operation is done in (bs[i] & 1<<j)
      """

      header = np.array(header)
      bandSel = header[9:9+64]
      bands_present = []

      for i in range(64):
        for j in range(0, 8):
          if (bandSel[i] & 1<<j):
            bands_present.append(8*i+j)

      return bands_present

    def get_band_bit(self, bands):
      """ Inverse of self.get_bands_present 
      Give a band number between 0-511, returns 
      index of 64-element bandSel array and the value. 

      E.g. band 200 should return (25, 1) meaning 
      bandSel[25] = 1
      """
      index, vals = bands//8, 2**(bands%8)

      return index, vals

    def pack_header(self, header):
        assert type(header) != str, "Already in hex"
        s = struct.Struct(self.fmt)
        packed_header_data = s.pack(*header)

        return packed_header_data

    def get_band_array(self, header_full, nframe=np.inf):
        sb_arr = []
        n = min(len(header_full), nframe)

        for ii in xrange(n//1):
            if ii>0:
                sb2 = self.get_bands_present(header_full[ii*1])[0]
                if sb!=sb2:
                    print(sb, sb2, ii)
            sb = self.get_bands_present(header_full[ii*1])[0]
            sb_arr.append(sb)

        sb_arr = np.array(sb_arr)

        return sb_arr

    def change_band(self, header, bands=None):
        """ header should be tuple / array
        """

        if bands is None:
          return header

        if type(bands)==list:
          bands = np.array(bands)
        if type(header) in (list, tuple):
          header = np.array(header)


        index, vals = self.get_band_bit(bands)

        header[9:9+64] *= 0
        header[9:9+64][index] = vals 

        # # Replace sub-band number 
        # values = header[0:21]+(newval,)+header[22:]

        # packed_header_data = self.pack_header(values)
        # #header_ = header[:offs] + packed_data + header[offs+header_len:]

        return header

    def print_one_frame(self, header):
        """ Imitating Alexander's cpp code
        """
        print "-------------------------------"
        print "Station ID:       %12d" % header[0]
        print "RSP ID:           %12d" % header[1]
        print "RCU ID:           %12d" % header[2]
        print "Sample Freq:      %12d" % header[3]
        print "Seq Nr:           %12d" % header[4]
        print "Time:             %12d" % header[5]
        print "Band Nr:          %12d" % 0 # For now
        print "Slice Nr:         %12d" % header[6]
        print "NSamples/fr:      %12d" % header[7]
        print "NFreq Bands:      %12d" % header[8]
        print "Band(s) present:  %12d" % self.get_bands_present(header)[0]
        print "-------------------------------"

    def print_frames(self, fn, print_data=False, nframe=1000):
        f = open(fn)

        for ii in range(nframe):
            try:
                header, data, crc = self.read_data(f)
                self.print_one_frame(header)

                if print_data==True:
                    print data
            except:
                return 

    def write_to_file(self, header_list, data_list, crc32_list, fnout):
        """ Take a list of headers, list of frames, and crc32 list 
            and write to file fnout 

            example: 
        """
        f = open(fnout, 'w+')
        nframe = len(header_list)

        for ii in np.arange(nframe)[:]:
            h = self.pack_header(header_list[ii])
            f.write(h)
            f.write(data_list[ii])
            f.write(crc32_list[ii])

        print "Done writing to %s" % fnout


    def write_for_new_image(self, fn, fnout='', 
                nband=2, nframe=10, RCU_ID=None, 
                subbands=None, RSP_ID=None, crc16=None):
        g = open(fnout, 'a')

        header_list = []
        data_list = []
        crc32_list = []

        for ii in subbands:
            f = open(fn, 'r+')

            for ff in range(nframe):
#            for ii in subbands:
                r_tup = self.read_data(f)
                try:
                  r_tup = self.read_data(f)
                  if r_tup is None:
                    continue
                except:
                  continue 
                h, d, crc32 = r_tup

                data_list.append(d)
                crc32_list.append(crc32)
                h = self.alter_header(h, RCU_ID=RCU_ID, \
                    band=ii, RSP_ID=RSP_ID, crc16=crc16)

                header_list.append(h)
                h = self.pack_header(h)

                g.write(h)  
                g.write(d)
                g.write(crc32)

            f.close()

        return header_list, data_list, crc32_list

class TBBh5_Reader():

  def __init__(self, fn):
    self.fn = fn
    self.nsubband_full = 512
    self.delay_constant = 4148.808
    self.time_per_sample = 5.12e-6 
    h5_attrs = [u'DOC_VERSION',
               u'PROJECT_CO_I',
               u'FILETYPE',
               u'NOTES',
               u'DOC_NAME',
               u'FILENAME',
               u'OBSERVATION_NOF_STATIONS',
               u'OBSERVATION_FREQUENCY_MAX',
               u'OBSERVATION_END_MJD',
               u'FILTER_SELECTION',
               u'OBSERVATION_STATIONS_LIST',
               u'OBSERVATION_START_UTC',
               u'FILEDATE',
               u'TELESCOPE',
               u'ANTENNA_SET',
               u'OBSERVATION_END_UTC',
               u'PROJECT_PI',
               u'TARGETS',
               u'OBSERVATION_ID',
               u'OPERATING_MODE',
               u'PROJECT_CONTACT',
               u'SYSTEM_VERSION',
               u'OBSERVATION_NOF_BITS_PER_SAMPLE',
               u'GROUPTYPE',
               u'CLOCK_FREQUENCY',
               u'OBSERVATION_FREQUENCY_CENTER',
               u'PROJECT_ID',
               u'OBSERVATION_FREQUENCY_MIN',
               u'OBSERVATION_FREQUENCY_UNIT',
               u'CLOCK_FREQUENCY_UNIT',
               u'PROJECT_TITLE',
               u'OBSERVATION_START_MJD',
               u'NOF_STATIONS']

  def read_h5(self):
    f = h5py.File(self.fn,'r')

    return f, f.attrs 

  def get_stations_groups(self):
    f, attrs = self.read_h5()
    stations = attrs['OBSERVATION_STATIONS_LIST']
    stations_groups = []
    stations_list = []

    for stat in stations:
      try:
        statname = 'STATION_%s' % stat[:5]
        stations_groups.append(f[statname])
        stations_list.append(statname)
      except:  
        continue

    if len(stations_list)==0:
        print("Could not find the correct station")
 
    return stations_groups, stations_list

  def get_rcus_present(self, station_group):

    dipole_names = station_group.keys()

    rcus = []

    for dn in dipole_names:
      rcus.append(dn[-2:])

    return rcus, dipole_names

  def print_summary_data(self, fn):
    station_groups, stations_list = self.get_stations_groups()
    for station_group in station_groups:
      rcus, dipole_names = self.get_rcus_present(station_group)
      print("===========DIPOLES for %s==========" % station_group.name)
      for dd in dipole_names:
        print(dd)
        tarr = []
        print("SUB-BANDS:")
        for sb in np.sort(station_group[dd].keys()):
          t0_int = station_group[dd][sb].attrs['TIME']
          slice0 = station_group[dd][sb].attrs['SLICE_NUMBER']
          t0 = t0_int + self.time_per_sample*slice0
          tarr.append(t0)
          #print(station_group[dd][sb].attrs.items())
#        for sb in station_group.items()[-1][-1].keys():
          print("     %s t0: %f" % (sb, t0))
        print("     max delay: %f" % np.abs(tarr[-1]-tarr[0]))
      print("===========DIPOLES for %s==========\n\n" % station_group.name)

  def get_time_arr(self, fn):
    station_groups, stations_list = self.get_stations_groups()
    tarr = []
    for station_group in station_groups:
      rcus, dipole_names = self.get_rcus_present(station_group)
      for dd in dipole_names:
        for sb in np.sort(station_group[dd].keys()):
          t0_int = station_group[dd][sb].attrs['TIME']
          slice0 = station_group[dd][sb].attrs['SLICE_NUMBER']
          t0 = t0_int + self.time_per_sample*slice0
          tarr.append(t0)

    return tarr

  def voltage_to_intensity(self, data_arr_volt):
    """ Assume last axis is time (Re, Im, Re, ...)
    """

    data_arr_int = data_arr_volt[..., ::2]**2 + data_arr_volt[..., 1::2]**2

    return data_arr_int

  def station_data(self, station_name, rebin=True):
    """ Construct voltage array of all data 
    in given station. The output array will be 
    (nrcu, nsubband, nsamples). The dipole/SB map 
    is also returned
    """
    f, attrs = self.read_h5()
    stations_group = f[station_name]
    rcus, dipole_names = self.get_rcus_present(stations_group)

    ind = np.argsort(rcus)
    dipole_names = np.array(dipole_names)[ind]
    rcus = np.array(rcus)[ind]

    nrcu = len(rcus)

    data_rcu = []
    dipole_sb_map = []
    t0_alldipoles = []

    for zz, dipole_name in enumerate(dipole_names[:10]):
      print("Dipole:%s Number:%d" % (dipole_name, zz))
      dipole_groups = stations_group[dipole_name]
      t0_sb = []

      for sb_name in np.sort(dipole_groups):
        data_r = dipole_groups[sb_name]['real']
        data_i = dipole_groups[sb_name]['imag']        
        
        nsample = len(data_r)

        try:
          t0_int = dipole_groups[sb_name].attrs['TIME']
          slice0 = dipole_groups[sb_name].attrs['SLICE_NUMBER']
          t0 = t0_int + self.time_per_sample*slice0
          tarr = np.linspace(t0, t0+nsample*self.time_per_sample, nsample)
          t0_sb.append(t0)
        except:
          print("Could not get time data")

        if rebin:
          data_intensity = data_r**2 + data_i**2
          data_intensity = data_intensity.reshape(-1, 10).mean(-1)
          data_rcu.append(data_intensity)
          print(self.time_per_sample)
          self.time_per_sample *= 10
          print(self.time_per_sample)
        else:
          data_complex = np.empty([2*nsample])
          data_complex[::2] = data_r
          data_complex[1::2] = data_i 
          data_rcu.append(data_complex)

        dipole_sb_map.append([dipole_name, sb_name])
        
      t0_alldipoles.append(t0_sb)

    rcu_set = list(set(dipole_names))
    t0_alldipoles = np.concatenate(t0_alldipoles)
    np.save('output_data', data_rcu)
#    print(data_rcu.shape, len(dipole_sb_map))    
#    exit()
    if rebin:
      nsample /= 2

    data_arr = self.construct_fullband_arr(data_rcu, dipole_sb_map, 
                                           nsample, nrcu, rcu_set, t0_alldipoles)

    np.save('fullarr', data_arr)

    #data_rcu = np.concatenate(data_rcu)  
    #data_rcu = data_rcu.reshape(nrcu, -1, 2*nsample)

    return data_arr, dipole_sb_map, t0_alldipoles

  def construct_array_sweep(self):
    station_name = list(set(self.get_stations_groups()[1]))
    data, mapping, t0_alldipoles = T.station_data(station_name[0])

  def construct_fullband_arr(self, data_rcu, dipole_sb_map, 
                             nsample, nrcu, rcu_set, t0_alldipoles):
    """ Take list of voltage data arrays len(data_rcu)=N_dipole_SBs, 
    data_rcu[0].shape = 
    """


    t0_min = t0_alldipoles.min()
    t0_max = t0_alldipoles.max()
    print(t0_min, t0_max, self.time_per_sample)
    offset = int((t0_max - t0_min)/self.time_per_sample)

    data_arr = np.empty([nrcu, self.nsubband_full, 2*(offset+nsample)])

    rcu_set = np.sort(rcu_set)

    for ii in range(len(dipole_sb_map)):
      sb = int(dipole_sb_map[ii][1][-3:])
      nsample_ii = len(data_rcu[ii])
      rcu_ind = np.where(rcu_set==dipole_sb_map[ii][0])[0][0]
      offset_ii = int((t0_alldipoles[ii] - t0_min)/self.time_per_sample)
      data_arr[rcu_ind, sb, offset_ii:offset_ii+nsample_ii] = data_rcu[ii]

    return data_arr

def compare_data(fnh5, fndat):

  Traw = TBB_Rawdata('new')

  h, d, crc32_full = Traw.read_all_data(fndat, nframe=10000000)
  data_arr_dat, rcus = Traw.construct_rcu_arr(d, h)

  Th5 = TBBh5_Reader(fnh5)

  station_group = Th5.get_stations_groups()[0][0]
  data_arr_h5, rcu_map, t0_alldipoles = Th5.station_data(station_group)

  return data_arr_h5, data_arr_dat


def TBB_Writer_attrs():
    h5_attrs = [u'DOC_VERSION',
             u'PROJECT_CO_I',
             u'FILETYPE',
             u'NOTES',
             u'DOC_NAME',
             u'FILENAME',
             u'OBSERVATION_NOF_STATIONS',
             u'OBSERVATION_FREQUENCY_MAX',
             u'OBSERVATION_END_MJD',
             u'FILTER_SELECTION',
             u'OBSERVATION_STATIONS_LIST',
             u'OBSERVATION_START_UTC',
             u'FILEDATE',
             u'TELESCOPE',
             u'ANTENNA_SET',
             u'OBSERVATION_END_UTC',
             u'PROJECT_PI',
             u'TARGETS',
             u'OBSERVATION_ID',
             u'OPERATING_MODE',
             u'PROJECT_CONTACT',
             u'SYSTEM_VERSION',
             u'OBSERVATION_NOF_BITS_PER_SAMPLE',
             u'GROUPTYPE',
             u'CLOCK_FREQUENCY',
             u'OBSERVATION_FREQUENCY_CENTER',
             u'PROJECT_ID',
             u'OBSERVATION_FREQUENCY_MIN',
             u'OBSERVATION_FREQUENCY_UNIT',
             u'CLOCK_FREQUENCY_UNIT',
             u'PROJECT_TITLE',
             u'OBSERVATION_START_MJD',
             u'NOF_STATIONS']

    return h5_attrs

def lofar_usage_ICD_20227():
    ICD_attrs = ['GROUPTYPE','FILENAME','FILEDATE','FILETYPE','TELESCOPE','PROJECT_ID','PROJECT_TITLE',
          'PROJECT_PI','PROJECT_CO_I','PROJECT_CONTACT','OBSERVATION_ID','OBSERVATION_START_UTC',
          'OBSERVATION_START_MJD','OBSERVATION_END_MJD','OBSERVATION_END_UTC','OBSERVATION_NOF_STATIONS', 'OBSERVATION_STATIONS_LIST',
          'OBSERVATION_FREQUENCY_MIN','OBSERVATION_FREQUENCY_CENTER','OBSERVATION_FREQUENCY_MAX',
          'OBSERVATION_FREQUENCY_UNIT','OBSERVATION_NOF_BITS_PER_SAMPLE','CLOCK_FREQUENCY','CLOCK_FREQUENCY_UNIT',
          'ANTENNA_SET','FILTER_SELECTION','TARGETS','SYSTEM_VERSION','PIPELINE_NAME','PIPELINE_VERSION',
          'DOC_NAME', 'DOC_VERSION','NOTES']

    # In recent .h5 files written by the TBB_Writer it seems 
    # PIPELINE_NAME, PIPELINE_VERSION are in ICD keys, but not in 
    # file

    # OPERATING_MODE
    # NOF_STATIONS

    return ICD_attrs 





