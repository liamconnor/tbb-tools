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

        return rcus, data_rcu

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

  def get_stations_present(self):
    f, attrs = self.read_h5()
    stations = attrs['OBSERVATION_STATIONS_LIST']

    for stat in stations:
      try:
        stations_list.append(f[stat])
      except:
        continue

    return stations_list

  def get_rcus_present(self):
    pass





  # Goal: Create array like that produced by TBB_Rawdata.construct_rcu_arr



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




