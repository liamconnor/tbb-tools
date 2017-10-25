"""
Code to read and manipulate TBB packets. 
The data are organized in 2040 byte frames, each 
of which starts with an 88 byte header of format 
'BBBBIIIHH64BHH'. That is followed by a 1948 byte 
data frame. 
"""


import numpy as np
import struct

offs=0
nfram=2040
fmt = 'BBBBIIIHH64BHH'
header_len = 88 # number of bytes in header
data_len = 1948 # number of bytes in data packet
fn = 'sb_20170711_094130_1310.dat'
fn = 'sb_20170711_094130_1310_copy.dat'
fn = './test4.dat'

# 013001010 is station 13, RSP 001, RCU 010 

class TBB_rawdata():
    """ Class to deal with the raw data from the 
    transient buffer boards
    """
    def __init__(self):
        self.filename = None
        self.fmt = 'BBBBIIIHH64BHH'
        self.nframe = 2040
        self.header_len = 88 # number of bytes in header
        self.data_len = 1948 # number of bytes in data packet
        self.filename = None

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
        data_bi = f.read(self.nframe)

        if len(data_bi)==0:
            return

        header = self.parse_header(data_bi[:self.header_len])
        data = np.fromstring(data_bi[self.header_len:\
                            self.header_len+self.data_len],\
                            dtype=np.int16)
        final4 = data_bi[-4:]

        return data, header, final4

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
      this opreation is done in (bs[i] & 1<<j)
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

    def write_for_new_image(self, fn, fnout='', 
                nband=2, nframe=10, RCU_ID=None, 
                subbands=None, RSP_ID=None, crc16=None):
        g = open(fnout, 'a')

        # all data frames for a single subband, then 
        # all dipoles 

        # if subbands is None:
        #   subbands = range(1)
        # if RCU_ID is None:
        #   RCU_ID = range(1)
        # if nframes is None:
        #   nframes = range(1)

#        for xx in RCU_ID:

        H = []
        D = []
        for ii in subbands:
            f = open(fn, 'r+')

            for ff in range(nframe):
#            for ii in subbands:
                try:
                  r_tup = self.read_data(f)
                  print r_tup[1][:7]
                except:
                  continue 
                d, h, final4 = r_tup

                H.append(h)
                D.append(D)

                h = self.alter_header(h, RCU_ID=RCU_ID, band=ii, RSP_ID=RSP_ID, crc16=crc16)
                h = self.pack_header(h)

                g.write(h)  
                g.write(d)
                g.write(final4)

            f.close()

        return H, D

t = TBB_rawdata()
H, D = t.write_for_new_image(fn, fnout='./out_4dip.dat', nframe=10000, \
                      RCU_ID=None, subbands=[100, 200, 300], \
                      RSP_ID=None, crc16=None)
# os.system('rm -rf out3.dat')
# t.write_for_new_image(fn, fnout='./out_crc16err.dat', nframe=10, RCU_ID=[10], subbands=[16], RSP_ID=1, crc16=21874)
# os.system('scp out_crc16err.dat lofar:~/')

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


# def doit(d, a):

#     for ii in range(len(d)):
#         if d[ii]==a[0]:
#             if d[ii+1]==a[1]:
#                 if d[ii+2]==a[2]:
#                     if d[ii+3]==a[3]:
#                         print d[ii:ii+8], ii
#                         return ii

# def parse_header(header_str, fmt):
#     header = struct.unpack(fmt, header_str)

#     return header 

# def read_data(f):
#     """ Read in one frame of data, length nfram bytes
#     and return data, header, and final four bytes (whose 
#         purpose I don't yet know)
#     """
#     data_bi = f.read(nfram)

#     if len(data_bi)==0:
#         return

#     header = parse_header(data_bi[:header_len], fmt)
#     data = np.fromstring(data_bi[header_len:header_len+data_len], dtype=np.int16)
#     final4 = data_bi[-4:]

#     return data, header, final4

# def pack_header(header, fmt):
#     s = struct.Struct(fmt)
#     packed_header_data = s.pack(*header)

#     return packed_header_data

# def change_band(header, fmt, newval=8):
#     # Replace sub-band number 
#     values = header[0:21]+(newval,)+header[22:]

#     packed_header_data = pack_header(values, fmt)
#     #header_ = header[:offs] + packed_data + header[offs+header_len:]

#     return packed_header_data

# def write_for_new_image(fn, fnout='', nband=2, nframe=10):
#     f = open(fn, 'r+')
#     g = open(fnout, 'a')

#     antennas = [1, 2]

#     for xx in antennas:

#     for ii in range(nframe):
#         r_tup = read_data(f)
#         d, h, final4 = r_tup
#         h = change_band(h, fmt, newval=8)
#         print d[0],ii
#         g.write(h)
#         g.write(d)
#         g.write(final4)

#     f = open(fn, 'r+')

#     for ii in range(nframe):
#         r_tup = read_data(f)
#         d, h, final4 = r_tup
#         h = change_band(h, fmt, newval=16)
#         g.write(h)
#         g.write(d)
#         g.write(final4)


# f = open(fn, 'r')

# D, H = [], []

# fnout = 'sb_20170711_094130_1310_output_twoband_oneheader.dat'
# g = open(fnout, 'a')

# while True:
#     r_tup = read_data(f)

#     if r_tup==None:
#         break
#     else:
#         d, h, final4 = r_tup
# #       print d
# #       h_ = change_band(h, fmt, newval=16)
#         h = change_band(h, fmt, newval=16)

# #       print parse_header(h, fmt)
# #       print parse_header(h_, fmt)

#         g.write(h)
#         g.write(d)
#         g.write(final4)

#         print h 
#         print d 
#         print final4

# #       g.write(h_)
#         # g.write(d)
#         # g.write(final4)   


# fn = 'sb_20170711_094130_1310.dat'

# f = open(fn, 'r')
# data = f.read(2400)
# struct.unpack(fmt,data[offs:offs+88])

# d = struct.unpack(fmt,data[offs:offs+88])
# dd = np.fromstring(data[offs+88:offs+nfram], dtype=np.int16)
# f.close()

# # Replace sub-band number
# values = d[0:21]+(8,)+d[22:]
# s = struct.Struct(fmt)
# packed_data = s.pack(*values)
# data_ = data[:offs] + packed_data + data[offs+88:]

# assert len(data_)==len(data)
# fnout = 'sb_20170711_094130_1310_output_test.dat'
# g = open(fnout, 'a')
# g.write(data_)
# g.close()


