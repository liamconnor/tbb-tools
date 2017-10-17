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
        data = np.fromstring(data_bi[self.header_len:self.header_len+self.data_len], dtype=np.int16)
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
            h = self.change_band(header, newval=band)[0]
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

    def pack_header(self, header):
        assert type(header) != str, "Already in hex"
        s = struct.Struct(self.fmt)
        packed_header_data = s.pack(*header)

        return packed_header_data

    def change_band(self, header, newval=8):
        """ header should be tuple / array
        """
        # Replace sub-band number 
        values = header[0:21]+(newval,)+header[22:]

        packed_header_data = self.pack_header(values)
        #header_ = header[:offs] + packed_data + header[offs+header_len:]

        return values, packed_header_data

    def write_for_new_image(self, fn, fnout='', 
                nframes=10, nband=2, nframe=10, RCU_ID=[10], 
                subbands=[16], RSP_ID=1, crc16=None):
        g = open(fnout, 'a')

        # all data frames for a single subband, then 
        # all dipoles 

        for xx in RCU_ID:
            for ii in subbands:
                f = open(fn, 'r+')

                for ff in range(nframes):
                    r_tup = self.read_data(f)
                    d, h, final4 = r_tup
                    h = self.alter_header(h, RCU_ID=xx, band=ii, RSP_ID=RSP_ID, crc16=crc16)
                    print h
                    h = self.pack_header(h)
                    g.write(h)  
                    g.write(d)
                    g.write(final4)

                f.close()

t = TBB_rawdata()

os.system('rm -rf out3.dat')
t.write_for_new_image(fn, fnout='./out_crc16err.dat', nframe=10, RCU_ID=[10], subbands=[16], RSP_ID=1, crc16=21874)
os.system('scp out_crc16err.dat lofar:~/')


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


def doit(d, a):

    for ii in range(len(d)):
        if d[ii]==a[0]:
            if d[ii+1]==a[1]:
                if d[ii+2]==a[2]:
                    if d[ii+3]==a[3]:
                        print d[ii:ii+8], ii
                        return ii

def parse_header(header_str, fmt):
    header = struct.unpack(fmt, header_str)

    return header 

def read_data(f):
    """ Read in one frame of data, length nfram bytes
    and return data, header, and final four bytes (whose 
        purpose I don't yet know)
    """
    data_bi = f.read(nfram)

    if len(data_bi)==0:
        return

    header = parse_header(data_bi[:header_len], fmt)
    data = np.fromstring(data_bi[header_len:header_len+data_len], dtype=np.int16)
    final4 = data_bi[-4:]

    return data, header, final4

def pack_header(header, fmt):
    s = struct.Struct(fmt)
    packed_header_data = s.pack(*header)

    return packed_header_data

def change_band(header, fmt, newval=8):
    # Replace sub-band number 
    values = header[0:21]+(newval,)+header[22:]

    packed_header_data = pack_header(values, fmt)
    #header_ = header[:offs] + packed_data + header[offs+header_len:]

    return packed_header_data

def write_for_new_image(fn, fnout='', nband=2, nframe=10):
    f = open(fn, 'r+')
    g = open(fnout, 'a')

    antennas = [1, 2]

    for xx in antennas:

    for ii in range(nframe):
        r_tup = read_data(f)
        d, h, final4 = r_tup
        h = change_band(h, fmt, newval=8)
        print d[0],ii
        g.write(h)
        g.write(d)
        g.write(final4)

    f = open(fn, 'r+')

    for ii in range(nframe):
        r_tup = read_data(f)
        d, h, final4 = r_tup
        h = change_band(h, fmt, newval=16)
        g.write(h)
        g.write(d)
        g.write(final4)


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


