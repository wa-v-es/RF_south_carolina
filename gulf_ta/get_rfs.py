import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client
import numpy as np
import datetime
import os
import sys
from rf import read_rf, rfstats,RFStream
from rf import iter_event_data, IterMultipleComponents
from tqdm import tqdm
########
def signoise(self):#, winsig, winnoise, relative='onset'):
        """
        Determine signal noise ratio by dividing the maximum in the two windows.
        """
        st = self.stats
        self_copy = self.copy()
        self_copy.detrend().taper(max_percentage=0.05)
        self_copy.filter('bandpass', freqmin=0.1, freqmax=1)#,corners=2, zerophase=True)
        winsig=[-5,25] #signal window
        winnoise=[-45,-15] #noise window
        # if relative in ('onset', 'middle'):
            # if relative == 'onset':
        rel_time = getattr(st, 'onset')
            # else:
            #     rel_time = st.starttime + (st.endtime - st.starttime) / 2
        winsig0 = rel_time - st.starttime + winsig[0]
        winsig1 = rel_time - st.starttime + winsig[1]
        winnoise0 = rel_time - st.starttime + winnoise[0]
        winnoise1 = rel_time - st.starttime + winnoise[1]
        #
        t = np.arange(self.stats.npts) * 1. / st.sampling_rate
        datasig = self_copy.data[(t >= winsig0) * (t <= winsig1)]
        datanoise = self_copy.data[(t >= winnoise0) * (t <= winnoise1)]
        # ipshell()
        try:
            st.signoise = max(abs(datasig)) / max(abs(datanoise))
            return st.signoise
        except:
            st.signoise = 0
            return st.signoise
########
client = Client("IRIS")
time_start_a=UTCDateTime('2010-03-01')
time_start_b=UTCDateTime('2012-05-01')
time_end_b=UTCDateTime('2014-01-01')
time_end_a=UTCDateTime('2013-01-01')

# inventory = client.get_stations(network="TA",station='*',maxlatitude=36,startafter=time_start_a,startbefore=time_start_b,
    # endbefore=time_end_b, minlongitude=-95,maxlongitude=-84)#,level='stations')startbefore=time_start_b,
catfile='rf_profile_events.xml'
if not os.path.exists(catfile):
    catalog = client.get_events(starttime=time_start_a, endtime=time_end_b,minmagnitude=5.9,latitude=30,longitude=-90, minradius=30, maxradius=95)
    catalog.write(catfile, 'QUAKEML')
catalog = read_events(catfile)
print('# of events:',len(catalog))
for st in ['446A','346A','246A','146A','Z46A','Y46A','X46A','W46A']:
# for st in ['446A']:

    data = os.path.join('rf_data/{}'.format(st), '')
    # catfile = data+'rf_profile_events.xml'
    datafile = data + 'rf_profile_data.h5'
    rffile = data + 'rf_profile_rfs.h5'
    # station= data+'perm_st.txt'
    inventory = client.get_stations(network="TA",station=st,level='response')
    if not os.path.exists(data):
        os.makedirs(data)
    if not os.path.exists(datafile):
        stream = RFStream()
        with tqdm() as pbar:
            try:
                for s in iter_event_data(catalog, inventory, client.get_waveforms, pbar=pbar):#,**kw):
                    s.filter('bandpass', freqmin=0.03, freqmax=5)
                    if s[0].stats.sampling_rate > 19:
                        # SNR=signoise(s.select(component='Z')[0])
                        # if SNR >= 1.5:
                        s.detrend('linear')
                        s.detrend('demean').taper(max_percentage=0.05)
                        stream.extend(s)
            except:
                print('an event didnt work')
        stream.write(datafile, 'H5')
        print('Len of data per compo after SNR:',len(stream)/3)
