import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client
from obspy.core.trace import Trace
import numpy as np
import scipy.fft as fft
from scipy import signal
from scipy.signal import hilbert, correlate
import datetime
import os
import sys
import rf
from rf import read_rf, rfstats,RFStream
from rf import iter_event_data, IterMultipleComponents
from tqdm import tqdm
from scipy import signal
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from math import cos, sin, asin, sqrt, radians
from matplotlib import cm
import func_plotrf as mm
from importlib import reload
###
def max_p_onset(self):
    """
    Determine max amplitude of r-rf around P-onset.
    """

    onset=self.stats.onset
    self_copy = self.copy()
    self_copy.trim(onset-5,onset+ 35)
    self_copy.taper(max_percentage=0.05)
    st = self_copy.stats
    rel_time = getattr(st, 'onset')

    winP_st = rel_time - st.starttime - 0.1    # P onset -1
    winP_end = rel_time - st.starttime + 1.8  # P onset +1
    t = np.arange(st.npts) * 1. / st.sampling_rate
    max_abs=max(abs(self_copy.data[(t >= winP_st) * (t <= winP_end)])) # gets max amp around P (-1,1 sec)
    max_=max(self_copy.data[(t >= winP_st) * (t <= winP_end)])
    self.stats.max_abs_P=max_abs
    self.stats.max_P=max_
    return max_abs,max_

#####
def apply_reverberation_filter(cha_data):
    result_stream = []
    # station_file=open("P_delay_5G_noNan.txt","r")
    # for line1 in station_file:
    #     RF_tps=float(line1.split()[2])
    #     name= line1.split()[3]
    #     if name==cha_data[0].stats.station:
    #         cha_data[0].stats.update({'pdelay':RF_tps})
    # station_file.close()
    auto_c=RFStream()
    for i, tr in enumerate(cha_data):
        lead_time = tr.stats.onset - tr.stats.starttime
        auto_c_temp=Trace()
        relative_time = tr.times() - lead_time
        mask = np.array((relative_time < 0) | (relative_time > 5.0))
        loc = np.argmax(np.ma.masked_array(tr.data, mask=mask))
        dt = relative_time[loc] # delay of P rf

        data = correlate(tr.data, tr.data, mode='full') #autocorrelation!
        data /= np.max(data)
        data = data[len(data) // 2:]
        auto_c_temp.data=data

        r0 = -(np.min(data)) #
        Dt = np.argmin(data) * 1. / tr.stats.sampling_rate # Dt in sec..

        # building our sinosouidal decay in time domain..
        tr_copy = tr.copy()
        resonance_filter = np.zeros(len(tr.data))
        resonance_filter[0] = 1
        resonance_filter[int(Dt * tr.stats.sampling_rate)] = r0 # read eq 5 Cunningham & Lekic (2019)
        tr_copy.data = np.convolve(tr_copy.data, resonance_filter, mode='full')
        tr_copy.data = tr_copy.data[:len(tr_copy.data) // 2 + 1]

        assert tr.data.shape == tr_copy.data.shape, 'Input/output length mismatch detected in ' \
                                                    'reverberation removal routine'

        tr_copy.stats.update({'t1_offset':dt, # delay of P rf
                              't2_offset':Dt - dt, # Dt is from autoc
                              't3_offset':Dt,
                              'r0':r0})
        result_stream.append(tr_copy)
        tr_ccc=tr_copy.copy()
        auto_c_temp.stats=tr_ccc.stats
        auto_c_temp.stats.starttime=tr_ccc.stats.onset # setting startime as onset for autoC
        t=auto_c_temp.stats.starttime
        auto_c.append(auto_c_temp.trim(t,t+30))
    # end for

    return rf.RFStream(result_stream),auto_c
###''
def quality_check(rfs):
    stream_rf = RFStream()
    for tr in rfs:
        # stream3c.taper(max_percentage=0.05)
        # stream3c.filter('bandpass', freqmin=0.1, freqmax=1)
        onset=tr.stats.onset
        tr.trim(onset-5,onset+ 45)
        tr.normalize()
        max_rf_r_abs=max(abs(tr.data))
        max_rf_r=max(tr.data)
        # cc_conv=stream3c.select(component='R')[0].stats.cc_conv
        tr.stats.max_rf_r_abs=max_rf_r_abs
        (max_abs_P,max_P)=max_p_onset(tr)
        # if max_abs_P==max_rf_r_abs:# and max_P >0:# and cc_conv>0.4: # if max of Direct P in R-rf is greater than rest of signal
        if max_P==max_rf_r: #and max_rf_r_abs <1:# and cc_conv>0.4 :# and max_P >0:# and cc_conv>0.4:
        # if max_rf_r_abs <1:# and cc_conv>0.4 :# and max_P >0:# and cc_conv>0.4:
            stream_rf.append(tr)
    return stream_rf
###
def plot_auto(auto_c,station):
    fig = plt.figure(figsize=(4,3))
    # plt.style.use('grayscale')
    ax1 = fig.add_axes([0.1, 0.07, 0.85, 0.85]) #[left, bottom, width, height]
    ax1.set_ylim(26, 95)
    ax1.set_xlim(0,10)
    ax1.set_yticks(np.linspace(30,90,3))
    for i,tr in enumerate(auto_c):
        dist=tr.stats.distance
        auto=tr
        time = np.arange(auto.stats.npts) * auto.stats.delta
        auto.data /= np.max(np.abs(auto.data))
        r0=auto.stats.r0
        Dt=auto.stats.t3_offset#*auto.stats.sampling_rate

        ax1.plot(time,5*auto.data + dist,  lw=0.15, color='black',alpha=.5)
        ax1.plot(Dt,5*-r0 + dist,marker='o', color='slateblue',markersize=1,alpha=.8)

        # ax1.fill_between(time, i+1, auto.data + i+1, lw=0.55,
                          # color=cm.viridis(auto.stats.tpdelay/1.32), where=(auto.data < 0),alpha=.8)
        # plt.setp(auto.stats.station,fontsize=3)
    ax1.xaxis.set_minor_locator(MultipleLocator(.25))
    ax1.yaxis.set_minor_locator(MultipleLocator(10))
    ax1.xaxis.set_major_locator(MultipleLocator(1))
    plt.setp(ax1.get_xticklabels(),fontsize=8)

    ax1.grid(which='major', axis='x',color='LightGrey', linestyle='--',linewidth=.45,alpha=.75)
    ax1.set_ylabel('Distance')
    ax1.set_xlabel('Time')
    ax1.set_title('AutoC for {}'.format(station),fontsize=11)
    plt.savefig('RFstacks/figs/autoc/{}_autoC.pdf'.format(station,'R'),bbox_inches='tight', pad_inches=0.2)
    # fig.show()
    plt.close()
####
def curate_rvr(rfs,dt_max,st):
    stream_rever_rmv_cur=RFStream()
    for tr in rfs:
        # tr=stream3c.select(component='R')
        # if tr.stats.max_P < 1:#used as another quality check for thick sediment stations
        dt_rf=tr.stats.t1_offset
        if 1.9*dt_rf < tr.stats.t3_offset < 3.1*dt_rf and dt_rf < dt_max:# and tr.stats.max_P < 1:#and tr.stats.distance :# change here for different stations based on Dt in auto_c plots
            stream_rever_rmv_cur.append(tr)
    print('length after removing large DT=', len(stream_rever_rmv_cur))
    kw = {'trim': (-2.5, 25),'fig_width':6, 'fillcolors': ('steelblue', 'gainsboro'), 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}
    stream_rever_rmv_cur.sort(['distance']).plot_rf(**kw)
    plt.savefig('RFstacks/figs/{}_q_rmvDt.pdf'.format(st),bbox_inches='tight', pad_inches=0.2)
    # plt.show()
    plt.close()
    stream_rever_rmv_cur.write('RFstacks/{}/{}_rvr_rmvDt'.format(st,st), 'H5')
    # return stream_rever_rmv_cur

#####
station_names = ["E30", "Z57A",
"E31", "Z58A",
"TEEBA", "Z59A",
"155A", "Y58A",
"156A", "Y59A",
"157A", "Z55A",
"158A", "Z56A"
]

# station_names = [ "Y58A", "Y59A"]
# station_names = [
    # "W58A", "X58A", "Y57A", "HAW"]
# st_name='E31'
for st_name in station_names:
    rfs_raw=read_rf('RFstacks/{}/radials_{}_updated.h5'.format(st_name,st_name),'H5')

    rfs_quality= quality_check(rfs_raw)
    print('No. of RF before/after =',len(rfs_raw),' / ',len(rfs_quality))
    ####

    kw = {'trim': (-5, 25),'fig_width':6, 'fillcolors': ('seagreen', 'white'), 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}
    rfs_quality.select(component='R').sort(['back_azimuth']).plot_rf(**kw)
    plt.savefig('RFstacks/figs/{}_q.pdf'.format(st_name),bbox_inches='tight', pad_inches=0.2)
    plt.close()
    # sys.exit()

    rfs_quality_rvr,auto_c= apply_reverberation_filter(rfs_quality)

    plot_auto(auto_c,st_name)

    # sys.exit()
    # rf_z57.normalize()
    kw = {'fig_width':6, 'fillcolors': ('maroon', 'white'), 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}
    # kw = {'trim': (-5, 25),'fig_width':6, 'fillcolors': ('seagreen', 'white'), 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}
    # kw = {'trim': (-5, 25),'fig_width':6, 'fillcolors': ('mistyrose', 'white'), 'trace_height': 0.1, 'show_vlines': 'True','scale':2.5}

    rfs_quality_rvr.select(component='R').sort(['back_azimuth']).plot_rf(**kw)
    plt.savefig('RFstacks/figs/{}_q_rvr.pdf'.format(st_name),bbox_inches='tight', pad_inches=0.2)
    plt.close()
    # teeba_q_rvr.select(component='R').sort(['back_azimuth']).plot_rf(**kw)

    # Z57: .5 < dt <1 ; 1.25< Dt < 1.75
    # teeba: .6 < dt <1.5 ; 1.75< Dt < 2.25
    curate_rvr(rfs_quality_rvr,1.75,st_name)

# curate_rvr(teeba_q_rvr,.6,1.5,1.75,2.25,'teeba')

# plt.show()

# plt.close()
# stream_rf.write(rffile_pdelay, 'H5')
