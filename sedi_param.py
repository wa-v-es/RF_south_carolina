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

#####
def get_stack(rfs):
    rfs.trim2(-1,10,reftime='onset')
    stack_st=rfs.stack()
    if len(stack_st) >1:
        stack_st=Stream(stack_st[0])

    npts = stack_st[0].stats.npts
    delta = stack_st[0].stats.delta
    times = np.arange(npts) * delta
    trace = -1*stack_st[0].data / np.max(np.abs(stack_st[0].data)) * 0.6  # scale for visual spacing

    return trace,times

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
st_name='E31'

fig = plt.figure(figsize=(15, 12)) # for PA
ax1 = fig.add_axes([0.07, 0.2, 0.83, 0.2]) #[left, bottom, width, height]
ax2 = fig.add_axes([0.07, 0.4,0.83, 0.2])
ax3 = fig.add_axes([0.07, 0.6, 0.83, 0.2])
ax4 = fig.add_axes([0.07, 0.8, 0.83, 0.15])

ax2.sharex(ax1)
ax3.sharex(ax2)
ax4.sharex(ax3)

# fig = plt.figure(figsize=(6, 9)) #for SA

# ax1 = fig.add_axes([0.07, 0.2, 0.83, 0.35]) #[left, bottom, width, height]
plt.ion()
# plt.style.use('fivethirtyeight')
plt.style.use('seaborn-v0_8-poster')


station_dt_means = []
for st_name in station_names:
    rfs = read_rf('RFstacks/{}/{}_rvr_rmvDt.h5'.format(st_name,st_name),'H5')
    dt_vals = [tr.stats.t1_offset for tr in rfs]
    mean_dt = np.mean(dt_vals)
    station_dt_means.append((st_name, mean_dt))

# Sort by mean dt_st
sorted_stations = sorted(station_dt_means, key=lambda x: x[1])  # ascending
sorted_station_names = [s[0] for s in sorted_stations]

i=0
tick_labels=[]
for i, st_name in enumerate(sorted_station_names, start=1):

    rfs_raw=read_rf('RFstacks/{}/radials_{}_updated.h5'.format(st_name,st_name),'H5')
    rfs_rvr_corr=read_rf('RFstacks/{}/{}_rvr_rmvDt.h5'.format(st_name,st_name),'H5')

    raw_tr,raw_times=get_stack(rfs_raw)
    stack_corr,times_corr=get_stack(rfs_rvr_corr)

    # Plot vertically: time on y-axis, offset on x-axis
    # ax4.plot(raw_tr + i, raw_times, color='maroon', linewidth=1, zorder=5)
    ax4.plot(stack_corr + i, times_corr, color='royalblue', linewidth=1.1, zorder=5)

    # ax4.fill_betweenx(times, i, trace + i, where=(trace < 0), facecolor='lightcoral', alpha=0.8, linewidth=0.0, zorder=5)

    print('No. of RF =',len(rfs_rvr_corr),'for station',st_name)
    ####
    r0_st, dt_st, Dt_st = [], [], []
    for tr in rfs_rvr_corr:
        ax1.scatter(i,tr.stats.t1_offset,s=40,marker='o',facecolor='cadetblue',alpha=.5,linewidth=.75,zorder=10)
        ax2.scatter(i,tr.stats.t3_offset,s=40,marker='o',facecolor='indianred',alpha=.5,linewidth=.75,zorder=10)
        ax3.scatter(i,tr.stats.r0,s=40,marker='o',facecolor='seagreen',alpha=.5,linewidth=.75,zorder=10)

        r0_st.append(tr.stats.r0)
        dt_st.append(tr.stats.t1_offset)
        Dt_st.append(tr.stats.t3_offset)

    ax1.errorbar(i,np.mean(dt_st), yerr=0,ecolor='brown',marker='d', alpha=.95,markerfacecolor='khaki', markeredgecolor='dimgray',markersize=9,linestyle='none',zorder=90)
    ax2.errorbar(i,np.mean(Dt_st), yerr=0,ecolor='brown',marker='d', alpha=.95,markerfacecolor='khaki', markeredgecolor='dimgray',markersize=9,linestyle='none',zorder=90)
    ax3.errorbar(i,np.mean(r0_st), yerr=0,ecolor='brown',marker='d', alpha=.95,markerfacecolor='khaki', markeredgecolor='dimgray',markersize=9,linestyle='none',zorder=90)



ax1.set_xticks(range(1, len(sorted_station_names) + 1))
ax1.set_xticklabels(sorted_station_names,rotation=45, ha='center')
ax2.set_xticklabels(sorted_station_names,rotation=45, ha='center',alpha=0)
ax3.set_xticklabels(sorted_station_names,rotation=45, ha='center',alpha=0)
ax4.set_xticklabels(sorted_station_names,rotation=45, ha='center',alpha=0)


ax1.grid(color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
ax2.grid(color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
ax3.grid(color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)
ax4.grid(color='dimgrey', linestyle='--',linewidth=.65,alpha=.75)



ax1.set_xlabel('Seismic stations',labelpad=15)
ax1.set_ylabel('$\delta t$ (s)')
ax2.set_ylabel('$\Delta t$ (s)')
ax3.set_ylabel('$r_{0}$')
ax4.set_ylabel('time (s)')

# ax2.set_xticks([])
# ax3.set_xticks([])
# plt.show()

# plt.close()
plt.savefig('RFstacks/sedi_param_stacks_corr.png',dpi=600,bbox_inches='tight', pad_inches=0.1)
