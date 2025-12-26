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
import cartopy.feature as cfeature

#######################################


client = Client("IRIS")
# Client.get_stations(starttime=None, endtime=None, startbefore=None, startafter=None,
# endbefore=None, endafter=None, network=None, station=None, location=None,
# channel=None, minlatitude=None, maxlatitude=None, minlongitude=None,
# maxlongitude=None, latitude=None, longitude=None, minradius=None,
# maxradius=None, level=None, includerestricted=None, includeavailability=None,
#  updatedafter=None, matchtimeseries=None, filename=None, format=None, **kwargs)
time_start_a=UTCDateTime('2010-03-01')
time_start_b=UTCDateTime('2012-05-01')

time_end_b=UTCDateTime('2014-01-01')
time_end_a=UTCDateTime('2013-01-01')

st_list=('446A','346A','246A','146A','Z46A','Y46A','X46A','W46A')
inventory = client.get_stations(network="TA",station='*',maxlatitude=36,startafter=time_start_a,startbefore=time_start_b,
endbefore=time_end_b, minlongitude=-95,maxlongitude=-84)#,level='stations')startbefore=time_start_b,
# inventory.plot()
start_times=[]
end_times=[]
for network in inventory:
    # for station in network:
    for station in network:
        # sys.exit()
        if station.code in st_list:
            start_times.append(station.start_date)
            end_times.append(station.end_date)
            print(f"Station: {station.code}, Start: {station.start_date}, End: {station.end_date}")
            print('#days active:{:.1f}\n'.format((station.end_date-station.start_date)/86400))
###
sys.exit()
plt.rcParams.update({'font.size': 5})
fig = inventory.plot(projection='local',resolution='h',marker='^',
    fillstyle='none',size=15,label=True,color='darkcyan',continent_fill_color='ivory')

gl = ax.gridlines(draw_labels=True,linewidth=0.3,color='gray',alpha=0.5,linestyle='--')
gl.top_labels = False
gl.right_labels = False
# plt.grid(True)
plt.show()
# plt.savefig('TA_2012_active.png', dpi=300,bbox_inches='tight', pad_inches=0.1)
######
###
# ax = fig.axes[0]   # get the Cartopy GeoAxes
# states = cfeature.NaturalEarthFeature(
#     category='cultural',
#     name='admin_1_states_provinces_lines',
#     scale='50m',
#     facecolor='none')
# ax.add_feature(states, edgecolor='0.4', linewidth=0.3)
