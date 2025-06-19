import numpy as np
import matplotlib.pyplot as plt
import obspy
from obspy import read, Stream, UTCDateTime,read_events
import os
import sys
import rf
from rf import read_rf, rfstats,RFStream
from rf import iter_event_data, IterMultipleComponents
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import shapely.geometry as sgeom
from shapely.ops import unary_union
import cartopy.feature as cfeature
from matplotlib.path import Path
from shapely.prepared import prep
import rasterio
from rasterio.warp import reproject, calculate_default_transform, Resampling
from matplotlib.colors import ListedColormap
from cmcrameri import cm
import xarray as xr
import pandas as pd
#####################
def mask_ocean(dst_array,image_clean):

    land_geom = cfeature.NaturalEarthFeature('physical', 'land', '50m').geometries()
    land_union = unary_union(list(land_geom))
    land_prepared = prep(land_union)

    # Create mask using lon/lat mesh
    lon_grid, lat_grid = np.meshgrid(
        np.linspace(extent[0], extent[1], dst_array.shape[1]),
        np.linspace(extent[2], extent[3], dst_array.shape[0])
    )

    # Flatten and check each point
    mask = np.array([
        not land_prepared.contains(sgeom.Point(x, y))
        for x, y in zip(lon_grid.ravel(), lat_grid.ravel())
    ]).reshape(dst_array.shape)

    # Apply land mask (mask ocean)
    masked_array = np.ma.masked_where(mask, image_clean)

    return masked_array,mask

moho_file='CPAB_sedi_boyd/CPABasementDepth_v220517.tif'

with rasterio.open(moho_file) as src:
    # Set target CRS
    dst_crs = 'EPSG:4326'

    # Calculate transform and shape for the output image
    transform, width, height = calculate_default_transform(
        src.crs, dst_crs, src.width, src.height, *src.bounds)

    # Create destination array
    dst_array = np.empty((height, width), dtype=np.float32)

    # Reproject
    reproject(
        source=src.read(1),
        destination=dst_array,
        src_transform=src.transform,
        src_crs=src.crs,
        dst_transform=transform,
        dst_crs=dst_crs,
        resampling=Resampling.bilinear
    )


proj = ccrs.Stereographic(central_longitude=-90, central_latitude=90, true_scale_latitude=33)
# plt.clf()
cmap = cm.batlow
# Transparent colours
###
colA = cmap(np.arange(cmap.N))
colA[:,-1] = 0.35 + 0.5 * np.linspace(-1.0, 1.0, cmap.N)**2.0
#adjusts the opacity based on the quadratic curve, setting values between 0.35 and 0.85.

cmapA = ListedColormap(colA)

# sys.exit()
plt.ion()
fig = plt.figure(figsize=(15, 8), facecolor=None)
ax1 = plt.subplot(1, 1, 1, projection=proj)
ax1.set_extent([-83.5,-78.7,32.2,35.2], crs=ccrs.PlateCarree())
#
image_clean = np.ma.masked_where(dst_array < 0.0001, dst_array)
extent = [
    transform[2], transform[2] + transform[0] * width,
    transform[5] + transform[4] * height, transform[5]
]

img=ax1.imshow(image_clean, origin='upper', extent=extent,transform=ccrs.PlateCarree(), cmap=cmapA, vmin=0,vmax=900)
ax1.coastlines(resolution="10m",color="#111111", linewidth=0.5, zorder=99)
ax1.add_feature(cfeature.OCEAN.with_scale('10m'),alpha=1,facecolor='xkcd:light grey',zorder=90)

ax1.add_feature(cartopy.feature.STATES.with_scale('10m'),linewidth=0.5,edgecolor='dimgray')

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=.8, color='gray', alpha=0.5, linestyle='--',rotate_labels=False,
        x_inline=False, y_inline=False)


gl.xlocator = mticker.FixedLocator([-83,-82,-81,-80,-79])
gl.ylocator = mticker.FixedLocator([33,34,35])

# gl.xlines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}
cbar=plt.colorbar(img,ax=ax1,shrink=0.5, extend='max', drawedges=False, pad=0.02 )
cbar.set_label("Sediment thickness (m)",fontsize=12)
# fig.savefig('sedi_sc.png', dpi=400,bbox_inches='tight', pad_inches=0.1)


### plot sc stations and get sedi thickness
station_names = ["E30", "Z57A",
"E31", "Z58A",
"TEEBA", "Z59A",
"155A", "Y58A",
"156A", "Y59A",
"157A", "Z55A",
"158A", "Z56A"
]
###
# df['Longitude'] = df['Longitude'].round(2)
# df['Latitude'] = df['Latitude'].round(2)
df = pd.read_csv('CPAB_sedi_boyd/CPABasementDepth_v220517_p01.csv', delim_whitespace=True)
df.columns = df.columns.str.strip().str.replace(',', '', regex=False)

with open('CPAB_sedi_boyd/CPABasementDepth_v220517_p01.csv', 'r') as f:
    lines = f.readlines()

cleaned_lines = [line.replace(',', '').strip() for line in lines if line.strip()]
split_lines = [line.split() for line in cleaned_lines]

df = pd.DataFrame(split_lines, columns=['Longitude', 'Latitude', 'Thickness'])

df['Longitude'] = pd.to_numeric(df['Longitude'], errors='coerce')
df['Latitude'] = pd.to_numeric(df['Latitude'], errors='coerce')
df['Thickness'] = pd.to_numeric(df['Thickness'], errors='coerce')

# Drop rows with NaNs
df = df.dropna(subset=['Longitude', 'Latitude', 'Thickness'])

# next command takes time!
thickness_lookup = {(row['Longitude'], row['Latitude']): row['Thickness']
    for _, row in df.iterrows()
}

station_thickness = {}

for line in open('stations_SC.txt','r'):
    line=line.split()
    if line[0] in station_names:
        ax1.plot(float(line[3]), float(line[2]), marker='^',markersize=9, linestyle='None', markerfacecolor='none', markeredgecolor='navy', transform=ccrs.PlateCarree(),zorder=99)
        ax1.text(float(line[3]) + 0.05, float(line[2])-.04, line[0], fontsize=10, color='navy', transform=ccrs.PlateCarree(), zorder=99)
        station_lon = round(float(line[3]), 2)
        station_lat = round(float(line[2]), 2)

        thickness = thickness_lookup.get((station_lon, station_lat), None)
        station_thickness[line[0]] = thickness

        # Interpolate thickness value
        # thickness = thickness_ds.interp(lon=float(line[3]), lat=float(line[2]), method='nearest').item()

        print(f"Station {line[0]}: Thickness = {thickness:.2f} meters")
        thick=int(thickness)
        ax1.text(float(line[3]) - 0.4, float(line[2]), "{}(m)".format(thick), fontsize=9, color='maroon', transform=ccrs.PlateCarree(), zorder=99)


plt.show()
