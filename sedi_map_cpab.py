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
import rasterio
from rasterio.warp import reproject, calculate_default_transform, Resampling
from matplotlib.colors import ListedColormap
from cmcrameri import cm

#####################
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
ax1.set_extent([-101,-73,25.2,41], crs=ccrs.PlateCarree())
#
# with rasterio.open(moho_file) as src:
#     image = src.read(1)
#     extent = [src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]

image_clean = np.ma.masked_where(dst_array < 0.0001, dst_array)
extent = [
    transform[2], transform[2] + transform[0] * width,
    transform[5] + transform[4] * height, transform[5]
]

img=ax1.imshow(image_clean, origin='upper', extent=extent,transform=ccrs.PlateCarree(), cmap=cmapA, vmin=0,vmax=8000)
ax1.coastlines(resolution="10m",color="#111111", linewidth=0.5, zorder=99)
ax1.add_feature(cartopy.feature.STATES.with_scale('10m'),linewidth=0.5,edgecolor='gray')

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=.8, color='gray', alpha=0.5, linestyle='--',rotate_labels=False,
        x_inline=False, y_inline=False)


gl.xlocator = mticker.FixedLocator([-100,-90,-80,-70])
gl.ylocator = mticker.FixedLocator([25,30,35,40])

# gl.xlines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 14}
gl.ylabel_style = {'size': 14}
cbar=plt.colorbar(img,ax=ax1,shrink=0.5, extend='max', drawedges=False, pad=0.02 )
cbar.set_label("Sediment thickness (m)",fontsize=12)
# fig.savefig('sedi_cpa.png', dpi=400,bbox_inches='tight', pad_inches=0.1)
plt.show()
