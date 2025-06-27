import os
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from math import cos, sin, asin, sqrt, radians
from matplotlib import cm
import numpy as np
from scipy.optimize import curve_fit
#
def compaction_model(z, v0, a, p):
    return v0 + a * z**p

##
thickness_dict = {}
for line in open("sedi_sc_stations_all.txt", "r"):
    parts = line.split()
    station = parts[1]
    thickness = round(float(parts[4])/1000,3)
    thickness_dict[station] = thickness

print(thickness_dict)

output = []
for line in open("rf_sc_dt_DT_r0.txt", "r"):
    parts = line.split()
    net, sta = parts[0], parts[1]
    floats = list(map(float, parts[2:]))  # convert remaining to float
    thick = thickness_dict.get(sta, float('nan'))  # NaN if not found
    floats.append(thick)
    vs=round(2 * floats[5] / floats[3], 3)
    floats.append(vs) ## Vs= 2*D / DT
    vp= round(thick * vs /(thick - vs*floats[2]),3) ## Vp= D*Vs / d - dt Vs
    floats.append(vp)
    floats.append(round(vp/vs,2))
    output.append([net, sta] + floats)

# format network station lat long dt Dt r0 basement_depth Vs(from DT) Vp Vp/Vs


### next bit fits a compaction model to Vs and basement depth..
vs= np.array([row[8] for row in output])
vp= np.array([row[9] for row in output])
vp_vs= np.array([row[10] for row in output])
depth= np.array([row[7] for row in output])


# conclusion from 27th June.
# Vs makes sense. Vp doesn't. Why? Vp/vs is screwed as a result.
sys.exit()
# Initial guess for [v0, a, p]
initial_guess = [0.05, .1, 0.6]
bounds = ([0.001, 0, 0.1], [1.0, 2, 1.0])

# Fit the model
params, cov = curve_fit(compaction_model, depth, vs,bounds=bounds, p0=initial_guess,maxfev=10000)

# Unpack parameters
v0, a, p = params
print(f"Fitted parameters:\nv0 = {v0:.4f}, a = {a:.6f}, p = {p:.4f}")

# Plot for visualization (optional)
depth_fit = np.linspace(min(depth), max(depth), 100)
vs_fit = compaction_model(depth_fit, *params)

plt.scatter(depth, vs, label='Data', color='cadetblue')
plt.plot(depth_fit, vs_fit, label='Fit: v0 + a*z^p', color='salmon')
plt.xlabel('Thickness z (km)')
plt.ylabel('Vs (km/s)')
plt.legend()
plt.title('Compaction Curve Fit')
plt.grid(True)
plt.text(0.9, 0.15, f'v₀ = {v0:.3f}, a = {a:.3f}, p = {p:.3f}', transform=plt.gca().transAxes, ha='right', va='bottom')
# plt.show()
plt.savefig('vs_compaction.png',dpi=600,bbox_inches='tight', pad_inches=0.1)
