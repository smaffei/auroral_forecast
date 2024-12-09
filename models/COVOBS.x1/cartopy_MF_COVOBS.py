#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 15:02:23 2019
s
@author: Stefano Maffei

"""



import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/nfs/see-fs-01_users/earsmaf/Python/SpecialFunctions/SphericalHarmonics/')
from SH_library import DthSchmidtSH
from SH_library import DphSchmidtSH
from SH_library import SchmidtSH
from scipy.interpolate import BSpline
from mpl_toolkits.basemap import Basemap

# decide the times for which interpolate the splines
#times = np.linspace(1840,2020,100)
times = [2019]

nlats = 80
nlons = 100
"""
if we use cartopy (bugged!) we need to arrange lat and lon in the following way
lats = np.linspace(-np.pi / 2, np.pi / 2, nlats)
lons = np.linspace(0, 2 * np.pi, nlons)
"""
lats = np.linspace(np.pi/2,-np.pi/2,num=100)
lons = np.linspace(-np.pi,np.pi,num=100)
lons, lats = np.meshgrid(lons, lats)
theta = -lats+np.pi/2

# figures are saved here
figfolder = '/nfs/see-fs-01_users/earsmaf/geomag/models/COVOBS.x1/figures/'
# spline coefficients
fcovobs = open('/nfs/see-fs-01_users/earsmaf/geomag/models/COVOBS.x1/COV-OBS.x1-int', 'r')
header = fcovobs.readline()
degree, num_splines, spline_degree = [int(x) for x in fcovobs.readline().split()]
knots, coeffs = [],[]
while len(knots) < num_splines + spline_degree:
    x = [list(map(float, fcovobs.readline().split() ))]
    knots = np.append(knots,x)
while len(coeffs) < degree * (degree + 2) * num_splines:
    x = [list(map(float, fcovobs.readline().split() ))]
    coeffs = np.append(coeffs,x)
fcovobs.close()

coeffs = np.reshape(coeffs,((degree * (degree + 2)),num_splines),order ='F')

g10 = BSpline(knots,coeffs[0,:],spline_degree)
g11 = BSpline(knots,coeffs[1,:],spline_degree)
h11 = BSpline(knots,coeffs[2,:],spline_degree)


# colatidude and dipole intensity
m1 = np.sqrt(g10(times)**2 + g11(times)**2 + h11(times)**2)
colatitude = -(90.0 - 180/np.pi * np.arccos(np.divide(g10(times),m1)))


Br_a = np.zeros((len(times), theta.shape[0], theta.shape[1]) ) 
Br_c = np.zeros((len(times), theta.shape[0], theta.shape[1]) )
Bt_a = np.zeros((len(times), theta.shape[0], theta.shape[1]) )
Bp_a = np.zeros((len(times), theta.shape[0], theta.shape[1]) ) 

r_cmb = 3485.0e3
r_e   = 6371.0e3

l=1
m=0
sincos = 'c'
for ic in range(coeffs.shape[0]):

    dthSH = DthSchmidtSH(l,m,theta,lons,sincos)
    dphSH = DphSchmidtSH(l,m,theta,lons,sincos)
    SH    = SchmidtSH(l,m,theta,lons,sincos)
    coeff_B = BSpline(knots,coeffs[ic,:],spline_degree)
    for it in range(len(times)):
        Br_a[it] = Br_a[it] + (l+1) * (r_e/r_e)**(l+2) * coeff_B(times[it]) * SH 
    
        Br_c[it] = Br_c[it] + (l+1) * (r_e/r_cmb)**(l+2) * coeff_B(times[it]) * SH 
                    
        Bt_a[it] = Bt_a[it] - (r_e/r_e)**(l+2) * coeff_B(times[it]) * dthSH 
    
        Bp_a[it] = Bp_a[it] - (r_e/r_e)**(l+2) * coeff_B(times[it]) * np.divide(dphSH,np.sin(theta))

    # updte indices
    if sincos == 'c':
        if m==0:
            m=1
        elif m>0:
            sincos = 's'
    elif sincos == 's':
        if m<l:
            sincos='c'
            m = m+1
        elif m == l:
            l=l+1
            sincos = 'c'
            m=0   

F_a = np.sqrt(Br_a**2 + Bt_a**2 + Bp_a**2)

lats = np.rad2deg(lats)
lons = np.rad2deg(lons)

########
########
# plots
########
########

#sinlge plot
files = []

it =-1

# geomagnetic dipole
dipole_colatitude = -(90.0 - 180/np.pi * np.arccos(np.divide(g10(times[it]),m1[it])))
dipole_lon        =  180/np.pi * np.arctan(np.divide( g11(times[it]),h11(times[it]) ) ) 

# find magnetic poles
inclination = np.arctan(np.divide(-Br_a[it],np.sqrt(Bt_a[it]**2+Bp_a[it]**2)))
idxN = np.nanargmin(np.abs(inclination - np.pi/2))
idxS = np.nanargmin(np.abs(inclination + np.pi/2))

N_pole_colatitude = lats.flat[idxN]
N_pole_lon        = lons.flat[idxN]
S_pole_colatitude = lats.flat[idxS]
S_pole_lon        = lons.flat[idxS]

# find minimum intensity (SAA)
idxSAA = np.nanargmin(F_a[it])
SAA_colatitude = lats.flat[idxSAA]
SAA_lon        = lons.flat[idxSAA]


####################
# Br at cmb
####################

plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)

# draw lat/lon grid lines every 30 degrees.
#map1.drawmeridians(np.arange(0,360,30),linewidth=1)
#map1.drawparallels(np.arange(-90,90,30),linewidth=1)
bounds = [-1000, -750, -500, -250, 0, 250, 500, 750, 1000]

xs, ys = map1(lons, lats)
cf = map1.contourf(xs,ys,Br_c[it]/1000, cmap='coolwarm',
                   vmax=bounds[-1], vmin=bounds[0], levels =bounds, extend ='both')
cs = map1.contour(xs,ys,Br_c[it]/1000, 23,
                  colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')
'''
# check this again, it does not work....
# geomagnetic poles location
plt.scatter([dipole_lon,dipole_lon-180],[dipole_colatitude,-dipole_colatitude],
           s=100, c='w',edgecolor='k',alpha=1, zorder=3,marker='*')
# magnetic poles location
plt.scatter([N_pole_lon,S_pole_lon],[N_pole_colatitude,S_pole_colatitude],
           s=100, c='r',edgecolor='k',alpha=1, zorder=3,marker='*')
'''
plt.colorbar(cf,fraction=0.03, pad=0.04)
plt.title(times[it], y=1.08)
plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)

# saving the plot
fname = figfolder +'Br_CMB_%07d.png' % it
plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
files.append(fname)



####################
# plot to check things
####################

fig = plt.figure(figsize=(11, 6))
plt.contourf(lons, lats, Br_c[it,:,:]/1000, # breaks down with contourlevels >=24
            cmap='coolwarm')


####################
# surface intensity
####################

plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)

# draw lat/lon grid lines every 30 degrees.
#map1.drawmeridians(np.arange(0,360,30),linewidth=1)
#map1.drawparallels(np.arange(-90,90,30),linewidth=1)
bounds = [18000, 24000, 30000, 36000, 42000, 48000, 54000, 60000, 66000, 72000]

xs, ys = map1(lons, lats)
cf = map1.contourf(xs,ys,F_a[it], cmap='jet',
                   vmax=bounds[-1], vmin=bounds[0], levels =bounds, extend = 'both')
cs = map1.contour(xs,ys,F_a[it], 23,
                  colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')
# find border of SAA
idxBorder = np.argmin(np.abs(cs.levels -32000))
cs.collections[idxBorder].set_color('white')
cs.collections[idxBorder].set_linewidth(3)
cs.collections[idxBorder].set_alpha(1)

'''
# geomagnetic poles location
plt.scatter([dipole_lon,dipole_lon-180],[dipole_colatitude,-dipole_colatitude],
           s=100, c='w',edgecolor='k',alpha=1, zorder=3,marker='*',
           transform=ccrs.PlateCarree())
# magnetic poles location
plt.scatter([N_pole_lon,S_pole_lon],[N_pole_colatitude,S_pole_colatitude],
           s=100, c='r',edgecolor='k',alpha=1, zorder=3,marker='*',
           transform=ccrs.PlateCarree())
# SAA location
plt.scatter([SAA_lon],[SAA_colatitude],
           s=100, c='r',edgecolor='k',alpha=1, zorder=3,marker='o',
           transform=ccrs.PlateCarree())
'''    
plt.colorbar(cf,fraction=0.03, pad=0.04)
plt.title(times[it], y=1.08)
plt.text(420e5, 9e6, r'$nT$' , fontsize=16)

# saving the plot
fname = figfolder +'F_surface_%07d.png' % it
plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
files.append(fname)
