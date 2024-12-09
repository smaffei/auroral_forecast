#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 15:02:23 2019
s
@author: Stefano Maffei

read in the CHAOS 6 model coefficients (stored in the shc format) and plot
"""



import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs
from mpl_toolkits.basemap import Basemap
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/nfs/see-fs-01_users/earsmaf/Python/SpecialFunctions/SphericalHarmonics/')
from SH_library import DthSchmidtSH
from SH_library import DphSchmidtSH
from SH_library import SchmidtSH
import os

# inputs for animation
yps = 10 # years per second of animation 


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

figname = 'CHAOS-6-x9_core'
folder = '/nfs/see-fs-01_users/earsmaf/geomag/models/CHAOS-6-x9/'
figfolder = '/nfs/see-fs-01_users/earsmaf/geomag/models/CHAOS-6-x9/figures/'
Bcoeffs = open(folder+'CHAOS-6-x9_core.shc', 'r')
junk1 = Bcoeffs.readline()
junk2 = Bcoeffs.readline()
junk3 = Bcoeffs.readline()
# read the header
header = Bcoeffs.readline()
x = [list(map(float, header.split() ))]
lmaxB = x[0][1]
# read times
times = Bcoeffs.readline()
x = [list(map(float, times.split() ))]
times = x[0][:]
line = []
coeffs_B = []
while True:
    line = Bcoeffs.readline()
    x = [list(map(float, line.split() ))]
    l = x[0][0]
    m = x[0][1]
    coeffs_B = np.append(coeffs_B,x[0][:])
    if m == -lmaxB: break # end of file

# Here is a matrix with rows the coefficients g/h evolution in time 
# and in columns all coefficients for a given time
# First two columns are the l and m coefficients. 
# For m negative, it's the sine coefficient
# For m positive or zero, it's the cosine coefficient
coeffs_B = np.reshape(coeffs_B,(len(coeffs_B)/len(x[0][:]), len(x[0][:]))) 
n_instants = len(x[0][2:-1])

Bcoeffs.close()

Br_a = np.zeros((n_instants, theta.shape[0], theta.shape[1]))
Br_c = np.zeros((n_instants, theta.shape[0], theta.shape[1]))
Bt_a = np.zeros((n_instants, theta.shape[0], theta.shape[1])) 
Bp_a = np.zeros((n_instants, theta.shape[0], theta.shape[1]))

r_cmb = 3485.0e3
r_e   = 6371.0e3

it = 0
# There is a problem at higher l with the core field. Too much power on small scales.
# is this a problem of how the .shc file was generated? 
# Is this litospheric field leaking in the low degrees?
# for now we can limit the upper l for the cmb field
#for ic in range(len(coeffs_B)): # up to l = 20
#for ic in range(223): # up to l = 14
for ic in range(194): # up to l = 13 (as used in the 2016 paper and in the Treatise)

    l, m = coeffs_B[ic][0:2]
    betalm = coeffs_B[ic][2:-1]
    
    if m<0:
        dthSH_s = DthSchmidtSH(l,-m,theta,lons,'s')
        dphSH_s = DphSchmidtSH(l,-m,theta,lons,'s')
        SH_s    = SchmidtSH(l,-m,theta,lons,'s')
        for it in range(n_instants):
            Br_a[it] = Br_a[it] + (l+1) * (r_e/r_e)**(l+2) * betalm[it] * SH_s
            Br_c[it] = Br_c[it] + (l+1) * (r_e/r_cmb)**(l+2) * betalm[it] * SH_s
            Bt_a[it] = Bt_a[it] - (r_e/r_e)**(l+2) * betalm[it] * dthSH_s
            Bp_a[it] = Bp_a[it] - (r_e/r_e)**(l+2) * betalm[it] * np.divide(dphSH_s,np.sin(theta))
    else:
        dthSH_c = DthSchmidtSH(l,m,theta,lons,'c')
        dphSH_c = DphSchmidtSH(l,m,theta,lons,'c')
        SH_c    = SchmidtSH(l,m,theta,lons,'c')
        for it in range(n_instants):
            Br_a[it] = Br_a[it] + (l+1) * (r_e/r_e)**(l+2) * betalm[it] * SH_c
            Br_c[it] = Br_c[it] + (l+1) * (r_e/r_cmb)**(l+2) * betalm[it] * SH_c
            Bt_a[it] = Bt_a[it] - (r_e/r_e)**(l+2) * betalm[it] * dthSH_c
            Bp_a[it] = Bp_a[it] - (r_e/r_e)**(l+2) * betalm[it] * np.divide(dphSH_c,np.sin(theta))
     
#    print("time = " +str(times[it]) + "; l = " + str(l) + "; m= " + str(m) + "; beta = " + str(betalm[it]))
                
F_a = np.sqrt(Br_a**2 + Bt_a**2 + Bp_a**2)

g10 = coeffs_B[0][2:-1]
g11 = coeffs_B[1][2:-1]
h11 = coeffs_B[2][2:-1]
m1  = np.sqrt(g10**2 + g11**2 + h11**2)

# geomagnetic dipole
dipole_colatitude = -(90.0 - 180/np.pi * np.arccos(np.divide(g10,m1)))
dipole_lon        =  180/np.pi * np.arctan(np.divide(g11,h11)) 



lats = np.rad2deg(lats)
lons = np.rad2deg(lons)

# find magnetic poles
inclination = np.arctan(np.divide(-Br_a,np.sqrt(Bt_a**2+Bp_a**2)))
idxN = np.zeros(n_instants)
idxS = np.zeros(n_instants)
idxSAA = np.zeros(n_instants)
N_pole_colatitude = np.zeros(n_instants)
N_pole_lon        = np.zeros(n_instants)
S_pole_colatitude = np.zeros(n_instants)
S_pole_lon        = np.zeros(n_instants)
SAA_colatitude = np.zeros(n_instants)
SAA_lon        = np.zeros(n_instants)

for it in range(n_instants):
    # find poles locations
    idxN[it] = np.nanargmin(np.abs(inclination[it] - np.pi/2))
    idxS[it] = np.nanargmin(np.abs(inclination[it] + np.pi/2))
    N_pole_colatitude[it] = lats.flat[int(idxN[it])]
    N_pole_lon[it]        = lons.flat[int(idxN[it])]
    S_pole_colatitude[it] = lats.flat[int(idxS[it])]
    S_pole_lon[it]        = lons.flat[int(idxS[it])]
    # find minimum intensity (SAA)
    idxSAA[it] = np.nanargmin(F_a[it])
    SAA_colatitude[it] = lats.flat[int(idxSAA[it])]
    SAA_lon[it]        = lons.flat[int(idxSAA[it])]



#################################################################################
# write to file
# (I m doing this beacue readspline.f does not work on the latest CHAOS release.
# should check from the original spline coefficients)
#################################################################################
it = 219

Bdat = open(folder+'CHAOS-6-x9_core_'+str(times[it])+'.dat', 'w+')

Bdat.write("CHAOS-6-x9_core Gauss coefficients for "+str(times[it])+"\n")
#lmax = -1 + np.sqrt(coeffs_B.shape[0] +1)
for ic in range(coeffs_B.shape[0]):
    l = coeffs_B[ic,0]
    m = coeffs_B[ic,1]
    if m ==0:
        g = coeffs_B[ic,it +2]
        h = 0
        Bdat.write('%3d%3d%18f%18f\n' % (l,0,g,h))
    elif m >0:
        g = coeffs_B[ic,it +2]
    elif m<0:
        h = coeffs_B[ic,it +2]
        Bdat.write('%3d%3d%18f%18f\n' % (l,-m,g,h))
        
Bdat.close()



#######
#######
# plots
#######
#######

######################
# cmb radial field
######################
# plot to check things

fig = plt.figure(figsize=(11, 6))
cf = plt.contourf(lons, lats, Br_a[it], 
            cmap='coolwarm')
plt.colorbar(cf)

fig = plt.figure(figsize=(11, 6))
cf = plt.contourf(lons, lats, Br_c[it], 
            cmap='coolwarm')
plt.colorbar(cf)


it = 219
plt.figure(figsize=(11,6))
#Mollweide projectionfrom scipy.interpolate import griddata
map1 = Basemap(projection='moll',lat_0=0,lon_0=0,resolution='l')

# draw coastlines, country boundaries, fill continents.
map1.drawcoastlines(linewidth=0.50)


xs, ys = map1(lons, lats)
cf = map1.contourf(xs,ys,Br_c[it]/1000, cmap='coolwarm')
cs = map1.contour(xs,ys,Br_c[it]/1000, 23,
                  colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')

plt.colorbar(cf,fraction=0.03, pad=0.04)
plt.title(times[it], y=1.08)
plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)


# Animation

files = []

for it in range(n_instants):
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
    
    plt.colorbar(cf,fraction=0.03, pad=0.04)
    plt.title(times[it], y=1.08)
    plt.text(420e5, 9e6, r'$\mu T$' , fontsize=16)

    # saving the plot
    fname = folder +'Br_CMB_%07d.png' % it
    plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
    files.append(fname)
    plt.close()
    
"""
# with Cartopy

files = []

it =99

fig = plt.figure(figsize=(11, 6))
#    plt.cla()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
#ax.gridlines(linewidth=1, alpha=0.5, linestyle='-')
Brmax = np.amax(abs(Br_a[it]))
cf = ax.contourf(lons, lats, Br_a[it], # breaks down with contourlevels >=24
            transform=ccrs.PlateCarree(),
            cmap='coolwarm')
ax.contour(lons, lats, Br_a[it],23, 
            transform=ccrs.PlateCarree(),
            colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')

plt.colorbar(cf,fraction=0.03, pad=0.04)
#cf.clim(vmin=-Brmax, vmax=Brmax)

plt.text(24300000, 0, r'$\mu T$' , fontsize=16)

ax.coastlines()
ax.set_global()

# geomagnetic poles location
plt.scatter([dipole_lon[it],dipole_lon[it]-180],[dipole_colatitude[it],-dipole_colatitude[it]],
           s=100, c='w',edgecolor='k',alpha=1, zorder=3,marker='*',
           transform=ccrs.PlateCarree())
# magnetic poles location
plt.scatter([N_pole_lon[it],S_pole_lon[it]],[N_pole_colatitude[it],S_pole_colatitude[it]],
           s=100, c='r',edgecolor='k',alpha=1, zorder=3,marker='*',
           transform=ccrs.PlateCarree())

plt.title(times[it], y=1.08)
plt.show()
fname = folder +'Br_CMB_%07d.png' % it
plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
files.append(fname)




"""




####################
# Surface intensity
###################

"""
# with cartopy

files = []

it = 99

fig = plt.figure(figsize=(11, 6))
#    plt.cla()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mollweide())
#ax.gridlines(linewidth=1, alpha=0.5, linestyle='-')
cf = ax.contourf(lons, lats, F_a[it], # breaks down with contourlevels >=24
            transform=ccrs.PlateCarree(),
            cmap='jet')
cs = ax.contour(lons, lats, F_a[it],23, 
            transform=ccrs.PlateCarree(),
            colors='k',linewidths=0.5,alpha = 0.3,linestyles = 'solid')

# find border of SAA
idxBorder = np.argmin(np.abs(cs.levels -32000))
cs.collections[idxBorder].set_color('white')
cs.collections[idxBorder].set_linewidth(3)
cs.collections[idxBorder].set_alpha(1)

plt.colorbar(cf,fraction=0.03, pad=0.04)
#cf.clim(vmin=-Brmax, vmax=Brmax)

plt.text(24300000, 0, r'$\mu T$' , fontsize=16)

ax.coastlines()
ax.set_global()


plt.title(times[it], y=1.08)
plt.show()
fname = folder +'F_surface_%07d.png' % it
plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
files.append(fname)

"""
# with basemap


# animation (can be modified to only plot one, of course)
files = []

for it in range(n_instants):
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
    
    plt.colorbar(cf,fraction=0.03, pad=0.04)
    plt.title(times[it], y=1.08)
    plt.text(420e5, 9e6, r'$nT$' , fontsize=16)

    # saving the plot
    fname = folder +'F_surface_%07d.png' % it
    plt.savefig(fname,  bbox_inches='tight',pad_inches=0.7,dpi=400)
    files.append(fname)
    plt.close()


dt = times[1]-times[0] #assuming a constant timestep
fps = int(yps/dt)
fps = 20
os.system('ffmpeg -r ' +str(fps)+' -i '+folder+'F_surface_%07d.png -vcodec libx264 -y -an ' +folder+ 'F_surface.mp4 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"')
#os.system('ffmpeg -r ' +str(fps)+' -i -pattern_type glob '+folder+'F_surface_%07d.png -c:v libx264 -pix_fmt yuv420p -s 1920x1080 ' +folder+ 'F_surface_mac.mp4')