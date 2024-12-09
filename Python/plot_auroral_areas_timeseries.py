#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 13:32:59 2021

@author: stefanomaffei

To plot the time-series calculated via auroal_areas_timeseries.py
"""

import sys
import aacgmv2
import aacgm_functions # my module to implement aacgmv2
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as c
import cartopy.crs as ccrs
from cartopy.feature import ShapelyFeature
#import matplotlib.patches as mpatches
from shapely import geometry
import math
import os
from multiprocessing import Process, Pool, Queue
#from skimage import measure
import pandas as pd
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
import matplotlib.ticker

class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format
            
            
class Labeloffset():
    def __init__(self,  ax, label="", axis="y"):
        self.axis = {"y":ax.yaxis, "x":ax.xaxis}[axis]
        self.label=label
        ax.callbacks.connect(axis+'lim_changed', self.update)
        ax.figure.canvas.draw()
        self.update(None)

    def update(self, lim):
        fmt = self.axis.get_major_formatter()
        self.axis.offsetText.set_visible(False)
        self.axis.set_label_text(self.label + " "+ fmt.get_offset() )




# change plot fonts globally   
#plt.rcParams["font.family"] = "Times"
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
plt.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
plt.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.close('all')


# define linestiles:
igrf_LS = '--'
mpg_LS  = '-'
ipgp_LS = ':'

# base folder for script and figures
folder_base = './'
folder_out = folder_base+'aacgmv2_coords/'
folder_models = folder_base+'models/'
folder_figures = folder_base+'figures/'


# this should contain 1900-2020
csv_magmodel = folder_base+'areas_magmodel_lon_res0.5.csv'
df_magmodel = pd.read_csv(csv_magmodel, index_col=0)


# cap area

FS = 20



fig,ax = plt.subplots(figsize=(8,5))
ax.set_xlabel('year',fontsize=FS)
ax.set_ylabel('Area [km$^2$]', fontsize=FS)
#ax.set_ylim([77,81])
#ax.set_xlim([1900,2100])
ax.tick_params('x',labelsize=FS)

ax.plot(df_magmodel['year'],df_magmodel['cap_areaN[km^2]'],'-',color='tab:blue',label=r'North')

ax.plot(df_magmodel['year'],df_magmodel['cap_areaS[km^2]'],'-',color='tab:orange',label=r'South')


ax.legend(fontsize=15,loc='center left')

ax.tick_params('y',labelsize=FS)
#ax.yaxis.get_offset_text().set_fontsize(FS)
ax.yaxis.set_major_formatter(OOMFormatter(6, "%1.1f"))
lo = Labeloffset(ax, label="Area [km$^2$]", axis="y")

#ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
#ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
plt.title(r'Polar caps area',fontsize=FS)
plt.savefig('figures/cap_areas_magmodel.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()

# zone area
fig,ax = plt.subplots(figsize=(8,5))
ax.set_xlabel('year',fontsize=FS)
ax.set_ylabel('Area [km$^2$]', fontsize=FS)
#ax.set_ylim([77,81])
#ax.set_xlim([1900,2100])
ax.tick_params('x',labelsize=FS)

ax.plot(df_magmodel['year'],df_magmodel['zone_areaN[km^2]'],'-',color='tab:blue',label=r'North')

ax.plot(df_magmodel['year'],df_magmodel['zone_areaS[km^2]'],'-',color='tab:orange',label=r'South')
    
ax.legend(fontsize=15,loc='center left')

ax.tick_params('y',labelsize=FS)
#ax.yaxis.get_offset_text().set_fontsize(FS)
ax.yaxis.set_major_formatter(OOMFormatter(6, "%1.1f"))
lo = Labeloffset(ax, label="Area [km$^2$]", axis="y")
#ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
#ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
plt.title(r'Auroral zones area',fontsize=FS)
plt.savefig('figures/zone_areas_magmodel.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()




# centroid latitude
fig,ax = plt.subplots(figsize=(8,5))
ax.set_xlabel('year',fontsize=FS)
ax.set_ylabel('Latitude [deg]', fontsize=FS)
#ax.set_ylim([77,81])
#ax.set_xlim([1900,2100])
ax.tick_params('x',labelsize=FS)

ax.plot(df_magmodel['year'],df_magmodel['centroid_latsN[deg]'],'-',color='tab:blue',label=r'North')

ax.plot(df_magmodel['year'],df_magmodel['centroid_latsS[deg]'],'-',color='tab:orange',label=r'South')
    
ax.legend(fontsize=15,loc='center left')

ax.tick_params('y',labelsize=FS)
#ax.yaxis.get_offset_text().set_fontsize(FS)
#ax.yaxis.set_major_formatter(OOMFormatter(6, "%1.1f"))
lo = Labeloffset(ax, label="Latitude [deg]", axis="y")
#ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
#ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
plt.title(r'Auroral zones centroid latitude',fontsize=FS)
plt.savefig('figures/zone_centroid_lats_magmodel.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()


# centroid abs latitude
fig,ax = plt.subplots(figsize=(8,5))
ax.set_xlabel('year',fontsize=FS)
ax.set_ylabel('Latitude [deg]', fontsize=FS)
#ax.set_ylim([77,81])
#ax.set_xlim([1900,2100])
ax.tick_params('x',labelsize=FS)

ax.plot(df_magmodel['year'],abs(df_magmodel['centroid_latsN[deg]']),'-',color='tab:blue',label=r'North')

ax.plot(df_magmodel['year'],abs(df_magmodel['centroid_latsS[deg]']),'-',color='tab:orange',label=r'South')
    
ax.legend(fontsize=15,loc='best')

ax.tick_params('y',labelsize=FS)
#ax.yaxis.get_offset_text().set_fontsize(FS)
#ax.yaxis.set_major_formatter(OOMFormatter(6, "%1.1f"))
lo = Labeloffset(ax, label="Latitude [deg]", axis="y")
#ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
#ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
plt.title(r'Auroral zones centroid latitude',fontsize=FS)
plt.savefig('figures/zone_centroid_abs_lats_magmodel.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()


# centroid lopngitude
fig,ax = plt.subplots(figsize=(8,5))
ax.set_xlabel('year',fontsize=FS)
ax.set_ylabel('Longitude [deg]', fontsize=FS)
#ax.set_ylim([77,81])
#ax.set_xlim([1900,2100])
ax.tick_params('x',labelsize=FS)

ax.plot(df_magmodel['year'],df_magmodel['centroid_lonsN[deg]'],'-',color='tab:blue',label=r'North')

ax.plot(df_magmodel['year'],df_magmodel['centroid_lonsS[deg]'],'-',color='tab:orange',label=r'South')
    
ax.legend(fontsize=15,loc='center left')

ax.tick_params('y',labelsize=FS)
#ax.yaxis.get_offset_text().set_fontsize(FS)
#ax.yaxis.set_major_formatter(OOMFormatter(6, "%1.1f"))
lo = Labeloffset(ax, label="Longitude [deg]", axis="y")
#ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
#ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
plt.title(r'Auroral zones centroid longitude',fontsize=FS)
plt.savefig('figures/zone_centroid_lons_magmodel.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()
