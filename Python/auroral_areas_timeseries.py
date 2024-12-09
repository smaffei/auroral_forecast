#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 16:29:48 2021

@author: stefanomaffei

Calculate surface area evolution for the different forecasts
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

# parameters
r_cmb = 3485.0e3
r_a   = 6371.0e3

# nominal latitudinal bounds for auroral zones:
latNpol_target = 70
latNeq_target  = 65
latSpol_target = -70
latSeq_target  = -65
# first guess for the search for the bisection algorithm
latNpol_i = 45
latNeq_i = 89
latSpol_i = -89
latSeq_i  = -45
# maximum number of iterations for the bisection algorithm
max_iterations = 50
# define target resolution for latitudinal bisection algorithm
target_res = 1e-9
# define target resolution for finding intersections of auroal boundaries
target_res_zero = 1e-9
# parameters for bisection code
lon_res_b = np.deg2rad(0.5)
lons_b = np.linspace(0,2*np.pi,num=int(2*np.pi/lon_res_b))


# lats and lons for plotting the coordinates
lats_lin = np.linspace(-np.pi/2,np.pi/2,num=180)
lons_lin = np.linspace(0,2*np.pi,num=360)

lats, lons = np.meshgrid(lats_lin,lons_lin)
theta = -lats+np.pi/2

Lmax = 13

# change plot fonts globally   
#plt.rcParams["font.family"] = "Times"
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
plt.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
plt.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.close('all')


#####################################
# UNCOMMENT THE MODEL YOU WANT TO USE


# the default magmodel that came with the software: gufm1 + igrf13
model_folder = '../models/magmodel/'
model_file   = 'magmodel_1590-2020.txt' # in this case I only want to look at 1900+
model_name    = 'magmodel'
# temporal limits and offsets
t_offset = 0
t_in = 1900
t_fin = 2020



# base folder for script and figures
folder_base = './'
folder_out = folder_base+'aacgmv2_coords/'
folder_models = folder_base+'models/'
folder_figures = folder_base+'figures/'

zone_areasN = np.array([])
zone_areasS = np.array([])
cap_areasN = np.array([])
cap_areasS = np.array([])
centroid_latsN  = np.array([])
centroid_latsS  = np.array([])
centroid_lonsN  = np.array([])
centroid_lonsS  = np.array([])
years = np.array([])


def plus(a,b): return [x+y for x,y in zip(a,b)]
def minus(a,b): return [x-y for x,y in zip(a,b)]
def cross(a,b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def dot(a,b): return sum([x*y for x,y in zip(a,b)])
def length(v): return math.sqrt(dot(v,v))
def normalized(v): l = length(v); return [1,0,0] if l==0 else [x/l for x in v]
def addVectorTimesScalar(accumulator, vector, scalar):
    for i in range(len(accumulator)): accumulator[i] += vector[i] * scalar
def angleBetweenUnitVectors(a,b):
    # http://www.plunk.org/~hatch/rightway.php
    if dot(a,b) < 0:
        return math.pi - 2*math.asin(length(plus(a,b))/2.)
    else:
        return 2*math.asin(length(minus(a,b))/2.)

def sphericalPolygonMoment(verts):
    moment = [0.,0.,0.]
    for i in range(len(verts)):
        a = verts[i]
        b = verts[(i+1)%len(verts)]
        addVectorTimesScalar(moment, normalized(cross(a,b)),
                                     angleBetweenUnitVectors(a,b) / 2.)
    return moment

# main function    
if __name__ == "__main__":   
    
    def lonlat_degrees_to_xyz(lon_degrees,lat_degrees):
        lon = lon_degrees*(math.pi/180)
        lat = lat_degrees*(math.pi/180)
        coslat = math.cos(lat)
        return [coslat*math.cos(lon), coslat*math.sin(lon), math.sin(lat)]
    
    # output file name
    csv_name = folder_base + 'areas_'+model_name+'_lon_res'+str(np.rad2deg(lon_res_b))+'.csv'
    
    for year in range(int(t_in),int(t_fin)+1,5): 
        print(' ')        
        print('* calculation areas for year '+str(year)+' *' )
        print(' ')        


         # initialise lats arrays: n x 2 where n is range(lons), and the other dimension contains the polward and equatorward boundaries
        latsN = np.zeros((lons_b.shape[0],2))
        latsS = np.zeros((lons_b.shape[0],2))
                
        # if the calculation goes on, update this index
        # northern zone
        # polar edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(model_folder,model_file, year-t_offset, 1, 1, np.rad2deg(lons_b), latNpol_target, target_res, latNpol_i, latNeq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        # equatorial edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(model_folder,model_file, year-t_offset, 1, 1, np.rad2deg(lons_b), latNeq_target, target_res, latNpol_i, latNeq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        
        # southern zone
        # polar edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(model_folder,model_file, year-t_offset, 1, 1, np.rad2deg(lons_b), latSpol_target, target_res, latSpol_i, latSeq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        # equatorial edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(model_folder,model_file, year-t_offset, 1, 1, np.rad2deg(lons_b), latSeq_target, target_res, latSpol_i, latSeq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
                         
        # read the calculated values
        latsN[:,0] = np.loadtxt(folder_out+'bisection_'+str(latNpol_target)+'_lats_'+model_file)
        latsN[:,1] = np.loadtxt(folder_out+'bisection_'+str(latNeq_target)+'_lats_'+model_file)
        os.system('rm '+folder_out+'bisection_'+str(latNpol_target)+'_lats_'+model_file)
        os.system('rm '+folder_out+'bisection_'+str(latNeq_target)+'_lats_'+model_file)
        
        latsS[:,0] = np.loadtxt(folder_out+'bisection_'+str(latSpol_target)+'_lats_'+model_file)
        latsS[:,1] = np.loadtxt(folder_out+'bisection_'+str(latSeq_target)+'_lats_'+model_file)
        os.system('rm '+folder_out+'bisection_'+str(latSpol_target)+'_lats_'+model_file)
        os.system('rm '+folder_out+'bisection_'+str(latSeq_target)+'_lats_'+model_file)
        
        
        ###########################
        # calculate areas
        areaNpol_b = aacgm_functions.spherical_polygon_area(np.flip(latsN[:,0]),np.flip(np.rad2deg(lons_b)),r_a/1000)
        areaNeq_b  = aacgm_functions.spherical_polygon_area(np.flip(latsN[:,1]),np.flip(np.rad2deg(lons_b)),r_a/1000)
        areaN = areaNeq_b - areaNpol_b
        
        areaSpol_b = aacgm_functions.spherical_polygon_area(latsS[:,0],np.rad2deg(lons_b),r_a/1000)
        areaSeq_b  = aacgm_functions.spherical_polygon_area(latsS[:,1],np.rad2deg(lons_b),r_a/1000)
        areaS = areaSeq_b - areaSpol_b      
        
        
        ###########################
        # calculate centroids
        
        # polar cap centroids
        # north
        vertsN = [lonlat_degrees_to_xyz(lons_b[i]*180/np.pi,latsN[i,0]) for i in range(latsN[:,0].shape[0]) ]
        #verts=np.flipud(verts)
        
        momentN = sphericalPolygonMoment(vertsN)
        
        centroidN = normalized(momentN)
        latcN = math.asin(centroidN[2])*180/np.pi
        loncN = math.atan2(centroidN[1],centroidN[0])*180/np.pi
        
        # south
        vertsS = [lonlat_degrees_to_xyz(lons_b[i]*180/np.pi,latsS[i,0]) for i in range(latsS[:,0].shape[0]) ]
        vertsS=np.flipud(vertsS)
        
        momentS = sphericalPolygonMoment(vertsS)
        
        centroidS = normalized(momentS)
        latcS = math.asin(centroidS[2])*180/np.pi
        loncS = math.atan2(centroidS[1],centroidS[0])*180/np.pi
        
        
        
        # and append values to the arrays
        zone_areasN = np.append(zone_areasN,areaN)
        cap_areasN = np.append(cap_areasN,areaNpol_b)
        zone_areasS = np.append(zone_areasS,areaS)                         
        cap_areasS = np.append(cap_areasS,areaSpol_b)
        
        centroid_latsN  = np.append(centroid_latsN,latcN)
        centroid_latsS  = np.append(centroid_latsS,latcS)
        centroid_lonsN  = np.append(centroid_lonsN,loncN)
        centroid_lonsS  = np.append(centroid_lonsS,loncS)
        
        years = np.append(years,year)
    
    ##################
    # save some output
    print('')
    print('* saving output of calculation * ')
    data = {'year': years,
            'cap_areaN[km^2]':cap_areasN,
            'zone_areaN[km^2]':zone_areasN,
            'cap_areaS[km^2]':cap_areasS,
            'zone_areaS[km^2]':zone_areasS,
            'centroid_latsN[deg]':centroid_latsN,
            'centroid_latsS[deg]':centroid_latsS,
            'centroid_lonsN[deg]':centroid_lonsN,
            'centroid_lonsS[deg]':centroid_lonsS,
           }
    
    df = pd.DataFrame (data, columns = ['year',
                                        'cap_areaN[km^2]',
                                        'zone_areaN[km^2]',
                                        'cap_areaS[km^2]',
                                        'zone_areaS[km^2]',
                                        'centroid_latsN[deg]',
                                        'centroid_latsS[deg]',
                                        'centroid_lonsN[deg]',
                                        'centroid_lonsS[deg]'])

    df.to_csv(csv_name) 
    
else:
   print("Something went wrong")