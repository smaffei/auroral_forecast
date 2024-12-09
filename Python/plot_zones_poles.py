#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 13:58:42 2022

@author: stefanomaffei

https://stackoverflow.com/questions/19897187#answer-38201499

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

#sys.path.append('/Users/stefanomaffei/Python//SpecialFunctions/SphericalHarmonics/')
import SH_library

# parameters
r_cmb = 3485.0e3
r_a   = 6371.0e3


# the list of locations I want the plots of
# Quebec City:
lat_QC = 46.83724124586545
lon_QC = -71.13847494276122    
# Leeds
lat_Leeds = 53.810380928219175 # in rad: 0.9391683189497871
lon_Leeds = -1.5494718924623057 # in rad : -0.027043386190574743
# Dunedin
lat_Dun = -45.813949257824476
lon_Dun = 170.43247246043904
# Edmonton
lat_Ed = 53.55954577320334
lon_Ed = -113.49274885945681
#Salekhard
lat_Sal=66.5333
lon_Sal=66.6333

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
# define target resolution for finding intersections of auroral boundaries
target_res_zero = 1e-9
# parameters for bisection code
lon_res_b = np.deg2rad(0.5)
lons_b = np.linspace(0,2*np.pi,num=int(2*np.pi/lon_res_b))


# marker size for city locations
MS=10
# font size for city location
FSL = 20



# lats and lons for plotting the coordinates
lats_lin = np.linspace(-np.pi/2,np.pi/2,num=180)
lons_lin = np.linspace(0,2*np.pi,num=360)

lats, lons = np.meshgrid(lats_lin,lons_lin)
theta = -lats+np.pi/2
    
# case and increase selection
l=1
m=0
cs = 'g'    


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
if __name__ == '__main__':

    def lonlat_degrees_to_xyz(lon_degrees,lat_degrees):
        lon = lon_degrees*(math.pi/180)
        lat = lat_degrees*(math.pi/180)
        coslat = math.cos(lat)
        return [coslat*math.cos(lon), coslat*math.sin(lon), math.sin(lat)]

    tic = time.time()  
        
    plt.close('all')

    # this was a good test...
    folder1   = '../models/IPGPforecast/2015/'
    filename1 = 'ipgpMF4aacgmv2_from1590.txt'
    filename2 = 'ipgpMF4aacgmv2_from1590_axial_dipole.txt'
    
    
    # for the full magmodel as background
    folder3   = '../models/magmodel/'
    filename3 = 'magmodel_1590-2020.txt'
    '''
    # for the dipolar magmodel as background
    folder3   = '../models/magmodel/'
    #filename3 = 'magmodel_1590-2020_dipole.txt'
    #filename3 = 'magmodel_1590-2020_axial_dipole.txt'
    filename3 = 'magmodel_1590-2020_axial_dipole_smallH.txt'
    '''
    folder_base = './'
    folder_out = folder_base+'aacgmv2_coords/'
    folder_models = folder_base+'models/'
    folder_figures = folder_base+'figures/'


    
    # for magmodel
    t2020 = 2020
    dtime2020 = dt.datetime(t2020, 1, 1)
    
    
    
    #################################
    # reference model: magmodel, 2020   
    ################################        
    header1, header2, header3, header4, MFcoeffs, SVcoeffs = aacgm_functions.read_magmodel_files(folder3+filename3)

    #############################
    # reference model auroral zones
    
    # initialise lats arrays: n x 2 where n is range(lons), and the other dimension contains the polward and equatorward boundaries
    latsN = np.zeros((lons_b.shape[0],2))
    latsS = np.zeros((lons_b.shape[0],2))
    
    # northern zone
    # polar edge
    my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t2020, 1, 1, np.rad2deg(lons_b), latNpol_target, target_res, latNpol_i, latNeq_i, max_iterations,folder_out))
    my_process3.start()
    my_process3.join() 
    # equatorial edge
    my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t2020, 1, 1, np.rad2deg(lons_b), latNeq_target, target_res, latNpol_i, latNeq_i, max_iterations,folder_out))
    my_process3.start()
    my_process3.join() 
    
    # southern zone
    # polar edge
    my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t2020, 1, 1, np.rad2deg(lons_b), latSpol_target, target_res, latSpol_i, latSeq_i, max_iterations,folder_out))
    my_process3.start()
    my_process3.join() 
    # equatorial edge
    my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t2020, 1, 1, np.rad2deg(lons_b), latSeq_target, target_res, latSpol_i, latSeq_i, max_iterations,folder_out))
    my_process3.start()
    my_process3.join() 
                     
    # read the calculated values
    latsN[:,0] = np.loadtxt(folder_out+'bisection_'+str(latNpol_target)+'_lats_'+filename3)
    latsN[:,1] = np.loadtxt(folder_out+'bisection_'+str(latNeq_target)+'_lats_'+filename3)
    os.system('rm '+folder_out+'bisection_'+str(latNpol_target)+'_lats_'+filename3)
    os.system('rm '+folder_out+'bisection_'+str(latNeq_target)+'_lats_'+filename3)
    
    latsS[:,0] = np.loadtxt(folder_out+'bisection_'+str(latSpol_target)+'_lats_'+filename3)
    latsS[:,1] = np.loadtxt(folder_out+'bisection_'+str(latSeq_target)+'_lats_'+filename3)
    os.system('rm '+folder_out+'bisection_'+str(latSpol_target)+'_lats_'+filename3)
    os.system('rm '+folder_out+'bisection_'+str(latSeq_target)+'_lats_'+filename3)
            
    # polar cap centroids
    # north
    vertsN = [lonlat_degrees_to_xyz(lons_b[i]*180/np.pi,latsN[i,0]) for i in range(latsN[:,0].shape[0]) ]
    #verts=np.flipud(verts)
    
    momentN = sphericalPolygonMoment(vertsN)
    print("momentN = "+str(momentN))
    print("centroidN unit direction = "+str(normalized(momentN))) 
    
    centroidN = normalized(momentN)
    latcN = math.asin(centroidN[2])*180/np.pi
    loncN = math.atan2(centroidN[1],centroidN[0])*180/np.pi
    
    # south
    vertsS = [lonlat_degrees_to_xyz(lons_b[i]*180/np.pi,latsS[i,0]) for i in range(latsS[:,0].shape[0]) ]
    vertsS=np.flipud(vertsS)
    
    momentS = sphericalPolygonMoment(vertsS)
    print("momentS = "+str(momentS))
    print("centroidS unit direction = "+str(normalized(momentS))) 
    
    centroidS = normalized(momentS)
    latcS = math.asin(centroidS[2])*180/np.pi
    loncS = math.atan2(centroidS[1],centroidS[0])*180/np.pi
    
    
    # dip poles
    _, _, _, _, MFcoeffs_magmodel, SVcoeffs_magmodel = aacgm_functions.read_magmodel_files(folder3+filename3)
    diff1 = abs(MFcoeffs_magmodel[:,0]-t2020)
    idx1 = np.where(diff1 == diff1.min())
    MFmat2020 = SH_library.lin2matCoeffs(MFcoeffs_magmodel[idx1[0][0],1:])
    #Br2020, Bt2020, Bp2020=SH_library.calcB(MFmat2020,theta,lons,r_a,r_a)
    #F2020 = np.sqrt(Br2020**2 + Bt2020**2 + Bp2020**2)
    #Incl2020 = np.arctan(-Br2020/np.sqrt(Bt2020**2 + Bp2020**2))*180/np.pi    
    
    lat_NP2020,lon_NP2020, lat_SP2020, lon_SP2020 = SH_library.find_poles_brute(MFmat2020)
    
    
    # geomagnetic poles

    aacgmv2.IGRF_COEFFS = folder3+filename3
    aacgmv2.wrapper.set_coeff_path(aacgmv2.IGRF_COEFFS,True)          
    m2 = aacgmv2.deprecated.igrf_dipole_axis(dtime2020)
    lat2 = np.pi/2 - np.arccos(m2[2])
    lon2 = math.atan2(m2[1],m2[0])
    
    
    ## plot zones 
    fig = plt.figure(figsize=(15,8))
    
    # north pole
    ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.Orthographic(0, 90))
    
    ax1.set_global() # important (apparently)
    ax1.coastlines()
    ax1.gridlines(crs=ccrs.PlateCarree(), 
                      linewidth=1, color='black', linestyle=':')    
        
    ax1.fill_between(np.rad2deg(lons_b), latsN[:,0], latsN[:,1], 
         alpha=0.5,color='red', linewidth =0,
         transform=ccrs.PlateCarree())
    
    ax1.scatter(loncN, latcN
                   ,s=70,edgecolor='k',color='red',alpha=1, zorder=3,marker='s'
                   ,transform=ccrs.PlateCarree())
        
    ax1.scatter(lon_NP2020, lat_NP2020,s=110
               ,edgecolor='k',color='blue',alpha=1, zorder=3,marker='v'
               ,transform=ccrs.PlateCarree())
    ax1.scatter(lon2*180./np.pi, lat2*180./np.pi
               ,s=110,edgecolor='k',color='green',alpha=1, zorder=3,marker='d'
               ,transform=ccrs.PlateCarree())
    
    # cities location
    
    ax1.plot(lon_Leeds,lat_Leeds,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Leeds+1.5,lat_Leeds-7,"Leeds",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    ax1.plot(lon_Ed,lat_Ed,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Ed+7,lat_Ed-37,"Edmonton",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    ax1.plot(lon_Sal,lat_Sal,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Sal+1,lat_Sal-3,"Salekhard",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    
    #ax1.text(-50, 40, '2020+', color='blue',fontsize=18,transform=ccrs.PlateCarree())
    ax1.text(95, 68, '2020', color='red',fontsize=18,transform=ccrs.PlateCarree())
    
    ax1.set_title('North pole',fontsize=20)
    
    
    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.Orthographic(0, -90))
    
    ax2.set_global() # important (apparently)
    ax2.coastlines()
    ax2.gridlines(crs=ccrs.PlateCarree(), 
                      linewidth=1, color='black', linestyle=':')    
    
    ax2.fill_between(np.rad2deg(lons_b), latsS[:,0], latsS[:,1], 
         alpha=0.5,color='red', linewidth =0,
         transform=ccrs.PlateCarree())                       
    ax2.scatter(loncS, latcS
                   ,s=70,edgecolor='k',color='red',alpha=1, zorder=3,marker='s'
                   ,transform=ccrs.PlateCarree(), label='centroid')
    ax2.scatter(lon_SP2020, lat_SP2020,s=110
               ,edgecolor='k',color='blue',alpha=1, zorder=3,marker='v'
               ,transform=ccrs.PlateCarree(), label='magnetic dip poles')
    ax2.scatter(180+lon2*180./np.pi, - lat2*180./np.pi
               ,s=110,edgecolor='k',color='green',alpha=1, zorder=3,marker='d'
               ,transform=ccrs.PlateCarree(), label='geomagnetic poles')

    # cities location
    ax2.plot(lon_Dun,lat_Dun,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax2.text(lon_Dun+10,lat_Dun-4,"Dunedin",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    
    #ax2.text(45,-85, '2020+', color='blue',fontsize=18,transform=ccrs.PlateCarree())
    ax2.text(-60,-68, '2020', color='red',fontsize=18,transform=ccrs.PlateCarree())
    
    ax2.set_title('South pole',fontsize=20)
    
    #ax2.legend(loc='center left', bbox_to_anchor=(-0.75, -0.05),ncol=3,fontsize=15)
    ax2.legend(loc='center left', bbox_to_anchor=(-0.3, 0.08),fontsize=15, frameon=False)
    plt.savefig(folder_figures+'polar_zones_poles_res_lon='+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)    
    
    plt.show(block=False)             
            
    
    
    
    
    # 2020 latitudes
    ref_fig_name_lat = folder_figures + 'AACGMlats_igrf13_2020.pdf'
    my_process = Process(target=aacgm_functions.print_aacgm_coords, args=(folder3,filename3,folder_out,t2020,lats,lons))
    my_process.start()
    my_process.join()   
    
    CGMlats2020 = np.loadtxt(folder_out+'CGMlats_'+str(t2020)+'_'+filename3)
    
    # plot latitude
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines(linewidth=0.5)
    ax.gridlines()   
    
    #Mollweide projectionfrom scipy.interpolate import griddata
    ax.set_facecolor('silver')
    cf = plt.contourf(np.rad2deg(lons),np.rad2deg(lats),CGMlats2020,cmap='bwr', 
                       levels=np.linspace(-90,90,19)
                       ,transform=ccrs.PlateCarree())
    
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsN[:,i], '--',
                 color='k',
                 transform=ccrs.PlateCarree())

        ax.plot(np.rad2deg(lons_b), latsS[:,i], '--',
                 color='k',
                 transform=ccrs.PlateCarree())

    clb = plt.colorbar(cf
                       , pad=0.04
                       , aspect = 75
                       , orientation='horizontal')
    clb.ax.tick_params(labelsize=15)
    clb.set_label('AACGM latitude',fontsize=20)


    #plt.title('2020, AACGM latitudes (IGRF13)',fontsize=20)
    plt.savefig(ref_fig_name_lat,bbox_inches='tight',pad_inches=0.1)
    plt.show(block=False)    
    
    
    
    ########################3
    # test angular distance 
    
    dists = [aacgm_functions.angular_distance(lon_Leeds,lat_Leeds,180*lons_b[i]/np.pi,latsN[i,1]) for i in range(lons_b.shape[0])]
    idx_min = np.argmin(dists)
    dist_min = dists[idx_min]