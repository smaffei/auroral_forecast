#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 16:40:01 2021

@author: stefanomaffei

Green's function approach to analyse the sensitivity of the polar caps centroids
and of the distance between some cities and the auroral zones
to different Gauss coefficient's time variations

Uses a bisection algorithm to find the edges of the auroral zones

Uses an algorithm from Bevis and Cambareri, 1987 (and used in Zossi, Fagre, Amit and Elias, 2020)
to compute areas of the zones

Based on the 2021_green_bisection_areas_full.py script

This script does not run the full kernel calculation, but it checks the convergence of the calculation

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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.special import roots_legendre, eval_legendre
from scipy.spatial.transform import Rotation as R

import SH_library


t1t = 2020


# parameters
r_cmb = 3485.0e3
r_c=r_cmb
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

names_loc = ["Salekhard","Leeds","Dunedin","Edmonton"]
lats_loc = np.array([lat_Sal,lat_Leeds,lat_Dun,lat_Ed])
lons_loc = np.array([lon_Sal,lon_Leeds,lon_Dun,lon_Ed])

lats_loc = np.reshape(lats_loc, (1, lats_loc.shape[0]))
lons_loc = np.reshape(lons_loc, (1, lons_loc.shape[0]))

# marker size for city locations
MS=10
# font size for city location
FSL = 20


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
# relative increase of each coefficients (see the single case experiments for this)
# for exponential dg inecrease: more stable
#rel_increase = np.array([-10,-1,-0.1,0.1,1,10])
#rel_increase = np.array([1])


##############################
# case and increase selection
##############################
l=2
m=0 # not sure it's necessary
cs = 'g'   
incr_rel = 150 # l=3
#incr_rel = 300 # l=3
#incr_rel = 100 # l=1,2
rel_increase = 1 # an additional handle: not necessary
Lmax = 13


# change plot fonts globally   
#plt.rcParams["font.family"] = "Times"
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
plt.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
plt.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'


# define things for centroid calculation
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

def lonlat_degrees_to_xyz(lon_degrees,lat_degrees):
    lon = lon_degrees*(math.pi/180)
    lat = lat_degrees*(math.pi/180)
    coslat = math.cos(lat)
    return [coslat*math.cos(lon), coslat*math.sin(lon), math.sin(lat)]

# main function    
if __name__ == "__main__":   
    
    plt.close('all')


    # for the dipolar magmodel as background
    folder3   = '../models/magmodel/'
    filename3 = 'magmodel_1590-2020.txt'
    #filename3 = 'magmodel_1590-2020_dipole.txt'
    #filename3 = 'magmodel_1590-2020_axial_dipole_smallH.txt'
    #filename3 = 'magmodel_1590-2020_axial_dipole.txt'

    folder_base = './'
    folder_out = folder_base+'aacgmv2_coords/'
    folder_models = folder_base+'models/'
    
    folder_green = folder_base+'Green_time_magmodel2020_bisection/'
    #folder_green = folder_base+'Green_time_magmodel2020axial_dipole_bisection/'
    folder_quantity = folder_green + 'auroral_zones/'
    
    if not os.path.exists(folder_quantity):
        os.makedirs(folder_quantity)
        
   
    # for magmodel
    t2015 = 2015
    t2020 = 2020
    t1970 = 1970
    dtime2015 = dt.datetime(t2015, 1, 1)        
    dtime2020 = dt.datetime(t2020, 1, 1)        
    
    # prepare bisecting longitude. If the reference zones are already calculated, load them        
    filename_latsN = folder_quantity+"latsN_"+str(t2020)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3
    filename_latsS = folder_quantity+"latsS_"+str(t2020)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3

    if os.path.exists(filename_latsN) and os.path.exists(filename_latsS):
        lons_b = np.loadtxt(filename_latsN)[:,0]*np.pi/180
    else:
        lons_b = np.linspace(0,2*np.pi,num=int(2*np.pi/lon_res_b))
    
    '''
    # lats and lons for plotting the coordinates
    lats_lin = np.linspace(-np.pi/2,np.pi/2,num=180)
    lons_lin = np.linspace(0,2*np.pi,num=360)
    
    lats, lons = np.meshgrid(lats_lin,lons_lin)
    theta = -lats+np.pi/2
    '''
    
    # quadrature based grid
    nleg = 180
    costh, weights = roots_legendre(nleg)
    theta_lin = np.arccos(costh)
    lats_lin = -theta_lin + np.pi/2

    lons_lin = np.linspace(0,2*np.pi,num=360)

    theta, lons = np.meshgrid(theta_lin,lons_lin)
    lats = -theta + np.pi/2

    #################################
    # reference model: magmodel, 2020   
    ################################        
    header1, header2, header3, header4, MFcoeffs, SVcoeffs = aacgm_functions.read_magmodel_files(folder3+filename3)
    
    # find temporal index in array
    idx_time = np.where(MFcoeffs[:,0]==t2020)
    idx_time = idx_time[0][0]
    
    # current model SV
    SVcoeffs3 =  (MFcoeffs[idx_time,1:] - MFcoeffs[idx_time-1,1:])/(MFcoeffs[idx_time,0] - MFcoeffs[idx_time-1,0])
    

    ########################################################
    # Green's function approach for the temporal variations   
    #######################################################    

    ##############################
    # initialise stuff

    #dg =   np.array([])
    #g  =   np.array([])
    
    #############################
    # reference model auroral zones

    if os.path.exists(filename_latsN) and os.path.exists(filename_latsS):
        print('loading reference auroral zone (IGRF)')
        latsN = np.loadtxt(filename_latsN)[:,1:3]
        latsS = np.loadtxt(filename_latsS)[:,1:3]
    else:
        # calculate them
        print('calculating reference auroral zone (IGRF)')
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
    
        np.savetxt(filename_latsN,np.hstack([lons_b[:,None]*180/np.pi, latsN]),header='lons lats(polar edge) lats(equatorward edge)')
        #os.system('rm '+folder_out+'bisection_'+str(latNpol_target)+'_lats_'+filename3)
        #os.system('rm '+folder_out+'bisection_'+str(latNeq_target)+'_lats_'+filename3)
        
        
        latsS[:,0] = np.loadtxt(folder_out+'bisection_'+str(latSpol_target)+'_lats_'+filename3)
        latsS[:,1] = np.loadtxt(folder_out+'bisection_'+str(latSeq_target)+'_lats_'+filename3)
        
        np.savetxt(filename_latsS,np.hstack([lons_b[:,None]*180/np.pi, latsS]),header='lons lats(polar edge) lats(equatorward edge)')
        #os.system('rm '+folder_out+'bisection_'+str(latSpol_target)+'_lats_'+filename3)
        #os.system('rm '+folder_out+'bisection_'+str(latSeq_target)+'_lats_'+filename3)
        
    
        
        
    ## polar cap centroids
    # north
    vertsN = [lonlat_degrees_to_xyz(lons_b[i]*180/np.pi,latsN[i,0]) for i in range(latsN[:,0].shape[0]) ]    
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
    
    ## southern edge distance from selected locations
    idx_min=np.zeros(lons_loc.shape)
    dist_min=np.zeros(lons_loc.shape)
    for loc in range(lats_loc.shape[1]):
        if lats_loc[0,loc]<0:
            latsZ = latsS
        else:
            latsZ = latsN
        dists = [aacgm_functions.angular_distance(lons_loc[0,loc],lats_loc[0,loc],180*lons_b[i]/np.pi,latsZ[i,1]) for i in range(lons_b.shape[0])] 
        idx_min[0,loc] = np.argmin(dists)
        dist_min[0,loc] = dists[int(idx_min[0,loc])]
    
    print('')    
    
    # name the folder for this coefficent
    #folder_coeff = folder_quantity+cs+'_'+str(l)+'_'+str(m)+'/'
    folder_coeff = folder_quantity+cs+'_'+str(l)+'/'
    
    if not os.path.exists(folder_coeff):
        os.makedirs(folder_coeff)
    
    tic = time.time()
    
    
    if cs=='g':
        cs_idx = 0
    else:
        cs_idx = 1    
    
    # calculate the index position for the given l,m,cs
    # not necessary in a for loop
    coeffs_in_l = lambda x : 2*x+1
    if m==0:
        m_pos = 1
    else:
        m_pos = 2*m    
    idx = int( sum(coeffs_in_l(np.linspace(1,l-1,l-1))) + m_pos + cs_idx )
    
    # coeffs to change
    if l==1:
        idx1=0
    else:
        idx1 = int( sum(coeffs_in_l(np.linspace(1,l-1,l-1))))
    
    idx2 = idx1+coeffs_in_l(l)
    
    #idx1=3
    #idx2=4
        
    # define increment
    increase = incr_rel * rel_increase

    print('case does not exist, calculating:')
    # if the calculation goes on, update this index
                        
    ####################################################
    # create mymodel where I only modify one coefficient
    # modify one coefficient
    print('perturbating coefficient')

    MFcoeffs_plus = 1*MFcoeffs
    # to change all l coeffs
    MFcoeffs_plus[-1,idx1+1:idx2+1] = MFcoeffs_plus[-1,idx1+1:idx2+1]+increase*SVcoeffs3[idx1:idx2]
    # to change a single coeff
    #MFcoeffs_plus[-1,idx] = MFcoeffs_plus[-1,idx]+dg[inotE]
    # save mymodel
    aacgm_functions.write_magmodel_files(header1, header2, header3, header4, MFcoeffs_plus, SVcoeffs, folder_quantity+'mymodel.txt')
    
    
    #############################
    # mymodel auroral zones
    filename_latsNplus = folder_coeff+"latsNplus_"+str(t2020)+"_dbeta_l_"+str(l)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3
    filename_latsSplus = folder_coeff+"latsSplus_"+str(t2020)+"_dbeta_l_"+str(l)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3

    if os.path.exists(filename_latsNplus) and os.path.exists(filename_latsSplus):
        print('loading perturbed auroral zone')
        latsNplus = np.loadtxt(filename_latsNplus)[:,1:3]
        latsSplus = np.loadtxt(filename_latsSplus)[:,1:3]
    else:
        print('calculating perturbed auroral zone')
    
        # initialise lats arrays: n x 2 where n is range(lons), and the other dimension contains the polward and equatorward boundaries
        latsNplus = np.zeros((lons_b.shape[0],2))
        latsSplus = np.zeros((lons_b.shape[0],2))
        
        # northern zone
        # polar edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder_quantity,'mymodel.txt', t2020, 1, 1, np.rad2deg(lons_b), latNpol_target, target_res, latNpol_i, latNeq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        # equatorial edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder_quantity,'mymodel.txt', t2020, 1, 1, np.rad2deg(lons_b), latNeq_target, target_res, latNpol_i, latNeq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        
        # southern zone
        # polar edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder_quantity,'mymodel.txt', t2020, 1, 1, np.rad2deg(lons_b), latSpol_target, target_res, latSpol_i, latSeq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        # equatorial edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder_quantity,'mymodel.txt', t2020, 1, 1, np.rad2deg(lons_b), latSeq_target, target_res, latSpol_i, latSeq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
                         
        # read the calculated values
        latsNplus[:,0] = np.loadtxt(folder_out+'bisection_'+str(latNpol_target)+'_lats_'+'mymodel.txt')
        latsNplus[:,1] = np.loadtxt(folder_out+'bisection_'+str(latNeq_target)+'_lats_'+'mymodel.txt')
        
        #np.savetxt(filename_latsNplus,np.hstack([lons_b[:,None]*180/np.pi, latsNplus]),header='lons lats(polar edge) lats(equatorward edge)')
        os.system('rm '+folder_out+'bisection_'+str(latNpol_target)+'_lats_'+'mymodel.txt')
        os.system('rm '+folder_out+'bisection_'+str(latNeq_target)+'_lats_'+'mymodel.txt')
        
        latsSplus[:,0] = np.loadtxt(folder_out+'bisection_'+str(latSpol_target)+'_lats_'+'mymodel.txt')
        latsSplus[:,1] = np.loadtxt(folder_out+'bisection_'+str(latSeq_target)+'_lats_'+'mymodel.txt')

        #np.savetxt(filename_latsSplus,np.hstack([lons_b[:,None]*180/np.pi, latsSplus]),header='lons lats(polar edge) lats(equatorward edge)')
        os.system('rm '+folder_out+'bisection_'+str(latSpol_target)+'_lats_'+'mymodel.txt')
        os.system('rm '+folder_out+'bisection_'+str(latSeq_target)+'_lats_'+'mymodel.txt')
        #os.system('rm '+folder_quantity+'mymodel.txt')
    

    ## polar cap centroids
    # north
    vertsNplus = [lonlat_degrees_to_xyz(lons_b[i]*180/np.pi,latsNplus[i,0]) for i in range(latsNplus[:,0].shape[0]) ]    
    momentNplus = sphericalPolygonMoment(vertsNplus)
    centroidNplus = normalized(momentNplus)
    latcNplus = math.asin(centroidNplus[2])*180/np.pi
    loncNplus = math.atan2(centroidNplus[1],centroidNplus[0])*180/np.pi
    
    # south
    vertsSplus = [lonlat_degrees_to_xyz(lons_b[i]*180/np.pi,latsSplus[i,0]) for i in range(latsSplus[:,0].shape[0]) ]
    vertsSplus=np.flipud(vertsSplus)    
    momentSplus = sphericalPolygonMoment(vertsSplus)    
    centroidSplus = normalized(momentSplus)
    latcSplus = math.asin(centroidSplus[2])*180/np.pi
    loncSplus = math.atan2(centroidSplus[1],centroidSplus[0])*180/np.pi    
    
    ## southern edge distance from selected locations
    idx_min_plus=np.zeros(lons_loc.shape)
    dist_min_plus=np.zeros(lons_loc.shape)
    for loc in range(lats_loc.shape[1]):
        if lats_loc[0,loc]<0:
            latsZ = latsSplus
        else:
            latsZ = latsNplus
        dists = [aacgm_functions.angular_distance(lons_loc[0,loc],lats_loc[0,loc],180*lons_b[i]/np.pi,latsZ[i,1]) for i in range(lons_b.shape[0])] 
        idx_min_plus[0,loc] = np.argmin(dists)
        dist_min_plus[0,loc] = dists[int(idx_min[0,loc])]
    
    
    #########################################################################
    ## find intersections of zones' edges: 
    # assumes small variations (the polar edges never intersects the equatorial ones)  
    
    # calculate reference area 
    areaNpol_b = aacgm_functions.spherical_polygon_area(np.flip(latsN[:,0]),np.flip(np.rad2deg(lons_b)),r_a/1000)
    areaNeq_b  = aacgm_functions.spherical_polygon_area(np.flip(latsN[:,1]),np.flip(np.rad2deg(lons_b)),r_a/1000)
    areaN = areaNeq_b - areaNpol_b
    
    areaSpol_b = aacgm_functions.spherical_polygon_area(latsS[:,0],np.rad2deg(lons_b),r_a/1000)
    areaSeq_b  = aacgm_functions.spherical_polygon_area(latsS[:,1],np.rad2deg(lons_b),r_a/1000)
    areaS  = areaSeq_b - areaSpol_b
    
    # calculate perturbed areas
    areaNpol_b = aacgm_functions.spherical_polygon_area(np.flip(latsNplus[:,0]),np.flip(np.rad2deg(lons_b)),r_a/1000)
    areaNeq_b  = aacgm_functions.spherical_polygon_area(np.flip(latsNplus[:,1]),np.flip(np.rad2deg(lons_b)),r_a/1000)
    areaNplus = areaNeq_b - areaNpol_b
    
    areaSpol_b = aacgm_functions.spherical_polygon_area(latsSplus[:,0],np.rad2deg(lons_b),r_a/1000)
    areaSeq_b  = aacgm_functions.spherical_polygon_area(latsSplus[:,1],np.rad2deg(lons_b),r_a/1000)
    areaSplus = areaSeq_b - areaSpol_b  
    

    ###########################################################################              
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
        
    ax1.fill_between(np.rad2deg(lons_b), latsNplus[:,0], latsNplus[:,1], 
         alpha=0.5,color='blue', linewidth =0,
         transform=ccrs.PlateCarree())
    
    ax1.scatter(loncN, latcN
                   ,s=70,edgecolor='k',color='red',alpha=1, zorder=3,marker='s'
                   ,transform=ccrs.PlateCarree())
        
    ax1.scatter(loncNplus, latcNplus
               ,s=70,edgecolor='k',color='blue',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    ax1.plot(lon_Leeds,lat_Leeds,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Leeds+1.5,lat_Leeds-7,"Leeds",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    ax1.plot(lon_Ed,lat_Ed,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Ed+7,lat_Ed-37,"Edmonton",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    ax1.plot(lon_Sal,lat_Sal,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Sal+1,lat_Sal-3,"Salekhard",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())


    
    #ax1.text(-50, 40, '2020+', color='blue',fontsize=18,transform=ccrs.PlateCarree())
    #ax1.text(95, 68, '2020', color='red',fontsize=18,transform=ccrs.PlateCarree())
    
    ax1.set_title('North pole',fontsize=20)
    
    
    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.Orthographic(0, -90))
    
    ax2.set_global() # important (apparently)
    ax2.coastlines()
    ax2.gridlines(crs=ccrs.PlateCarree(), 
                      linewidth=1, color='black', linestyle=':')    

    ax2.fill_between(np.rad2deg(lons_b), latsS[:,0], latsS[:,1], 
         alpha=0.5,color='red', linewidth =0,
         transform=ccrs.PlateCarree())       
    
    ax2.fill_between(np.rad2deg(lons_b), latsSplus[:,0], latsSplus[:,1], 
         alpha=0.5,color='blue', linewidth =0,
         transform=ccrs.PlateCarree())       
        
    ax2.scatter(loncS, latcS
                   ,s=70,edgecolor='k',color='red',alpha=1, zorder=3,marker='s'
                   ,transform=ccrs.PlateCarree())
    ax2.scatter(loncSplus, latcSplus
               ,s=70,edgecolor='k',color='blue',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    ax2.plot(lon_Dun,lat_Dun,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax2.text(lon_Dun+10,lat_Dun-4,"Dunedin",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())


    
    #ax2.text(45,-85, '2020+', color='blue',fontsize=18,transform=ccrs.PlateCarree())
    #ax2.text(-60,-70, '2020', color='red',fontsize=18,transform=ccrs.PlateCarree())
    
    ax2.set_title('South pole',fontsize=20)
    
    plt.savefig(folder_quantity+'polar_zones_dbetal_'+str(l)+'_incr_multiplier_'+str(increase)+'_res_lon='+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)    
    
    plt.show(block=False)        
    
    
    ###########################
    # what SV am I considering?
                
    SV_coeffs_mat = SH_library.lin2matCoeffs(MFcoeffs_plus[-1,1:] - MFcoeffs[-1,1:] )
    
    dBr_cdt, _, _ = SH_library.calcB(SV_coeffs_mat,theta,lons,r_a,r_cmb)
            
    fig = plt.figure(figsize=(10,6))
    Zs = dBr_cdt/1000
    Zmax = np.nanmax(abs(Zs))
    Zmin = -np.nanmax(abs(Zs))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    

    cs = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs,cmap='coolwarm' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    clb = plt.colorbar(cs
                       #,ticks=[20,30,40,50,60,70]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       )
    #clb.set_label(r'$\dot{B_c}$  [$\mu$T/year]', fontsize = 20)
    clb.set_label(str(incr_rel)+r'$ \ \times \  \dot{B_c}$  [$\mu$T/year]', fontsize = 20)
    #clb.ax.set_xticklabels(['20','30','40','50','60','70'])
    clb.ax.tick_params(labelsize=20)
    plt.title(r'$l=$'+str(l),fontsize=20)
    plt.savefig(folder_quantity+'CMB_SV_dbetal_'+str(l)+'_incr_multiplier_'+str(increase),bbox_inches='tight',pad_inches=0.1)    

    plt.show()
        
        
        
        
    #############################################
    toc = time.time()
    elapsed_time = toc-tic
    print('elapsed time for current Gauss coefficient: '+str(elapsed_time)+' seconds')
    # end of rel_increase
    

                       
  
   
    ##########
    # cleanup    
    #os.system('rm '+folder_out+'CGMlats_'+str(t2020)+'_'+filename3)
    #os.system('rm '+folder_out+'CGMlons_'+str(t2020)+'_'+filename3)
    #os.system('rm '+folder_out+'CGMmlt_'+str(t2020)+'_'+filename3)
    
    #os.system('rm '+folder_quantity+'CGMlats_'+str(t2020)+'_'+'mymodel.txt')
    #os.system('rm '+folder_quantity+'CGMlons_'+str(t2020)+'_'+'mymodel.txt')
    #os.system('rm '+folder_quantity+'CGMmlt_'+str(t2020)+'_'+'mymodel.txt')
    
    ########################
    # to check:
    ########################
    
    # plot the CGM latitudes
    '''
    
    my_process = Process(target=aacgm_functions.print_aacgm_coords, args=(folder3,filename3,folder_out,t2020,lats,lons))
    my_process.start()
    my_process.join()   
    
    CGMlats2020 = np.loadtxt(folder_out+'CGMlats_'+str(t2020)+'_'+filename3)    
    
    my_process = Process(target=aacgm_functions.print_aacgm_coords, args=(folder_quantity,'mymodel.txt',folder_out,t2020,lats,lons))
    my_process.start()
    my_process.join()   
    
    CGMlats_plus = np.loadtxt(folder_out+'CGMlats_'+str(t2020)+'_'+'mymodel.txt')    


    fig = plt.figure(figsize=(10,6))
    Zs = CGMlats_plus
    Zmax = np.nanmax(abs(Zs))
    Zmin = -np.nanmax(abs(Zs))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    

    cs = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs,cmap='coolwarm' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    clb = plt.colorbar(cs
                       #,ticks=[20,30,40,50,60,70]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       )
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsN[:,i], '--',
                 color='k',
                 transform=ccrs.PlateCarree())
        
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsS[:,i], '--',
                 color='k',
                 transform=ccrs.PlateCarree())
        
    #clb.set_label(r'$\dot{B_c}$  [$\mu$T/year]', fontsize = 20)
    clb.set_label(r'CGM latitudes ', fontsize = 20)
    #clb.ax.set_xticklabels(['20','30','40','50','60','70'])
    clb.ax.tick_params(labelsize=20)
    plt.title(r'$l=$'+str(l),fontsize=20)
    #plt.savefig(folder_quantity+'CMB_SV_dbetal_'+str(l)+'_incr_multiplier_'+str(increase),bbox_inches='tight',pad_inches=0.1)    

    plt.show()
    
    
    

    
    fig = plt.figure(figsize=(10,6))
    #Zs = CGMlats_plus
    Zs = CGMlats2020
    Zmax = np.nanmax(abs(Zs))
    Zmin = -np.nanmax(abs(Zs))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    

    cs = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs,cmap='coolwarm' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    clb = plt.colorbar(cs
                       #,ticks=[20,30,40,50,60,70]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       )
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsNplus[:,i], '--',
                 color='k',
                 transform=ccrs.PlateCarree())
        
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsSplus[:,i], '--',
                 color='k',
                 transform=ccrs.PlateCarree())
        
    #clb.set_label(r'$\dot{B_c}$  [$\mu$T/year]', fontsize = 20)
    clb.set_label(r'CGM latitudes (plus)', fontsize = 20)
    #clb.ax.set_xticklabels(['20','30','40','50','60','70'])
    clb.ax.tick_params(labelsize=20)
    plt.title(r'$l=$'+str(l),fontsize=20)
    #plt.savefig(folder_quantity+'CMB_SV_dbetal_'+str(l)+'_incr_multiplier_'+str(increase),bbox_inches='tight',pad_inches=0.1)    

    plt.show()
    
    
    ###########################################################################              
    ## plot zones 
    fig = plt.figure(figsize=(15,8))
    
    # north pole
    ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.Orthographic(0, 90))
    
    ax1.set_global() # important (apparently)
    ax1.coastlines()
    ax1.gridlines(crs=ccrs.PlateCarree(), 
                      linewidth=1, color='black', linestyle=':')    

    ax1.contourf(lons*180/np.pi,lats*180/np.pi,CGMlats2020, colors='red',alpha=0.4
                   , vmin = 65, vmax = 70
                   ,levels=np.linspace(65,70,41)
                   ,transform=ccrs.PlateCarree()
                   )
    ax1.contourf(lons*180/np.pi,lats*180/np.pi,CGMlats_plus
                 , colors='blue'
                 #,cmap = cMap
                 , alpha=0.2
                   , vmin = 65, vmax = 70
                   ,levels=np.linspace(65,70,41)
                   ,transform=ccrs.PlateCarree()
                   )
        
    #ax1.fill_between(np.rad2deg(lons_b), latsNplus[:,0], latsNplus[:,1], 
    #     alpha=0.5,color='blue', linewidth =0,
    #     transform=ccrs.PlateCarree())
    
    ax1.scatter(loncN, latcN
                   ,s=70,edgecolor='k',color='red',alpha=1, zorder=3,marker='s'
                   ,transform=ccrs.PlateCarree())
        
    ax1.scatter(loncNplus, latcNplus
               ,s=70,edgecolor='k',color='blue',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    ax1.plot(lon_Leeds,lat_Leeds,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Leeds+1.5,lat_Leeds-7,"Leeds",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    ax1.plot(lon_Ed,lat_Ed,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Ed+7,lat_Ed-37,"Edmonton",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    ax1.plot(lon_Sal,lat_Sal,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Sal+1,lat_Sal-3,"Salekhard",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())


    
    #ax1.text(-50, 40, '2020+', color='blue',fontsize=18,transform=ccrs.PlateCarree())
    #ax1.text(95, 68, '2020', color='red',fontsize=18,transform=ccrs.PlateCarree())
    
    ax1.set_title('North pole',fontsize=20)
    
    
    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.Orthographic(0, -90))
    
    ax2.set_global() # important (apparently)
    ax2.coastlines()
    ax2.gridlines(crs=ccrs.PlateCarree(), 
                      linewidth=1, color='black', linestyle=':')    

    ax2.contourf(lons*180/np.pi,lats*180/np.pi,CGMlats2020, colors='red',alpha=0.5
                   , vmin = -70, vmax = -65
                   ,levels=np.linspace(-70,-65,41)
                   ,transform=ccrs.PlateCarree()
                   )
    
    ax2.contourf(lons*180/np.pi,lats*180/np.pi,CGMlats_plus, colors='blue',alpha=0.5
                   , vmin = -70, vmax = -65
                   ,levels=np.linspace(-70,-65,41)
                   ,transform=ccrs.PlateCarree()
                   )   
    
    
        
    ax2.scatter(loncS, latcS
                   ,s=70,edgecolor='k',color='red',alpha=1, zorder=3,marker='s'
                   ,transform=ccrs.PlateCarree())
    ax2.scatter(loncSplus, latcSplus
               ,s=70,edgecolor='k',color='blue',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    ax2.plot(lon_Dun,lat_Dun,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax2.text(lon_Dun+10,lat_Dun-4,"Dunedin",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())


    
    #ax2.text(45,-85, '2020+', color='blue',fontsize=18,transform=ccrs.PlateCarree())
    #ax2.text(-60,-70, '2020', color='red',fontsize=18,transform=ccrs.PlateCarree())
    
    ax2.set_title('South pole',fontsize=20)
    
    #plt.savefig(folder_quantity+'polar_zones_dbetal_'+str(l)+'_incr_multiplier_'+str(increase)+'_res_lon='+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)    
    
    plt.show(block=False)        
    
    '''
    
    
    ########################
    # alternative implementations:
    ########################
    
    
    
    ############################################
    ### need to update the bisection function!!!
    
    
    # Initial, low-res guesses (I think it's unavoidable)
    lats_lin_lowres = np.linspace(-np.pi/2,np.pi/2,num=50)
    lons_lin_lowres = np.linspace(0,2*np.pi,num=45)
    
    lats_lowres, lons_lowres = np.meshgrid(lats_lin_lowres,lons_lin_lowres)
    theta_lowres = -lats_lowres+np.pi/2
    
    
    #my_process = Process(target=aacgm_functions.print_aacgm_coords, args=(folder_quantity,'mymodel.txt',folder_out,t2020,lats_lowres,lons_lowres))
    my_process = Process(target=aacgm_functions.print_aacgm_coords, args=(folder_quantity,'mymodel.txt',folder_out,t2020,lats_lowres,lons_lowres))
    my_process.start()
    my_process.join()   
    
    CGMlats_lowres = np.loadtxt(folder_out+'CGMlats_'+str(t2020)+'_'+'mymodel.txt')
    
    # make it general
    CGMlats = CGMlats_lowres
    lons1 =lons_lowres
    lats1 =lats_lowres
    
    # use previous calculations
    #CGMlats = CGMlats2020
    #lons1 =lons
    #lats1 =lats
    
    # find the zones corresponding to this low res case:
    zoneSpol = [] #  empty regular list
    zoneSeq  = [] #  empty regular list
    
    for i in range(lons1.shape[0]):
        CGMlats_prev = CGMlats[i,0]
        for j in range(lats1.shape[1]-1):
            # find crossings
            if (CGMlats_prev - latSpol_target)*(CGMlats[i,j+1] - latSpol_target)<=0:
                zoneSpol.append([lons1[i,j], lats1[i,j]])
                
            if (CGMlats_prev - latSeq_target)*(CGMlats[i,j+1] - latSeq_target)<=0:
                zoneSeq.append([lons1[i,j], lats1[i,j]])
        
            CGMlats_prev = CGMlats[i,j+1]
   
    
    zoneSpol = np.array(zoneSpol)            
    zoneSeq  = np.array(zoneSeq)            

    # we need to order these by nearest neighbours
    
    new_arr_pol = np.zeros((zoneSpol.shape[0]+1, zoneSpol.shape[1]))
    new_arr_pol[0] = 1*zoneSpol[0]
    # last element should overlap ( I think)
    new_arr_pol[-1] = 1*zoneSpol[0]
    # remove used element
    zoneSpol = np.delete(zoneSpol, 0, 0)
    
    
    for j in range(new_arr_pol.shape[0]-2):
        dists=np.zeros((zoneSpol.shape[0],))
        for i in range(dists.shape[0]):
            dists[i]=aacgm_functions.angular_distance(new_arr_pol[j,0],new_arr_pol[j,1],zoneSpol[i,0],zoneSpol[i,1])
        idx_closest = np.argmin(dists)
        new_arr_pol[j+1] = 1*zoneSpol[idx_closest]
        zoneSpol = np.delete(zoneSpol, idx_closest, 0)
            
    zoneSpol = 1* new_arr_pol
    
    
    # just copy-paste and modify for now
    new_arr_eq = np.zeros((zoneSeq.shape[0]+1, zoneSeq.shape[1]))
    new_arr_eq[0] = 1*zoneSeq[0]
    # last element should overlap ( I think)
    new_arr_eq[-1] = 1*zoneSeq[0]
    # remove used element
    zoneSeq = np.delete(zoneSeq, 0, 0)

    for j in range(new_arr_eq.shape[0]-2):
        dists=np.zeros((zoneSeq.shape[0],))
        for i in range(dists.shape[0]):
            dists[i]=aacgm_functions.angular_distance(new_arr_eq[j,0],new_arr_eq[j,1],zoneSeq[i,0],zoneSeq[i,1])
        idx_closest = np.argmin(dists)
        new_arr_eq[j+1] = 1*zoneSeq[idx_closest]
        zoneSeq = np.delete(zoneSeq, idx_closest, 0)
            
    zoneSeq = 1* new_arr_eq        
    

    # and this serves as first guess.
    # Now I need the second guess, to sandwhich the first guess.
    
    # given a point on the zone boundary, x_0, I want two points, x_1 and x_2, with small angular distance +-delta from x_0
    # x_0, x_1, x_2 and the centroid (x_c) all lie on the same great circle.
    
    # the parametric equation of a point on the great circle is (https://math.stackexchange.com/questions/1783746/equation-of-a-great-circle-passing-through-two-points)
    #
    # x = cos(t) x_c + sin(t) v_1
    #
    # and v_1 is found by orthogonalisation
    #
    # v = alpha x_c + beta x_0
    # with constaints
    # ||v||^2 = r_a^2
    # v.x_c = 0
    #
    # define d = x_c.x_0/r_a^2
    #
    # alpha = -beta * d
    # beta^2 = 1/(1-d^2) -> two solutions: v_1 and v_2
    
    # centroid:
    
    
    # south
    vertsSplus = [lonlat_degrees_to_xyz(zoneSpol[i,0]*180/np.pi,zoneSpol[i,1]*180/np.pi) for i in range(zoneSpol.shape[0]) ]
    #vertsSplus=np.flipud(vertsSplus) # apparently, comment it out......
    momentSplus = sphericalPolygonMoment(vertsSplus)    
    centroidSplus = normalized(momentSplus)
    latcSplus = math.asin(centroidSplus[2])*180/np.pi
    loncSplus = math.atan2(centroidSplus[1],centroidSplus[0])*180/np.pi    
    
    
    xc,yc,zc = aacgm_functions.spher2cart(r_a,latcSplus,loncSplus)
    
    
    # now I need to define the rotation matrices
    # rotation around z of angle lonc
    rz = R.from_matrix([[math.cos(loncSplus*np.pi/180),     math.sin(loncSplus*np.pi/180),  0],
                        [-math.sin(loncSplus*np.pi/180),    math.cos(loncSplus*np.pi/180),  0],
                        [0,                                 0,                              1]])
    # rotation around y of angle thetac
    ry = R.from_matrix([[math.cos((90-latcSplus)*np.pi/180),     0,                             -math.sin((90-latcSplus)*np.pi/180)],
                        [0,                                      1,                             0],
                        [math.sin((90-latcSplus)*np.pi/180),     0,                             math.cos((90-latcSplus)*np.pi/180)]])    
    
    # combined rotation:
    rtotal = ry * rz
    rinv   = rtotal.inv() 
    # as desired : rtotal.apply([xc,yc,zc]) is the positive z vector!
    
    
    # find first guesses, from the low-res auroral zone
         
    # polar border

    LATMIN_pol = 0
    LATMAX_pol = 0
    
    for i in range(zoneSpol.shape[0]):
        lon0 = zoneSpol[i,0]*180/np.pi
        lat0 = zoneSpol[i,1]*180/np.pi
        x0,y0,z0 = aacgm_functions.spher2cart(r_a,lat0,lon0)
        # transform this point to the alternative coordinate system
        X0,Y0,Z0=rtotal.apply([x0,y0,z0])
        
        _,LAT0,LON0 = aacgm_functions.cart2spher(X0,Y0,Z0)
        
        # check if max/min
        if LATMIN_pol==0:
            LATMIN_pol = LAT0
        if LATMAX_pol==0:
            LATMAX_pol = LAT0
            
        if LAT0<LATMIN_pol:
            LATMIN_pol = LAT0
        if LAT0>LATMAX_pol:
            LATMAX_pol = LAT0

    
    # perturb min/max
    LATMAX_pol = min([90, LATMAX_pol+5]) # just in case is verey close to the poles, but shold not happen in the transformed coordinates
    LATMIN_pol = LATMIN_pol-5
    
    
    # equatorial border
    LATMIN_eq = 0
    LATMAX_eq = 0
    
    for i in range(zoneSeq.shape[0]):
        lon0 = zoneSeq[i,0]*180/np.pi
        lat0 = zoneSeq[i,1]*180/np.pi
        x0,y0,z0 = aacgm_functions.spher2cart(r_a,lat0,lon0)
        # transform this point to the alternative coordinate system
        X0,Y0,Z0=rtotal.apply([x0,y0,z0])
        
        _,LAT0,LON0 = aacgm_functions.cart2spher(X0,Y0,Z0)
        
        # check if max/min
        if LATMIN_eq==0:
            LATMIN_eq = LAT0
        if LATMAX_eq==0:
            LATMAX_eq = LAT0
            
        if LAT0<LATMIN_eq:
            LATMIN_eq = LAT0
        if LAT0>LATMAX_eq:
            LATMAX_eq = LAT0
    
    # perturb min/max
    LATMAX_eq = min([90, LATMAX_eq+5]) # just in case is verey close to the poles, but shold not happen in the transformed coordinates
    LATMIN_eq = LATMIN_eq-5
    
    
    LONSB = np.linspace(0,2*np.pi,num=int(2*np.pi/lon_res_b))
    
    ### now run the bisection! ######
    
    
    
    # polar edge
    my_process3 = Process(target=aacgm_functions.lat_bisection_transformed_coords, args=(folder_quantity,'mymodel.txt', t2020, 1, 1, r_a, latcSplus, loncSplus, rtotal, np.rad2deg(lons_b), latSpol_target, target_res, LATMAX_pol, LATMIN_pol, max_iterations,folder_out))
    my_process3.start()
    my_process3.join() 
                     
    # read the calculated values
    coordsSplus_pol = np.loadtxt(folder_out+'bisection_transformed_coords_'+str(latSpol_target)+'_lats_lons_'+'mymodel.txt')
    # adjust direction
    coordsSplus_pol = np.flipud(coordsSplus_pol)

    # equatorial edge
    my_process3 = Process(target=aacgm_functions.lat_bisection_transformed_coords, args=(folder_quantity,'mymodel.txt', t2020, 1, 1, r_a, latcSplus, loncSplus, rtotal, np.rad2deg(lons_b), latSeq_target, target_res, LATMAX_eq, LATMIN_eq, max_iterations,folder_out))
    my_process3.start()
    my_process3.join() 
                     
    # read the calculated values
    coordsSplus_eq = np.loadtxt(folder_out+'bisection_transformed_coords_'+str(latSeq_target)+'_lats_lons_'+'mymodel.txt')
    # adjust direction
    coordsSplus_eq = np.flipud(coordsSplus_eq)
   
        
    
    #ax2.plot(np.rad2deg(lons_b), latsSplus[:,1],    

    # recalculate the area
    
    areaSpol_b = aacgm_functions.spherical_polygon_area(coordsSplus_pol[:,1],coordsSplus_pol[:,0],r_a/1000)
    areaSeq_b  = aacgm_functions.spherical_polygon_area(latsSplus[:,1],np.rad2deg(lons_b),r_a/1000)
    areaSplus_hr = areaSeq_b - areaSpol_b  


    
    # recalculate centroid:
    
    
    # south
    vertsSplus_hr = [lonlat_degrees_to_xyz(coordsSplus_pol[i,0],coordsSplus_pol[i,1]) for i in range(coordsSplus_pol.shape[0]) ]
    #vertsSplus=np.flipud(vertsSplus) # apparently, comment it out......
    momentSplus_hr = sphericalPolygonMoment(vertsSplus_hr)    
    centroidSplus_hr = normalized(momentSplus_hr)
    latcSplus_hr = math.asin(centroidSplus_hr[2])*180/np.pi
    loncSplus_hr = math.atan2(centroidSplus_hr[1],centroidSplus_hr[0])*180/np.pi    
    
    
    
    # brute force plot the zones
    
    Dlon  = np.deg2rad(0.5) # longitudinal res, multiply by 0.5 for the higher res
    DlatN = np.deg2rad(0.15) # lat res, northern cap, multiply by 0.5 for the higher res
    DlatS = DlatN # lat res, southern cap
    
    # longitude will be the same for both caps
    lonsNS_lin = np.arange(0,2*np.pi,Dlon)

    # define latitude for the northern Cap
    latNeq  = np.deg2rad(50)
    latNpol = np.deg2rad(85)        
    latsN_lin = np.arange(latNeq,latNpol,DlatN)
    
    # define latitude for the southern Cap
    latSeq  = np.deg2rad(-40)
    latSpol = np.deg2rad(-90)
    latsS_lin = np.arange(latSpol,latSeq,DlatS)
    
    #combine
    latsNS_lin = np.append(latsS_lin,latsN_lin)

    latsNS, lonsNS = np.meshgrid(latsNS_lin,lonsNS_lin)
    #latsNS, lonsNS = np.meshgrid(lats_lin,lons_lin)
    thetaNS = -latsNS+np.pi/2
    
    
    
    my_process = Process(target=aacgm_functions.print_aacgm_coords, args=(folder3,filename3,folder_out,t2020,latsNS,lonsNS))
    my_process.start()
    my_process.join()   
    
    CGMlats = np.loadtxt(folder_out+'CGMlats_'+str(t2020)+'_'+filename3)    
    
    my_process = Process(target=aacgm_functions.print_aacgm_coords, args=(folder_quantity,'mymodel.txt',folder_out,t2020,latsNS,lonsNS))
    my_process.start()
    my_process.join()   
    
    CGMlats_plus = np.loadtxt(folder_out+'CGMlats_'+str(t2020)+'_'+'mymodel.txt')   
    

    mask2020=np.array(CGMlats.copy())
    mask2020[np.isnan(CGMlats)] = 0 # get rid of nans
    mask2020[abs(CGMlats)>70] = 0 # mask to zero what is above 70
    mask2020[abs(CGMlats)<65] = 0 # mask to zero what is below 65
    mask2020[abs(mask2020)>=65] = 1 # mask to ones the rest
    
    mask2020plus=np.array(CGMlats_plus.copy())
    mask2020plus[np.isnan(CGMlats_plus)] = 0 # get rid of nans
    mask2020plus[abs(CGMlats_plus)>70] = 0 # mask to zero what is above 70
    mask2020plus[abs(CGMlats_plus)<65] = 0 # mask to zero what is below 65
    mask2020plus[abs(mask2020plus)>=65] = 1 # mask to ones the rest


    
    
    
    
    ###########################################################################              
    ## plot zones: only the difficult one is plotted brute force
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
        
    ax1.fill_between(np.rad2deg(lons_b), latsNplus[:,0], latsNplus[:,1], 
         alpha=0.5,color='blue', linewidth =0,
         transform=ccrs.PlateCarree())
    
    ax1.scatter(loncN, latcN
                   ,s=70,edgecolor='k',color='red',alpha=1, zorder=3,marker='s'
                   ,transform=ccrs.PlateCarree())
        
    ax1.scatter(loncNplus, latcNplus
               ,s=70,edgecolor='k',color='blue',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    ax1.plot(lon_Leeds,lat_Leeds,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Leeds+1.5,lat_Leeds-7,"Leeds",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    ax1.plot(lon_Ed,lat_Ed,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Ed+7,lat_Ed-37,"Edmonton",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    ax1.plot(lon_Sal,lat_Sal,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Sal+1,lat_Sal-3,"Salekhard",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())


    
    #ax1.text(-50, 40, '2020+', color='blue',fontsize=18,transform=ccrs.PlateCarree())
    #ax1.text(95, 68, '2020', color='red',fontsize=18,transform=ccrs.PlateCarree())
    
    ax1.set_title('North pole',fontsize=20)
    
    
    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.Orthographic(0, -90))
    
    ax2.set_global() # important (apparently)
    ax2.coastlines()
    ax2.gridlines(crs=ccrs.PlateCarree(), 
                      linewidth=1, color='black', linestyle=':')    

    ax2.fill_between(np.rad2deg(lons_b), latsS[:,0], latsS[:,1], 
         alpha=0.5,color='red', linewidth =0,
         transform=ccrs.PlateCarree())       
    
    cMap = c.ListedColormap(['blue'])
    ax2.pcolormesh(lonsNS*180./np.pi, latsNS*180./np.pi,np.ma.masked_array(mask2020plus, mask2020plus < 1)
                     ,transform=ccrs.PlateCarree()
                     ,cmap =cMap
                     ,alpha=0.5
                     ,linewidth=0, antialiased=True # this will hide the grid
                     )       
        
    ax2.scatter(loncS, latcS
                   ,s=70,edgecolor='k',color='red',alpha=1, zorder=3,marker='s'
                   ,transform=ccrs.PlateCarree())
    ax2.scatter(loncSplus_hr, latcSplus_hr
               ,s=70,edgecolor='k',color='blue',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    ax2.plot(lon_Dun,lat_Dun,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax2.text(lon_Dun+10,lat_Dun-4,"Dunedin",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())


    
    #ax2.text(45,-85, '2020+', color='blue',fontsize=18,transform=ccrs.PlateCarree())
    #ax2.text(-60,-70, '2020', color='red',fontsize=18,transform=ccrs.PlateCarree())
    
    ax2.set_title('South pole',fontsize=20)
    
    plt.savefig(folder_quantity+'polar_zones_dbetal_'+str(l)+'_incr_multiplier_'+str(increase)+'_res_lon='+str(np.rad2deg(lon_res_b))+'deg_bisection_brute.pdf',bbox_inches='tight',pad_inches=0.1)    
    
    plt.show(block=False)        
    
    
    
    #########################
    # to check:
    ## plot zones: only the difficult one is plotted brute force
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
        
    ax1.fill_between(np.rad2deg(lons_b), latsNplus[:,0], latsNplus[:,1], 
         alpha=0.5,color='blue', linewidth =0,
         transform=ccrs.PlateCarree())
    
    ax1.scatter(loncN, latcN
                   ,s=70,edgecolor='k',color='red',alpha=1, zorder=3,marker='s'
                   ,transform=ccrs.PlateCarree())
        
    ax1.scatter(loncNplus, latcNplus
               ,s=70,edgecolor='k',color='blue',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    ax1.plot(lon_Leeds,lat_Leeds,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Leeds+1.5,lat_Leeds-7,"Leeds",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    ax1.plot(lon_Ed,lat_Ed,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Ed+7,lat_Ed-37,"Edmonton",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())

    ax1.plot(lon_Sal,lat_Sal,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax1.text(lon_Sal+1,lat_Sal-3,"Salekhard",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())


    
    #ax1.text(-50, 40, '2020+', color='blue',fontsize=18,transform=ccrs.PlateCarree())
    #ax1.text(95, 68, '2020', color='red',fontsize=18,transform=ccrs.PlateCarree())
    
    ax1.set_title('North pole',fontsize=20)
    
    
    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.Orthographic(0, -90))
    
    ax2.set_global() # important (apparently)
    ax2.coastlines()
    ax2.gridlines(crs=ccrs.PlateCarree(), 
                      linewidth=1, color='black', linestyle=':')    

    ax2.fill_between(np.rad2deg(lons_b), latsS[:,0], latsS[:,1], 
         alpha=0.5,color='red', linewidth =0,
         transform=ccrs.PlateCarree())       
    
    cMap = c.ListedColormap(['blue'])
    ax2.pcolormesh(lonsNS*180./np.pi, latsNS*180./np.pi,np.ma.masked_array(mask2020plus, mask2020plus < 1)
                     ,transform=ccrs.PlateCarree()
                     ,cmap =cMap
                     ,alpha=0.5
                     ,linewidth=0, antialiased=True # this will hide the grid
                     )       
    
    ax2.plot(coordsSplus_pol[:,0], coordsSplus_pol[:,1],
         #alpha=0.5,
         color='blue', 
         #linewidth =0,
         transform=ccrs.PlateCarree())
    
    #ax2.plot(coordsSplus_eq[:,0], coordsSplus_eq[:,1],
    ax2.plot(np.rad2deg(lons_b), latsSplus[:,1],    
         #alpha=0.5,
         color='blue', 
         #linewidth =0,
         transform=ccrs.PlateCarree())
        
    ax2.scatter(loncS, latcS
                   ,s=70,edgecolor='k',color='red',alpha=1, zorder=3,marker='s'
                   ,transform=ccrs.PlateCarree())
    ax2.scatter(loncSplus_hr, latcSplus_hr
               ,s=70,edgecolor='k',color='blue',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    ax2.plot(lon_Dun,lat_Dun,'o',color='teal',ms =MS,transform=ccrs.PlateCarree())
    ax2.text(lon_Dun+10,lat_Dun-4,"Dunedin",fontsize=FSL,color='teal',transform=ccrs.PlateCarree())


    
    #ax2.text(45,-85, '2020+', color='blue',fontsize=18,transform=ccrs.PlateCarree())
    #ax2.text(-60,-70, '2020', color='red',fontsize=18,transform=ccrs.PlateCarree())
    
    ax2.set_title('South pole',fontsize=20)
    
    #plt.savefig(folder_quantity+'polar_zones_dbetal_'+str(l)+'_incr_multiplier_'+str(increase)+'_res_lon='+str(np.rad2deg(lon_res_b))+'deg_bisection_brute.pdf',bbox_inches='tight',pad_inches=0.1)    
    
    plt.show(block=False)        
    





    print('****************************************') 
    print('Let''s check:')   
    print('')
    print('Degree changed, l = '+str(l))   
    print('Actual annual SV multiplied by '+str(increase)+' times')   
    print('')
    print('A_N = '+str(areaN))
    print('A_S = '+str(areaS))
    print('(A_plus - A)_N, l = '+str(areaNplus-areaN)+'(rel change = '+str((areaNplus-areaN)/areaN)+')')   
    print('(A_plus - A)_S, l = '+str(areaSplus_hr-areaS)+'(rel change = '+str((areaSplus_hr-areaS)/areaS)+')')   
    print('')     
    print('for cities ')
    print(names_loc)
    print('dist =')
    print(dist_min)
    print('(dist_plus - dist), l = ')
    print(dist_min_plus-dist_min)
    print('relative change')
    print((dist_min_plus-dist_min)/dist_min)
    print('')     
    print('centroid latitudes ')
    print('lat_N = '+str(latcN))
    print('(lat_plus - lat)_N, l = '+str(latcNplus-latcN)+'(rel change = '+str((latcNplus-latcN)/latcN)+')')   
    print('lat_S = '+str(latcS))
    print('(lat_plus - lat)_S, l = '+str(latcSplus_hr-latcS)+'(rel change = '+str((latcSplus_hr-latcS)/latcS)+')')   
    print('****************************************')            
          

else:
   print("Something went wrong")