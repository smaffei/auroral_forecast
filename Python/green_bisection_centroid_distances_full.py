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


# old notes from previous script
# for fixed rel_increase = 0.001 I had probelsm (too small variation) for h_6^6 
#  and, checking, for h_13^13. Seemed like for l=m bigger variations are needed.
#  I then tested exponential increase on h_10^10, with seemently ok results

# exponential increase = math.exp(3.7206+0.296276*l) * [-1,1] seemed fine for all but g_8^8
#  which I then looked at separately and replaced in the values of G_N/S_plus/minus for math.exp(3.7206+0.296276*l) * [-0.1,0.1]
#  STILL NEED TO CHECK THE SURFACE ARE CHANGE
#  The reason must be related to sectorial harmonics not changin much the field at high latitudes.
# similar problem for g_13^13, for which I took math.exp(3.7206+0.296276*l) * [-10,10]
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
#Yakutsk
lat_Yak=62.0272
lon_Yak=129.7319

#names_loc = ["Quebec City","Leeds","Dunedin","Edmonton","Salekhard"]
#lats_loc = np.array([lat_QC,lat_Leeds,lat_Dun,lat_Ed,lat_Sal])
#lons_loc = np.array([lon_QC,lon_Leeds,lon_Dun,lon_Ed,lon_Sal])

#names_loc = ["Quebec City","Leeds","Dunedin","Edmonton"]
#lats_loc = np.array([lat_QC,lat_Leeds,lat_Dun,lat_Ed])
#lons_loc = np.array([lon_QC,lon_Leeds,lon_Dun,lon_Ed])

names_loc = ["Leeds","Dunedin","Edmonton","Salekhard"]
lats_loc = np.array([lat_Leeds,lat_Dun,lat_Ed,lat_Sal])
lons_loc = np.array([lon_Leeds,lon_Dun,lon_Ed,lon_Sal])

#names_loc = ["Salekhard","Yakutsk"]
#lats_loc = np.array([lat_Sal, lat_Yak])
#lons_loc = np.array([lon_Sal, lon_Yak])


lats_loc = np.reshape(lats_loc, (1, lats_loc.shape[0]))
lons_loc = np.reshape(lons_loc, (1, lons_loc.shape[0]))


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
# for dg = rel_increase*rms_g : not always stable
#rel_increase = np.array([-0.001, 0.001])
# for exponential dg inecrease: more stable
rel_increase = np.array([-1, 1])



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

    
    '''
    # for the full magmodel as background
    folder3   = '/Users/stefanomaffei/Python/aacgmv2/aacgmv2/'
    filename3 = 'magmodel_1590-2020.txt'
    '''
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
    t2020 = 2020
    t1970 = 1970
    dtime2020 = dt.datetime(t2020, 1, 1)        
    
    # prepare bisecting longitude. If the reference zones are already calculated, load them        
    filename_latsN = folder_quantity+"latsN_"+str(t2020)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3
    filename_latsS = folder_quantity+"latsS_"+str(t2020)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3

    if os.path.exists(filename_latsN) and os.path.exists(filename_latsS):
        lons_b = np.loadtxt(filename_latsN)[:,0]*np.pi/180
    else:
        lons_b = np.linspace(0,2*np.pi,num=int(2*np.pi/lon_res_b))
    
    # lats and lons for plotting the coordinates
    '''
    # uniform grid:
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
    
    # for the way I set up things the SVcoeffs might not make sense. Let's prepare explicit SV arrays:
    # magmodel 2020-2025 forecast
    _,_,_,_,_,SVcoeffs_igrf_20_25 = aacgm_functions.read_magmodel_files(folder3+filename3)
    SVcoeffs_igrf_20_25_mat = SH_library.lin2matCoeffs(SVcoeffs_igrf_20_25)

    
    # current model SV
    SVcoeffs3 =  (MFcoeffs[idx_time,1:] - MFcoeffs[idx_time-1,1:])/(MFcoeffs[idx_time,0] - MFcoeffs[idx_time-1,0])
    
    ########################################################
    # Green's function approach for the temporal variations   
    #######################################################    

    ##############################
    # csv file name with results
    # name for fixed rel_increase and dg = rel_increase*rms_g
    #csv_name = folder_quantity + 'full_green_results_target_res_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_rel_increase_'+str(rel_increase[1])+'.csv'
    # name for exponential increase 
    csv_centroid_name = folder_quantity + 'full_green_centroid_results_target_res_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_.csv'
    csv_dist_name = folder_quantity + 'full_green_cities_distance_results_target_res_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_.csv'
    
    ##############################
    # initialise stuff
    
    coeff_cs = [] # 'g' or 'h'
    coeff_l  = np.array([])
    coeff_m  = np.array([])    
    
    G_N_centroid = np.array([])
    G_S_centroid = np.array([])
    G_N_centroid_lon = np.array([])
    G_S_centroid_lon = np.array([])
    G_N_centroid_lat = np.array([])
    G_S_centroid_lat = np.array([])
    
    G_dist  =   np.array([])
    
    inotE = -1 # index of non-empty arrays

    dg =   np.array([])
    g  =   np.array([])
    
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
        os.system('rm '+folder_out+'bisection_'+str(latNpol_target)+'_lats_'+filename3)
        os.system('rm '+folder_out+'bisection_'+str(latNeq_target)+'_lats_'+filename3)
        
        
        latsS[:,0] = np.loadtxt(folder_out+'bisection_'+str(latSpol_target)+'_lats_'+filename3)
        latsS[:,1] = np.loadtxt(folder_out+'bisection_'+str(latSeq_target)+'_lats_'+filename3)
        
        np.savetxt(filename_latsS,np.hstack([lons_b[:,None]*180/np.pi, latsS]),header='lons lats(polar edge) lats(equatorward edge)')
        os.system('rm '+folder_out+'bisection_'+str(latSpol_target)+'_lats_'+filename3)
        os.system('rm '+folder_out+'bisection_'+str(latSeq_target)+'_lats_'+filename3)
        
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
    
    
    for l in range(1,Lmax+1,1):
    #for l in range(4,Lmax+1,1):
        for m in range(0,l+1,1):
        #for m in range(0,1,1):
            for cs in ['g','h']:
                if m==0:
                    if cs == 'h':
                        continue
                
                
                print('')    
                print('*******************************************************')    
                print('calculation for coefficient: '+cs+'_'+str(l)+'^'+str(m))
                print('*******************************************************')
                
                # name the folder for this coefficent
                folder_coeff = folder_quantity+cs+'_'+str(l)+'_'+str(m)+'/'
                
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
                
                ##################################
                # load panda dataframe in which to save results
                print('loading the .csv file to check for previously calculated results')
                rows_data =[]
                if os.path.exists(csv_centroid_name):
                    df_in_centroid = pd.read_csv(csv_centroid_name, index_col=0)
                if os.path.exists(csv_dist_name):
                    df_in_dist = pd.read_csv(csv_dist_name, index_col=0)
                    
                for ir in range(len(rel_increase)):
                    # define increment
                    ## use the rms of the coefficients at this l to define the increments
                    #rms_g = np.sqrt( sum( MFcoeffs[-1,idx:int( sum(coeffs_in_l(np.linspace(1,l,l))) +1)]**2 ) ) 
                    #increase = rms_g * rel_increase[ir]
                    ## use the experimental increase law I found in the Mathematica script:
                    increase = math.exp(3.7206+0.296276*l) * rel_increase[ir]
                    
                    # less sophisticated than the single case script: if file exists, load it and do not calculate anything
                    if os.path.exists(csv_centroid_name):
                        idx_E = np.where( (df_in_centroid['g/h']== cs) 
                                           & (df_in_centroid['l'].to_numpy().astype(int)== int(l)) 
                                           & (df_in_centroid['m'].to_numpy().astype(int)== int(m)) 
                                           #& (df_in['dg[nT]'].to_numpy().round(4) == increase.round(4))  # not for exponential-increase, since I manually corrected g_8^8
                                           & (df_in_centroid['res_lon[deg]'].to_numpy().round(4) == np.rad2deg(lon_res_b).round(4)) 
                                           )  
                        if np.size(idx_E): # if index array is not empty
                            print('case exists')
                            continue # skip current for loop iteration
                    
                    print('case does not exist, calculating:')
                    # if the calculation goes on, update this index
                    inotE = inotE+1
                    
                    # identify the coefficient
                    coeff_cs.append(cs)
                    coeff_l  = np.append(coeff_l,l)
                    coeff_m  = np.append(coeff_m,m)
                    
                    # and append values to the arrays
                    
                    G_N_centroid = np.append(G_N_centroid,0)
                    G_S_centroid = np.append(G_S_centroid,0)
                    G_N_centroid_lon = np.append(G_N_centroid_lon,0)
                    G_S_centroid_lon = np.append(G_S_centroid_lon,0)
                    G_N_centroid_lat = np.append(G_N_centroid_lat,0)
                    G_S_centroid_lat = np.append(G_S_centroid_lat,0)

                    if len(G_dist)==0:
                        G_dist = np.zeros(lons_loc.shape)
                    else:
                        G_dist  = np.append(G_dist,np.zeros(lons_loc.shape),axis=0)
                
                    dg = np.append(dg,increase)
                    g  = np.append(g,MFcoeffs[-1,idx]+dg[inotE])
                                        
                    ####################################################
                    # create mymodel where I only modify one coefficient
                    # modify one coefficient
                    print('perturbating coefficient')

                    MFcoeffs_plus = 1*MFcoeffs
                    MFcoeffs_plus[-1,idx] = MFcoeffs_plus[-1,idx]+dg[inotE]
                    # save mymodel
                    aacgm_functions.write_magmodel_files(header1, header2, header3, header4, MFcoeffs_plus, SVcoeffs, folder_quantity+'mymodel.txt')
                    
                    
                    #############################
                    # mymodel auroral zones
                    filename_latsNplus = folder_coeff+"latsNplus_"+str(t2020)+"_dg_"+str(dg[inotE])+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3
                    filename_latsSplus = folder_coeff+"latsSplus_"+str(t2020)+"_dg_"+str(dg[inotE])+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3
                
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
                        
                        np.savetxt(filename_latsNplus,np.hstack([lons_b[:,None]*180/np.pi, latsNplus]),header='lons lats(polar edge) lats(equatorward edge)')
                        os.system('rm '+folder_out+'bisection_'+str(latNpol_target)+'_lats_'+'mymodel.txt')
                        os.system('rm '+folder_out+'bisection_'+str(latNeq_target)+'_lats_'+'mymodel.txt')
                        
                        latsSplus[:,0] = np.loadtxt(folder_out+'bisection_'+str(latSpol_target)+'_lats_'+'mymodel.txt')
                        latsSplus[:,1] = np.loadtxt(folder_out+'bisection_'+str(latSeq_target)+'_lats_'+'mymodel.txt')

                        np.savetxt(filename_latsSplus,np.hstack([lons_b[:,None]*180/np.pi, latsSplus]),header='lons lats(polar edge) lats(equatorward edge)')
                        os.system('rm '+folder_out+'bisection_'+str(latSpol_target)+'_lats_'+'mymodel.txt')
                        os.system('rm '+folder_out+'bisection_'+str(latSeq_target)+'_lats_'+'mymodel.txt')
                    
                    

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
                    
                    ## equatorward edge distance from selected locations
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
                    
              
                    #############################################
                    #green's function calculation
                    
                    G_N_centroid[inotE]  =   aacgm_functions.angular_distance(loncN,latcN,loncNplus,latcNplus)  / dg[inotE] 
                    G_S_centroid[inotE]  =   aacgm_functions.angular_distance(loncS,latcS,loncSplus,latcSplus)  / dg[inotE] 
                    G_N_centroid_lat[inotE]  =   (latcNplus - latcN) / dg[inotE] 
                    G_S_centroid_lat[inotE]  =   (latcSplus - latcS)  / dg[inotE] 
                    G_N_centroid_lon[inotE]  =   (loncNplus - loncN)  / dg[inotE] 
                    G_S_centroid_lon[inotE]  =   (loncSplus - loncS)  / dg[inotE] 
                
                    G_dist[inotE,:]  =   (dist_min_plus-dist_min)  / dg[inotE] 

                    toc = time.time()
                    elapsed_time = toc-tic
                    print('elapsed time for current Gauss coefficient: '+str(elapsed_time)+' seconds')
                # end of rel_increase
                
                
    ##################
    # save some output
    print('')
    print('* saving output of calculation * ')
    
    if G_N_centroid.size != 0: # check that not all calculations have been skipped
        data = {'g/h': coeff_cs,
                'l': coeff_l,
                'm': coeff_m,
                'res_lon[deg]':  np.rad2deg(lon_res_b),
                'dg[nT]': dg,
                'g[nT]': g,
                'loncN[deg]':loncN,
                'latcN[deg]':latcN,
                'loncS[deg]':loncS,
                'latcS[deg]':latcS,
                'G_N_centroid[deg/nT]':G_N_centroid,
                'G_S_centroid[deg/nT]':G_S_centroid,
                'G_N_centroid_lon[deg/nT]':G_N_centroid_lon,
                'G_S_centroid_lon[deg/nT]':G_S_centroid_lon,
                'G_N_centroid_lat[deg/nT]':G_N_centroid_lat,
                'G_S_centroid_lat[deg/nT]':G_S_centroid_lat
               }
        
        df = pd.DataFrame (data, columns = ['g/h',
                                            'l',
                                            'm',
                                            'res_lon[deg]',
                                            'dg[nT]',
                                            'g[nT]',
                                            'loncN[deg]',
                                            'latcN[deg]',
                                            'loncS[deg]',
                                            'latcS[deg]',
                                            'G_N_centroid[deg/nT]',
                                            'G_S_centroid[deg/nT]',
                                            'G_N_centroid_lon[deg/nT]',
                                            'G_S_centroid_lon[deg/nT]',
                                            'G_N_centroid_lat[deg/nT]',
                                            'G_S_centroid_lat[deg/nT]'])
        
        # probably not the smartes way to go about it. Shold try and implement the list approach for each iteration of the for loop
        # https://stackoverflow.com/questions/10715965/create-pandas-dataframe-by-appending-one-row-at-a-time
        if os.path.exists(csv_centroid_name):
            df_out = df_in_centroid.append(df, ignore_index=True)
        else:
            df_out = df
            
        df_out.reset_index() # because index may be jumbled when appending
        df_out.to_csv(csv_centroid_name) 
           
                    
    if G_dist.size != 0: # check that not all calculations have been skipped
        data = {'g/h': coeff_cs,
                'l': coeff_l,
                'm': coeff_m,
                'res_lon[deg]':  np.rad2deg(lon_res_b),
                'dg[nT]': dg,
                'g[nT]': g                
                }
        
        df = pd.DataFrame (data, columns = ['g/h',
                                            'l',
                                            'm',
                                            'res_lon[deg]',
                                            'dg[nT]',
                                            'g[nT]'])
        # not sure this is correct...
        for loc in range(len(names_loc)):
            data_loc = {'dist_min, '+names_loc[loc]:dist_min[0,loc],
                        'G_dist[deg/nT], '+names_loc[loc]:G_dist[:,loc]}
            df_loc = pd.DataFrame (data_loc, columns = ['dist_min, '+names_loc[loc],
                                                        'G_dist[deg/nT], '+names_loc[loc]])
            df = df.join(df_loc)
            
        # probably not the smartes way to go about it. Shold try and implement the list approach for each iteration of the for loop
        # https://stackoverflow.com/questions/10715965/create-pandas-dataframe-by-appending-one-row-at-a-time
        if os.path.exists(csv_dist_name):
            df_out = df_in_dist.append(df, ignore_index=True)
        else:
            df_out = df
            
        df_out.reset_index() # because index may be jumbled when appending
        df_out.to_csv(csv_dist_name) 
           
                
                       
                    
    ###########
    # some plots
    print('* plotting * ')

    #reload csv for plotting
    df_in_centroid = pd.read_csv(csv_centroid_name, index_col=0)
    df_in_dist = pd.read_csv(csv_dist_name, index_col=0)
    
    df_plot_centroid = df_in_centroid[ (df_in_centroid["res_lon[deg]"].to_numpy().round(8)==np.rad2deg(lon_res_b).round(8)) ]
    df_plot_centroid = df_plot_centroid.sort_values(by='dg[nT]')

    df_plot_dist = df_in_dist[ (df_in_dist["res_lon[deg]"].to_numpy().round(8)==np.rad2deg(lon_res_b).round(8)) ]
    df_plot_dist = df_plot_dist.sort_values(by='dg[nT]')
    
    
    '''
    # reference areas  (they should all be the same)
    area_zone_N = df_plot['areaN[km^2]'][0]
    area_zone_S = df_plot['areaS[km^2]'][0]
    '''
    
    
    # background magnetic field
    MFcoeffs_lin = MFcoeffs[MFcoeffs[:,0]==t2020][0,1:]
    MFcoeffs_mat = SH_library.lin2matCoeffs(MFcoeffs_lin)    
    Br, Bt, Bp=SH_library.calcB(MFcoeffs_mat,theta,lons,r_a,r_a)
    F = np.sqrt(Br**2+Bt**2+Bp**2)

    

    # centroid latitude green's functions: matrix with positive values of dg
    G_N_centroid_lat_dgplus_matrix = np.zeros((Lmax,2*Lmax+1))
    G_N_centroid_lat_dgplus_matrix[:,:] = math.nan

    G_S_centroid_lat_dgplus_matrix = np.zeros((Lmax,2*Lmax+1))
    G_S_centroid_lat_dgplus_matrix[:,:] = math.nan
    
    G_dist_dgplus_matrix = np.zeros((Lmax,2*Lmax+1,len(names_loc)))
    G_dist_dgplus_matrix[:,:] = math.nan
    
    ## Green's function in Br (CMB)
    dlatcdBc_N = np.zeros(theta.shape)
    dlatcdBc_S = np.zeros(theta.shape)
    
    ddistdBc = np.zeros((theta.shape[0],theta.shape[1],len(names_loc)))
    
    ## to obtain the actual variations:
    SVcoeffs_mat = SH_library.lin2matCoeffs(SVcoeffs3)
    # d(surf)/dg * SV = d(surf)/dt
    dlatcdt_N = np.copy(G_N_centroid_lat_dgplus_matrix)
    dlatcdt_S = np.copy(G_S_centroid_lat_dgplus_matrix)

    ddistdt = np.copy(G_dist_dgplus_matrix)


    ## Green's function in Br (surf)
    dlatcdBr_N = np.zeros(theta.shape)
    dlatcdBr_S = np.zeros(theta.shape)
    
    dlatcdF_N = np.zeros(theta.shape)

    ddistdBr = np.zeros((theta.shape[0],theta.shape[1],len(names_loc)))
    
    
    ## to obtain the actual variations:
    SVcoeffs_mat = SH_library.lin2matCoeffs(SVcoeffs3)
    # d(surf)/dg * SV = d(surf)/dt
    dlatcdt_N = np.copy(G_N_centroid_lat_dgplus_matrix)
    dlatcdt_S = np.copy(G_S_centroid_lat_dgplus_matrix)

    ddistdt = np.copy(G_dist_dgplus_matrix)
    
    

    latc_STD_N = np.zeros((Lmax,2*Lmax+1))
    latc_STD_N[:,:] = math.nan    
    latc_STD_S = np.copy(latc_STD_N)
    
    latc_DIFF_N = np.zeros((Lmax,2*Lmax+1))
    latc_DIFF_N[:,:] = math.nan    
    latc_DIFF_S = np.copy(latc_DIFF_N)


    dist_STD = np.zeros((Lmax,2*Lmax+1,len(names_loc)))
    dist_STD[:,:] = math.nan    
    
    dist_DIFF = np.zeros((Lmax,2*Lmax+1,len(names_loc)))
    dist_DIFF[:,:] = math.nan  
    
    
    
    #df_plus_centroid = df_plot_centroid[df_plot_centroid['dg[nT]']>0]
    #df_plus_dist = df_plot_dist[df_plot_dist['dg[nT]']>0]
    
    # fill it with positive variations
    for l in range(1,Lmax+1,1):
        df_l_centroid = df_plot_centroid[df_plot_centroid['l']==l]
        df_l_dist = df_plot_dist[df_plot_dist['l']==l]
        
        SVcoeffs_l = SVcoeffs_mat[SVcoeffs_mat[:,0]==l]
        
        
        for m in range(0,l+1,1):
            df_m_centroid = df_l_centroid[df_l_centroid['m']==m]
            df_m_dist = df_l_dist[df_l_dist['m']==m]

            SVcoeffs_m = SVcoeffs_l[SVcoeffs_l[:,1]==m]
            

            for cs in ['g','h']:
                if cs =='g':
                    df_cs_centroid = df_m_centroid[df_m_centroid['g/h']=='g']
                    df_cs_dist = df_m_dist[df_m_dist['g/h']=='g']

                    df_plus_centroid = df_cs_centroid[df_cs_centroid['dg[nT]']>0]
                    df_plus_dist = df_cs_dist[df_cs_dist['dg[nT]']>0]
                    
                    SV = SVcoeffs_m[0,2]
                                        
                    G_N_centroid_lat_dgplus_matrix[l-1,Lmax+m] = df_cs_centroid['G_N_centroid_lat[deg/nT]'].mean()
                    G_S_centroid_lat_dgplus_matrix[l-1,Lmax+m] = df_cs_centroid['G_S_centroid_lat[deg/nT]'].mean()

                    for loc in range(len(names_loc)):
                        G_dist_dgplus_matrix[l-1,Lmax+m,loc] = df_cs_dist['G_dist[deg/nT], '+names_loc[loc]].mean()
                    
                    # d(surf)/dg * SV = d(surf)/dt
                    dlatcdt_N[l-1,Lmax+m] = G_N_centroid_lat_dgplus_matrix[l-1,Lmax+m]*SV
                    dlatcdt_S[l-1,Lmax+m] = G_S_centroid_lat_dgplus_matrix[l-1,Lmax+m]*SV

                    ddistdt[l-1,Lmax+m,:] = G_dist_dgplus_matrix[l-1,Lmax+m,:]*SV

                    
                    # Transform to Green's function in terms of Br at the CMB
                    SH = SH_library.SchmidtSH(l,m,theta,lons,'c')
                    dSHdth = SH_library.DthSchmidtSH(l,m,theta,lons,'c')
                    dSHdph = SH_library.DphSchmidtSH(l,m,theta,lons,'c')
                    
                    norm = 4*np.pi/(2*l+1)
                    
                    dgdBc_map =  (1/(l+1)) *  (r_cmb/r_a)**(l+2) * SH / norm
                    
                    
                    dlatcdBc_N = dlatcdBc_N + G_N_centroid_lat_dgplus_matrix[l-1,Lmax+m]*dgdBc_map
                    dlatcdBc_S = dlatcdBc_S + G_S_centroid_lat_dgplus_matrix[l-1,Lmax+m]*dgdBc_map
                    for loc in range(len(names_loc)):
                        ddistdBc[:,:,loc] = ddistdBc[:,:,loc] + G_dist_dgplus_matrix[l-1,Lmax+m,loc]*dgdBc_map

                    dgdBr_map =  (1/(l+1)) *  (r_a/r_a)**(l+2) * SH / norm
                    
                    
                    dlatcdBr_N = dlatcdBr_N + G_N_centroid_lat_dgplus_matrix[l-1,Lmax+m]*dgdBr_map
                    dlatcdBr_S = dlatcdBr_S + G_S_centroid_lat_dgplus_matrix[l-1,Lmax+m]*dgdBr_map
                    for loc in range(len(names_loc)):
                        ddistdBr[:,:,loc] = ddistdBr[:,:,loc] + G_dist_dgplus_matrix[l-1,Lmax+m,loc]*dgdBr_map
                    
                       
                else:
                    if m==0:
                        continue
                    else:
                        df_cs_centroid = df_m_centroid[df_m_centroid['g/h']=='h']
                        df_cs_dist = df_m_dist[df_m_dist['g/h']=='h']

                        df_plus_centroid = df_cs_centroid[df_cs_centroid['dg[nT]']>0]
                        df_plus_dist = df_cs_dist[df_cs_dist['dg[nT]']>0]
                        
                        SV = SVcoeffs_m[0,3]
                        
                        G_N_centroid_lat_dgplus_matrix[l-1,Lmax-m] = df_cs_centroid['G_N_centroid_lat[deg/nT]'].mean()
                        G_S_centroid_lat_dgplus_matrix[l-1,Lmax-m] = df_cs_centroid['G_S_centroid_lat[deg/nT]'].mean()
                        
                        for loc in range(len(names_loc)):
                            G_dist_dgplus_matrix[l-1,Lmax-m,loc] = df_cs_dist['G_dist[deg/nT], '+names_loc[loc]].mean()
                    
                        # d(surf)/dg * SV = d(surf)/dt
                        dlatcdt_N[l-1,Lmax-m] = G_N_centroid_lat_dgplus_matrix[l-1,Lmax-m]*SV
                        dlatcdt_S[l-1,Lmax-m] = G_S_centroid_lat_dgplus_matrix[l-1,Lmax-m]*SV
     
                        ddistdt[l-1,Lmax-m,:] = G_dist_dgplus_matrix[l-1,Lmax-m,:]*SV

                        # Transform to Green's function in terms of Br at the CMB
                        SH = SH_library.SchmidtSH(l,m,theta,lons,'s')
                        dSHdth = SH_library.DthSchmidtSH(l,m,theta,lons,'s')
                        dSHdph = SH_library.DphSchmidtSH(l,m,theta,lons,'s')
                        
                        norm = 4*np.pi/(2*l+1)
                        
                        dgdBc_map =  (1/(l+1)) *  (r_cmb/r_a)**(l+2) * SH / norm
                        

                        dlatcdBc_N = dlatcdBc_N + G_N_centroid_lat_dgplus_matrix[l-1,Lmax-m]*dgdBc_map
                        dlatcdBc_S = dlatcdBc_S + G_S_centroid_lat_dgplus_matrix[l-1,Lmax-m]*dgdBc_map
      
                        for loc in range(len(names_loc)):
                            ddistdBc[:,:,loc] = ddistdBc[:,:,loc] + G_dist_dgplus_matrix[l-1,Lmax-m,loc]*dgdBc_map
                        
                        dgdBr_map =  (1/(l+1)) *  (r_a/r_a)**(l+2) * SH / norm
                        

                        dlatcdBr_N = dlatcdBr_N + G_N_centroid_lat_dgplus_matrix[l-1,Lmax-m]*dgdBr_map
                        dlatcdBr_S = dlatcdBr_S + G_S_centroid_lat_dgplus_matrix[l-1,Lmax-m]*dgdBr_map
      
                        for loc in range(len(names_loc)):
                            ddistdBr[:,:,loc] = ddistdBr[:,:,loc] + G_dist_dgplus_matrix[l-1,Lmax-m,loc]*dgdBr_map
                     
     
        
    '''
    # is this good?
    # dFdc
    dlatcdF_N=np.zeros(theta.shape)
    Gr = np.zeros(theta.shape)
    Gth = np.zeros(theta.shape)
    Gph = np.zeros(theta.shape)
    GF = np.zeros(theta.shape)
    
    for Glon in range(len(lons_lin)):
        for Glat in range(len(lats_lin)):
            lon_G = lons_lin[Glon]
            lat_G = lats_lin[Glat]
            colat_G = np.pi/2-lat_G
            Br_G, Bt_G, Bp_G = SH_library.calcB(MFcoeffs_mat,np.array([lon_G*np.pi/180]),np.array([lat_G*np.pi/180]),r_a,r_a)
            F_G = np.sqrt(Br_G**2 + Bt_G**2 + Bp_G**2)
            H_G = np.sqrt(Bt_G**2 + Bp_G**2)
            for ilon in range(len(lons_lin)-1):
                for ilat in range(len(lats_lin)-1):
                    print(str(ilon)+ " " +str(ilat)+ " " +str(Glon)+ " " +str(Glat))
                    dph = abs(lons[ilon+1,ilat]-lons[ilon,ilat])
                    dth = abs(theta[ilon,ilat+1]-theta[ilon,ilat])
                    thav = (theta[ilon,ilat+1]+theta[ilon,ilat])/2
                    
                    
                    # check the result
                    #check= check + dth*dph*np.sin(thav) 
                    # check the green's functions

                    mu = np.cos(theta[ilon,ilat])*np.cos(np.deg2rad(colat_G))+np.sin(theta[ilon,ilat])*np.sin(np.deg2rad(colat_G))*np.cos(np.deg2rad(lon_G)-lons[ilon,ilat])
                    h = r_cmb/r_a
                    R = np.sqrt(r_a**2+r_cmb**2-2*r_a*r_cmb*mu)
                    f = R/r_a
                    
                    dNcdmu = (h/(4*np.pi)) * ( (1-2*h*mu+3*h**2)/f**3 + mu/(f*(f+h-mu)) -1/(1-mu) )
            
                    Gr[ilon,ilat]  =  h**2*(1-h**2)/f**3  /(4*np.pi)
                    Gth[ilon,ilat] = -dNcdmu * ( np.cos(np.deg2rad(colat_G))*np.sin(theta[ilon,ilat])*np.cos(np.deg2rad(lon_G)-lons[ilon,ilat]) - np.sin(np.deg2rad(colat_G))*np.cos(theta[ilon,ilat]) )
                    Gph[ilon,ilat] = dNcdmu * ( np.sin(theta[ilon,ilat])*np.sin(np.deg2rad(lon_G)-lons[ilon,ilat]) )
                
                    GF[ilon,ilat] = (1/F_G)*(Br_G*Gr[ilon,ilat] + Bt_G*Gth[ilon,ilat] + Bp_G*Gph[ilon,ilat])
                
                    dlatcdF_N[Glon,Glat] = dlatcdF_N[Glon,Glat] + (dlatcdBc_N[ilon,ilat]/GF[ilon,ilat])*dth*dph*np.sin(thav)   
                    
            
                
    '''
    ## ceontroid latitude change               
    fig,axs = plt.subplots(2,1,figsize=(8,8))
    FS = 12
    FT = 15
    FST = 18
    
    # find min and max of colorbar:
    Zmax = np.max([np.nanmax(G_N_centroid_lat_dgplus_matrix),np.nanmax(G_S_centroid_lat_dgplus_matrix)])    
    Zmin = np.min([np.nanmin(G_N_centroid_lat_dgplus_matrix),np.nanmin(G_S_centroid_lat_dgplus_matrix)])    
    # manual
    Zmax = 1e-3
    Zmin = -1e-3
    
    # NORTH
    ax = axs[0]
    Z = G_N_centroid_lat_dgplus_matrix 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale = 1, linthresh = 1e-8)
                    #,cmap = 'RdYlBu_r'
                    ,cmap = 'seismic'
                    )

    #clb = plt.colorbar(cf
    #               #,cax = ax
    #                )
    #clb.ax.tick_params(labelsize=FS)
    #clb.set_label(r'[km$^2$/nT]',fontsize=FS)

    ax.xaxis.set_ticks_position('bottom')

    m_ticks = np.arange(0,2*Lmax+1,1) # need the +1...simply obsciene
    m_labels = list(map(str, abs(np.arange(-Lmax,Lmax+1,1))))
    ax.set_xticks(m_ticks)
    ax.set_xticklabels(m_labels)

    l_ticks = np.arange(0,Lmax,1) # need the +1...simply obsciene
    l_labels = list(map(str, np.arange(1,Lmax+1,1)))
    ax.set_yticks(l_ticks)
    ax.set_yticklabels(l_labels)
    
    ax.tick_params(axis='x', which = 'major', labelsize=FS)
    ax.tick_params(axis='y', which = 'major', labelsize=FS)

    ax.set_xlabel(r'$\leftarrow \ h_l^m\qquad$ order, $m$ $\qquad g_l^m \ \rightarrow$',fontsize=FS) 
    ax.set_ylabel(r'degree, $l$',fontsize=FS)

    ax.tick_params(labeltop=True, labelright=True, top=True, right=True)
    
    
    ax.set_title(r'Northern Hemisphere, $\frac{\partial\lambda^c_N}{\partial\beta_l^m}$', fontsize=FT, y=1.13)    

    # SOUTH
    ax = axs[1]
    Z = G_S_centroid_lat_dgplus_matrix 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale =1, linthresh = 1e-8)
                    #,cmap = 'RdYlBu_r'
                    ,cmap = 'seismic'
                    )

    #clb = plt.colorbar(cf
    #               #,cax = ax
    #                )
    #clb.ax.tick_params(labelsize=FS)
    #clb.set_label(r'[km$^2$/nT]',fontsize=FS)

    ax.xaxis.set_ticks_position('bottom')

    m_ticks = np.arange(0,2*Lmax+1,1) # need the +1...simply obsciene
    m_labels = list(map(str, abs(np.arange(-Lmax,Lmax+1,1))))
    ax.set_xticks(m_ticks)
    ax.set_xticklabels(m_labels)

    l_ticks = np.arange(0,Lmax,1) # need the +1...simply obsciene
    l_labels = list(map(str, np.arange(1,Lmax+1,1)))
    ax.set_yticks(l_ticks)
    ax.set_yticklabels(l_labels)
    
    ax.tick_params(axis='x', which = 'major', labelsize=FS)
    ax.tick_params(axis='y', which = 'major', labelsize=FS)

    ax.set_xlabel(r'$\leftarrow \ h_l^m\qquad$ order, $m$ $\qquad g_l^m \ \rightarrow$',fontsize=FS) 
    ax.set_ylabel(r'degree, $l$',fontsize=FS)

    ax.tick_params(labeltop=True, labelright=True, top=True, right=True)
    
    
    ax.set_title(r'Southern Hemisphere, $\frac{\partial\lambda^c_S}{\partial\beta_l^m}$', fontsize=FT, y=1.13)    

    #clb = fig.colorbar(cf, ax=axs.ravel().tolist(), shrink=0.95, pad=0.1)
    cb_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
    clb = fig.colorbar(cf, cax=cb_ax)

    clb.ax.tick_params(labelsize=FS)
    clb.set_label(r'[deg/nT]',fontsize=FS)
    
    plt.suptitle('Sensitivity of centroid latitude change', fontsize=FST)
    plt.tight_layout()
                        
    plt.savefig(folder_quantity + 'G_centroid_lat_dgplus_matrix_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show(block=False)     
    

    ## city distances change               
    fig,axs = plt.subplots(2,2,figsize=(14,8))
    FS = 12
    FT = 15
    FST = 18
    
    # find min and max of colorbar:
    Zmax = np.nanmax(G_dist_dgplus_matrix)   
    Zmin = np.nanmin(G_dist_dgplus_matrix) 
    
    loc = -1
    for px in [0,1]:
        for py in [0,1]:
            loc = loc +1
            ax = axs[px,py]
            ax.set_position([0.05+px*0.45,0.05+(1-py)*0.45,0.34,0.4])

            Z = G_dist_dgplus_matrix[:,:,loc] 
            #Z = G_S_centroid_lat_dgplus_matrix 
        
            cf = ax.matshow(Z
                            ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale = 1, linthresh = 1e-6)
                            #,cmap = 'RdYlBu_r'
                            ,cmap = 'seismic'
                            )
        
            #clb = plt.colorbar(cf
            #               #,cax = ax
            #                )
            #clb.ax.tick_params(labelsize=FS)
            #clb.set_label(r'[km$^2$/nT]',fontsize=FS)
        
            ax.xaxis.set_ticks_position('bottom')
        
            m_ticks = np.arange(0,2*Lmax+1,1) # need the +1...simply obsciene
            m_labels = list(map(str, abs(np.arange(-Lmax,Lmax+1,1))))
            ax.set_xticks(m_ticks)
            ax.set_xticklabels(m_labels)
        
            l_ticks = np.arange(0,Lmax,1) # need the +1...simply obsciene
            l_labels = list(map(str, np.arange(1,Lmax+1,1)))
            ax.set_yticks(l_ticks)
            ax.set_yticklabels(l_labels)
            
            ax.tick_params(axis='x', which = 'major', labelsize=FS)
            ax.tick_params(axis='y', which = 'major', labelsize=FS)
        
            ax.set_xlabel(r'$\leftarrow \ h_l^m\qquad$ order, $m$ $\qquad g_l^m \ \rightarrow$',fontsize=FS) 
            ax.set_ylabel(r'degree, $l$',fontsize=FS)
        
            ax.tick_params(labeltop=True, labelright=True, top=True, right=True)
            
            #ax.set_title(r''+names_loc[loc]+'$\mathcal{D}(\delta\beta_l^m)$', fontsize=FT, y=1.13)    
            ax.set_title(names_loc[loc], fontsize=FT, y=1.13)    
        
        
    
    cb_ax = fig.add_axes([0.9,0.1,0.018,0.8])    
    #cb_ax = plt.subplot(1,3,3) 
    #cb_ax.set_position([0.9,0.1,0.03,0.8])
    
    clb = fig.colorbar(cf, cax=cb_ax)
    
    clb.ax.tick_params(labelsize=FS)
    clb.set_label(r'[deg/nT]',fontsize=FS)
    
    plt.suptitle('Sensitivity of city distance change', fontsize=FST)
    #plt.tight_layout()
                        
    plt.savefig(folder_quantity + 'G_dist_dgplus_matrix_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show(block=False)     
    

    
    LmaxSV = 13
    
    
    ## dlatc/dt               
    fig,axs = plt.subplots(2,1,figsize=(8,8))
    FS = 12
    FT = 15
    FST = 18
    
    
    
    # find min and max of colorbar:
    Zmax = np.max([np.nanmax(dlatcdt_N),np.nanmax(dlatcdt_S)])    
    Zmin = np.min([np.nanmin(dlatcdt_N),np.nanmin(dlatcdt_S)])    
    
    # NORTH
    ax = axs[0]
    Z = dlatcdt_N[:LmaxSV,Lmax-LmaxSV:-(Lmax-LmaxSV)] 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale = 1, linthresh = 1e-8)
                    #,cmap = 'RdYlBu_r'
                    ,cmap = 'seismic'
                    )

    ax.xaxis.set_ticks_position('bottom')

    m_ticks = np.arange(0,2*LmaxSV+1,1) # need the +1...simply obsciene
    m_labels = list(map(str, abs(np.arange(-LmaxSV,LmaxSV+1,1))))
    ax.set_xticks(m_ticks)
    ax.set_xticklabels(m_labels)

    l_ticks = np.arange(0,LmaxSV,1) # need the +1...simply obsciene
    l_labels = list(map(str, np.arange(1,LmaxSV+1,1)))
    ax.set_yticks(l_ticks)
    ax.set_yticklabels(l_labels)
    
    ax.tick_params(axis='x', which = 'major', labelsize=FS)
    ax.tick_params(axis='y', which = 'major', labelsize=FS)

    ax.set_xlabel(r'$\leftarrow \ h_l^m\qquad$ order, $m$ $\qquad g_l^m \ \rightarrow$',fontsize=FS) 
    ax.set_ylabel(r'degree, $l$',fontsize=FS)

    ax.tick_params(labeltop=True, labelright=True, top=True, right=True)
    
    
    ax.set_title(r'Northern Hemisphere, $\frac{\partial \lambda^c_N(\partial\beta_l^m)}{\partial t}$', fontsize=FT, y=1.13)    

    # SOUTH
    ax = axs[1]
    Z = dlatcdt_S[:LmaxSV,Lmax-LmaxSV:-(Lmax-LmaxSV)] 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale =1, linthresh = 1e-8)
                    #,cmap = 'RdYlBu_r'
                    ,cmap = 'seismic'
                    )

    #clb = plt.colorbar(cf
    #               #,cax = ax
    #                )
    #clb.ax.tick_params(labelsize=FS)
    #clb.set_label(r'[km$^2$/nT]',fontsize=FS)

    ax.xaxis.set_ticks_position('bottom')

    m_ticks = np.arange(0,2*LmaxSV+1,1) # need the +1...simply obsciene
    m_labels = list(map(str, abs(np.arange(-LmaxSV,LmaxSV+1,1))))
    ax.set_xticks(m_ticks)
    ax.set_xticklabels(m_labels)

    l_ticks = np.arange(0,LmaxSV,1) # need the +1...simply obsciene
    l_labels = list(map(str, np.arange(1,LmaxSV+1,1)))
    ax.set_yticks(l_ticks)
    ax.set_yticklabels(l_labels)
    
    ax.tick_params(axis='x', which = 'major', labelsize=FS)
    ax.tick_params(axis='y', which = 'major', labelsize=FS)

    ax.set_xlabel(r'$\leftarrow \ h_l^m\qquad$ order, $m$ $\qquad g_l^m \ \rightarrow$',fontsize=FS) 
    ax.set_ylabel(r'degree, $l$',fontsize=FS)

    ax.tick_params(labeltop=True, labelright=True, top=True, right=True)
    
    
    ax.set_title(r'Southern Hemisphere, $\frac{\partial \lambda^c_S(\partial\beta_l^m)}{\partial t}$', fontsize=FT, y=1.15)    

    #clb = fig.colorbar(cf, ax=axs.ravel().tolist(), shrink=0.95, pad=0.1)
    cb_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
    clb = fig.colorbar(cf, cax=cb_ax)

    clb.ax.tick_params(labelsize=FS)
    clb.set_label(r'[deg/year]',fontsize=FS)
    
    plt.suptitle('Ceontroid latitude yearly change', fontsize=FST)
    plt.tight_layout()
                        
    plt.savefig(folder_quantity + 'dlatcdt_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show(block=False)  
    
    

    ## same plot as the last two, but as a function of l only
    
    # dlatcdt
    FS = 20
    fig,ax = plt.subplots(figsize=(8,5))
    ax.set_xlabel(r'$l$',fontsize=FS)
    ax.set_ylabel('deg / year ', fontsize=FS)
    #ax.set_ylim([77,81])
    #ax.set_xlim([1900,2100])
    ax.tick_params('x',labelsize=FS)
    ax.axhline(y=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    
    # This seems broken, fix it sometime else
    #ax.plot(np.arange(1,LmaxSV+1,1),np.nansum(dlatcdt_N[:LmaxSV,Lmax-LmaxSV:-(Lmax-LmaxSV)],axis=1),'o-',color='tab:blue',label=r'North')
    
    #ax.plot(np.arange(1,LmaxSV+1,1),np.nansum(dlatcdt_S[:LmaxSV,Lmax-LmaxSV:-(Lmax-LmaxSV)],axis=1),'o-',color='tab:orange',label=r'South')
    
    
    ax.plot(np.arange(1,LmaxSV+1,1),np.nansum(dlatcdt_N[:LmaxSV,:],axis=1),'o-',color='tab:blue',label=r'North')
    
    ax.plot(np.arange(1,LmaxSV+1,1),np.nansum(dlatcdt_S[:LmaxSV,:],axis=1),'o-',color='tab:orange',label=r'South')
    
    ax.legend(fontsize=15,loc='upper right')
    
    ax.tick_params('y',labelsize=FS)
    #ax.yaxis.get_offset_text().set_fontsize(12)
    #ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    #ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    plt.title(r'Centroid latitude yearly change',fontsize=FS)
    plt.savefig(folder_quantity + 'dlatcdt_l_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()


    # ddistdt
    FS = 20
    fig,ax = plt.subplots(figsize=(8,5))
    ax.set_xlabel(r'$l$',fontsize=FS)
    ax.set_ylabel('deg / year ', fontsize=FS)
    #ax.set_ylim([77,81])
    #ax.set_xlim([1900,2100])
    ax.tick_params('x',labelsize=FS)
    ax.axhline(y=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    
    for loc in range(len(names_loc)):
        # this seems broken, fix it some other time
        #ax.plot(np.arange(1,LmaxSV+1,1),np.nansum(ddistdt[:LmaxSV,Lmax-LmaxSV:-(Lmax-LmaxSV),loc],axis=1),'o-',label=names_loc[loc])
        ax.plot(np.arange(1,LmaxSV+1,1),np.nansum(ddistdt[:LmaxSV,:,loc],axis=1),'o-',label=names_loc[loc])
        
    ax.legend(fontsize=15,loc='upper right')
    
    ax.tick_params('y',labelsize=FS)
    #ax.yaxis.get_offset_text().set_fontsize(12)
    #ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    #ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    plt.title(r'City distance yearly change',fontsize=FS)
    plt.savefig(folder_quantity + 'ddistdt_l_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()


    



    
    ## Green's functions for A_i as a function of Br(CMB)
    MS = 150
    FS = 20
    MS1 = 8
    
    fig = plt.figure(figsize=(10,6))
    
    # manually
    expmax=-5
    Zs = dlatcdBc_N/10**expmax

    Zmax = 5
    Zmin = -Zmax
    
    #Zmax = 10
    #Zmin = -10
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    
    
    cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                   ,cmap='PuOr_r' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    plt.scatter(loncN, latcN
               ,s=MS,edgecolor='k',color='r',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
        
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r'$\partial \lambda^c_N/\partial B_c  \times 10^{-5}$ [deg/nT]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dlatcdBc_N.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    
    
    
    fig = plt.figure(figsize=(10,6))
    
    # manually
    expmax=-5
    Zs = dlatcdBc_S/10**expmax

    Zmax = 5
    Zmin = -Zmax

    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    
    
    cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                   ,cmap='PuOr_r' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    plt.scatter(loncS, latcS
               ,s=MS,edgecolor='k',color='r',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r'$\partial \lambda^c_S/\partial B_c  \times 10^{-5}$ [km$^2$/nT]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dlatcdBc_S.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    



    ## Green's functions for A_i as a function of Br(surf)
    MS = 150
    
    fig = plt.figure(figsize=(10,6))
    
    expmax = -4
    
    Zs = dlatcdBr_N/10**expmax
    #Zmax = max(np.nanmax(abs(dAdBc_N)),np.nanmax(abs(dAdBc_S)))
    Zmax = 6
    Zmin = -Zmax
    #Zmax = 10
    #Zmin = -10
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    
    
    cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                   ,cmap='PuOr_r' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    plt.scatter(loncN, latcN
               ,s=MS,edgecolor='k',color='r',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
        
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r'$\partial \lambda^c_N/\partial B_r  \times 10^{-4}$ [deg/nT]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dlatcdBr_N.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    
    
    fig = plt.figure(figsize=(10,6))
    
    Zs = dlatcdBr_S/10**expmax
    
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    
    
    cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                   ,cmap='PuOr_r' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    plt.scatter(loncS, latcS
               ,s=MS,edgecolor='k',color='r',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r'$\partial \lambda^c_S/\partial B_r  \times 10^{-4}$ [km$^2$/nT]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dlatcdBr_S.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    
    '''
    ### try sensitivity to F
    
    fig = plt.figure(figsize=(10,6))
    Zs = dFdBc
    Zmax = np.nanmax(Zs)
    Zmin = -Zmax
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    
    
    cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                   ,cmap='PuOr_r' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    plt.scatter(loncS, latcS
               ,s=MS,edgecolor='k',color='r',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r'$\partial \lambda_c/\partial B_c$ [km$^2$/nT]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    #plt.savefig(folder_quantity + '/dlatcdBc_S.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    '''

    
    
    ## Green's functions for dist as a function of Br(CMB)
    fig = plt.figure(figsize=(10,6))
    
    
    # we manually adjust for the scientific notation 
       
    expmax = -5

    Z = ddistdBc/10**expmax

    #Zmax = np.nanmax(abs(Z))
    Zmax=6
    Zmin = -Zmax
    
    loc = -1
    for px in [0,1]:
        for py in [0,1]:
            loc = loc +1
            ax = fig.add_subplot(2, 2, loc+1, projection=ccrs.PlateCarree())
            ax.set_position([0.05+px*0.5,0.15+(1-py)*0.4,0.4,0.34])
    
            Zs = Z[:,:,loc]

            ax.set_global() # important (apparently)
            ax.coastlines()
            ax.gridlines()    
            
            cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                           ,cmap='PuOr_r' 
                           ,levels=np.linspace(Zmin,Zmax,41)
                           ,transform=ccrs.PlateCarree()
                           )

            ax.plot(lons_loc[:,loc],lats_loc[:,loc],'o',color='c',ms =MS1,transform=ccrs.PlateCarree())
            #ax.text(lons_loc[:,loc]-1,lats_loc[:,loc]-1,names_loc[loc],fontsize=10,color='brown',va='top',ha='right',transform=ccrs.PlateCarree())
    
            #clb.set_label(r'$\partial \lambda_c/\partial B_c$ [deg/nT]', fontsize = FS)
            #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
            #clb.ax.tick_params(labelsize=FS)
            ax.set_title(names_loc[loc], fontsize=FT, y=1.01)    

    cb_ax = fig.add_axes([0.1,0.1,0.8,0.02])    
    #cb_ax = plt.subplot(1,3,3) 
    #cb_ax.set_position([0.9,0.1,0.03,0.8])
    
    clb = fig.colorbar(cf, cax=cb_ax,orientation='horizontal')

    clb.ax.tick_params(labelsize=FS)
    #clb.set_label(r'$\partial d_j/\partial B_c  \times 10^{'+str(expmax)+'} $ [deg/nT]', fontsize = FS)
    clb.set_label(r'$\partial d_j/\partial B_c  \times 10^{-5} $ [deg/nT]', fontsize = FS)
    #clb.set_label(r'deg/ nT $\times 10^{'+str(expmax)+'} $',fontsize=FS)

    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/ddistdBc.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    
    
    ## Green's functions for dist as a function of Br(CMB)
    fig = plt.figure(figsize=(10,6))
    
    expmax = -3
    #Zmax = np.amax(abs(ddistdBr))/10**-3
    Zmax=1.6
    Zmin = -Zmax
    
    loc = -1
    for px in [0,1]:
        for py in [0,1]:
            loc = loc +1
            ax = fig.add_subplot(2, 2, loc+1, projection=ccrs.PlateCarree())
            ax.set_position([0.05+px*0.5,0.15+(1-py)*0.4,0.4,0.34])
    
            Zs = ddistdBr[:,:,loc]/10**-3

            ax.set_global() # important (apparently)
            ax.coastlines()
            ax.gridlines()    
            
            cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                           ,cmap='PuOr_r' 
                           ,levels=np.linspace(Zmin,Zmax,41)
                           ,transform=ccrs.PlateCarree()
                           )

            ax.plot(lons_loc[:,loc],lats_loc[:,loc],'o',color='c',ms =MS1,transform=ccrs.PlateCarree())
            #ax.text(lons_loc[:,loc]-1,lats_loc[:,loc]-1,names_loc[loc],fontsize=10,color='brown',va='top',ha='right',transform=ccrs.PlateCarree())
    
            #clb.set_label(r'$\partial \lambda_c/\partial B_c$ [deg/nT]', fontsize = FS)
            #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
            #clb.ax.tick_params(labelsize=FS)
            ax.set_title(names_loc[loc], fontsize=FT, y=1.01)    

    cb_ax = fig.add_axes([0.1,0.1,0.8,0.02])    
    #cb_ax = plt.subplot(1,3,3) 
    #cb_ax.set_position([0.9,0.1,0.03,0.8])
    
    clb = fig.colorbar(cf, cax=cb_ax,orientation='horizontal')

    clb.ax.tick_params(labelsize=FS)
    clb.set_label(r'$\partial d_j/\partial B_r  \times 10^{-3}$ [deg/nT]', fontsize = FS)
    
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/ddistdBr.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    

    
    #multiply the above by the actual SV from igrf13
    
    SVcoeffsmat = SH_library.lin2matCoeffs(SVcoeffs)
    dBcdt, _, _=SH_library.calcB(SVcoeffsmat,theta,lons,r_a,r_cmb)
    dBrdt, _, _=SH_library.calcB(SVcoeffsmat,theta,lons,r_a,r_a)
    

    check = 0
    dAdtN_cmb =0
    dAdtS_cmb =0
    dAdtS_eq =0
    dAdtS_N =0
    dAdtS_S =0
    dAdtN_eq =0
    dAdtN_N =0
    dAdtN_S =0
    

    '''    
    for iph in range(len(lons[:,0])-1): #I should technically wrap around here
        for ith in range(len(theta[0,:])-1):
            dph = abs(lons[iph+1,ith]-lons[iph,ith])
            dth = abs(theta[iph,ith+1]-theta[iph,ith])
            thav = (theta[iph,ith+1]+theta[iph,ith])/2
            dAdtN_cmb = dAdtN_cmb + dBcdt[iph,ith]*dAdBc_N[iph,ith]*dth*dph*np.sin(thav)
            dAdtS_cmb = dAdtS_cmb + dBcdt[iph,ith]*dAdBc_S[iph,ith]*dth*dph*np.sin(thav)
            
            if theta[iph,ith] < np.pi/3:
                dAdtS_N = dAdtS_N + dBcdt[iph,ith]*dAdBc_S[iph,ith]*dth*dph*np.sin(thav)
                dAdtN_N = dAdtN_N + dBcdt[iph,ith]*dAdBc_N[iph,ith]*dth*dph*np.sin(thav)
            elif theta[iph,ith] > 2*np.pi/3:
                dAdtS_S = dAdtS_S + dBcdt[iph,ith]*dAdBc_S[iph,ith]*dth*dph*np.sin(thav)
                dAdtN_S = dAdtN_S + dBcdt[iph,ith]*dAdBc_N[iph,ith]*dth*dph*np.sin(thav)
            elif theta[iph,ith] < 2*np.pi/3 and theta[iph,ith] > np.pi/3:
                dAdtS_eq = dAdtS_eq + dBcdt[iph,ith]*dAdBc_S[iph,ith]*dth*dph*np.sin(thav)
                dAdtN_eq = dAdtN_eq + dBcdt[iph,ith]*dAdBc_N[iph,ith]*dth*dph*np.sin(thav)
            check= check + dth*dph*np.sin(thav)
    '''
    
    dlatcdtN_profile =np.zeros(theta[0,:-1].shape)
    dlatcdtS_profile =np.zeros(theta[0,:-1].shape)
    
    ddistdt_profile =np.zeros((theta[0,:-1].shape[0],len(names_loc)))

    dlatcdtN_surf_profile =np.zeros(theta[0,:-1].shape)
    dlatcdtS_surf_profile =np.zeros(theta[0,:-1].shape)
    
    ddistdt_surf_profile =np.zeros((theta[0,:-1].shape[0],len(names_loc)))
    
    for ith in range(len(theta[0,:])-1):
        for iph in range(len(lons[:,0])-1):
            dph = abs(lons[iph+1,ith]-lons[iph,ith])
            thav = (theta[iph,ith+1]+theta[iph,ith])/2  
            dlatcdtN_profile[ith] = dlatcdtN_profile[ith] + dBcdt[iph,ith]*dlatcdBc_N[iph,ith]*dph*np.sin(thav)
            dlatcdtS_profile[ith] = dlatcdtS_profile[ith] + dBcdt[iph,ith]*dlatcdBc_S[iph,ith]*dph*np.sin(thav)
            
            ddistdt_profile[ith,:] = ddistdt_profile[ith,:] + dBcdt[iph,ith]*ddistdBc[iph,ith,:]*dph*np.sin(thav)
    
            dlatcdtN_surf_profile[ith] = dlatcdtN_surf_profile[ith] + dBrdt[iph,ith]*dlatcdBr_N[iph,ith]*dph*np.sin(thav)
            dlatcdtS_surf_profile[ith] = dlatcdtS_surf_profile[ith] + dBrdt[iph,ith]*dlatcdBr_S[iph,ith]*dph*np.sin(thav)
            
            ddistdt_surf_profile[ith,:] = ddistdt_surf_profile[ith,:] + dBrdt[iph,ith]*ddistdBr[iph,ith,:]*dph*np.sin(thav)
        

    dlatcdtN_tot=0
    dlatcdtS_tot=0
    dlatcdtN_surf_tot=0
    dlatcdtS_surf_tot=0
    
    ddist_tot=np.zeros((len(names_loc)))


    # assume uniform phi grid
    dph = abs(lons_lin[1]-lons_lin[0])
    for ilon in range(len(lons_lin)-1):
        for ilat in range(len(theta_lin)):
            # check the actual time diff
            dlatcdtN_tot      = dlatcdtN_tot      + dBcdt[ilon,ilat]*dlatcdBc_N[ilon,ilat]* dph * weights[ilat]
            dlatcdtN_surf_tot = dlatcdtN_surf_tot + dBrdt[ilon,ilat]*dlatcdBr_N[ilon,ilat]* dph * weights[ilat]
            dlatcdtS_tot      = dlatcdtS_tot      + dBcdt[ilon,ilat]*dlatcdBc_S[ilon,ilat]* dph * weights[ilat]
            dlatcdtS_surf_tot = dlatcdtS_surf_tot + dBrdt[ilon,ilat]*dlatcdBr_S[ilon,ilat]* dph * weights[ilat]
            
            ddist_tot         = ddist_tot         + dBcdt[ilon,ilat]*ddistdBc[ilon,ilat,:]* dph * weights[ilat]
            
            
            
        
 
    # cmb
    fig = plt.figure(figsize=(10,6))
    
    #expmax = -1
    
    Zs = dBcdt*dlatcdBc_N
    
    #Zmax = max(np.nanmax(abs(Zs)),np.nanmax(abs(Zs)))
    Zmax=0.3
    #Zmax = 10.5
    Zmin = -Zmax
    #Zmax = 10
    #Zmin = -10
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    
    
    cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                   ,cmap='PuOr_r' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    plt.scatter(loncN, latcN
               ,s=MS,edgecolor='k',color='r',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    #clb.set_label(r'$\partial \lambda^c_N/\partial t   \times 10^{-1}$ (CMB) [deg/yr]', fontsize = FS)
    clb.set_label(r'$\partial \lambda^c_N/\partial t$ (CMB) [deg/yr]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dlatcdt_N.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    
    

    
    fig = plt.figure(figsize=(10,6))
    
    #Zs = dBcdt*dlatcdBc_S/10**expmax
    Zs = dBcdt*dlatcdBc_S
    #Zmax = max(np.nanmax(abs(Zs)),np.nanmax(abs(Zs)))
    #Zmax = 10.5
    Zmin = -Zmax
    #Zmax = 10
    #Zmin = -10
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    
    
    cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                   ,cmap='PuOr_r' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    plt.scatter(loncS, latcS
               ,s=MS,edgecolor='k',color='r',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r'$\partial \lambda^c_S/\partial t$ (CMB) [deg/yr]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dlatcdt_S.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()   
    
    
    fig = plt.figure(figsize=(10,6))
    Zmax = 0.2
    Zmin = -Zmax
    loc = -1
    for px in [0,1]:
        for py in [0,1]:
            loc = loc +1
            ax = fig.add_subplot(2, 2, loc+1, projection=ccrs.PlateCarree())
            ax.set_position([0.05+px*0.5,0.15+(1-py)*0.4,0.4,0.34])
            Zs = dBcdt*ddistdBc[:,:,loc]
            
            ax.set_global() # important (apparently)
            ax.coastlines()
            ax.gridlines()    
            
            cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                           ,cmap='PuOr_r' 
                           ,levels=np.linspace(Zmin,Zmax,41)
                           ,transform=ccrs.PlateCarree()
                           )
            ax.plot(lons_loc[:,loc],lats_loc[:,loc],'o',color='c',ms =MS1,transform=ccrs.PlateCarree())

            #clb.set_label(r'$\partial \lambda_c/\partial B_c$ [deg/nT]', fontsize = FS)
            #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
            #clb.ax.tick_params(labelsize=FS)
            ax.set_title(names_loc[loc], fontsize=FT, y=1.01)    

    cb_ax = fig.add_axes([0.1,0.1,0.8,0.02])    
    #cb_ax = plt.subplot(1,3,3) 
    #cb_ax.set_position([0.9,0.1,0.03,0.8])
    
    clb = fig.colorbar(cf, cax=cb_ax,orientation='horizontal')

    clb.ax.tick_params(labelsize=FS)
    clb.set_label(r'$\partial d_j/\partial t$ (CMB) [deg/year]', fontsize = FS)
    
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/ddistdt.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    
    
    
    # surface
    fig = plt.figure(figsize=(10,6))
    
    expmax=-2
    Zs = dBrdt*dlatcdBr_N/10**expmax
    
    #Zmax = max(np.nanmax(abs(Zs)),np.nanmax(abs(Zs)))
    Zmax=5
    #Zmax = 10.5
    Zmin = -Zmax
    #Zmax = 10
    #Zmin = -10
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    
    
    cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                   ,cmap='PuOr_r' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    plt.scatter(loncN, latcN
               ,s=MS,edgecolor='k',color='r',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r"$\partial \lambda^c_N/\partial t   \times 10^{-2}$ (Earth's surface) [deg/yr]", fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dlatcdt_surf_N.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    

    
    fig = plt.figure(figsize=(10,6))
    Zs = dBrdt*dlatcdBr_S/10**expmax
    #Zmax = max(np.nanmax(abs(Zs)),np.nanmax(abs(Zs)))
    #Zmax = 10.5
    #Zmin = -Zmax
    #Zmax = 10
    #Zmin = -10
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()    
    
    cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                   ,cmap='PuOr_r' 
                   ,levels=np.linspace(Zmin,Zmax,41)
                   ,transform=ccrs.PlateCarree()
                   )
    plt.scatter(loncS, latcS
               ,s=MS,edgecolor='k',color='r',alpha=1, zorder=3,marker='s'
               ,transform=ccrs.PlateCarree())
    
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r"$\partial \lambda^c_S/\partial t  \times 10^{-2}$ (Earth's surface) [deg/yr]", fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dlatcdt_surf_S.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()   
    
    
    fig = plt.figure(figsize=(10,6))
    Zmax = 0.08
    Zmin = -Zmax
    loc = -1
    for px in [0,1]:
        for py in [0,1]:
            loc = loc +1
            ax = fig.add_subplot(2, 2, loc+1, projection=ccrs.PlateCarree())
            ax.set_position([0.05+px*0.5,0.15+(1-py)*0.4,0.4,0.34])
            Zs = dBrdt*ddistdBr[:,:,loc]
            
            ax.set_global() # important (apparently)
            ax.coastlines()
            ax.gridlines()    
            
            cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
                           ,cmap='PuOr_r' 
                           ,levels=np.linspace(Zmin,Zmax,41)
                           ,transform=ccrs.PlateCarree()
                           )
            ax.plot(lons_loc[:,loc],lats_loc[:,loc],'o',color='c',ms =MS1,transform=ccrs.PlateCarree())

            #clb.set_label(r'$\partial \lambda_c/\partial B_c$ [deg/nT]', fontsize = FS)
            #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
            #clb.ax.tick_params(labelsize=FS)
            ax.set_title(names_loc[loc], fontsize=FT, y=1.01)    

    cb_ax = fig.add_axes([0.1,0.1,0.8,0.02])    
    #cb_ax = plt.subplot(1,3,3) 
    #cb_ax.set_position([0.9,0.1,0.03,0.8])
    
    clb = fig.colorbar(cf, cax=cb_ax,orientation='horizontal')

    clb.ax.tick_params(labelsize=FS)
    clb.set_label(r"$\partial d_j/\partial t$ (Earth's surface) [deg/year]", fontsize = FS)
    
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/ddistdt_surf.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    

    
    
    # integrated in longitude
    
    #cmb
    fig,ax = plt.subplots(figsize=(5,8))
    ax.set_ylabel(r'lat [deg]',fontsize=FS)
    ax.set_xlabel('deg / year ', fontsize=FS)
    #ax.set_ylim([77,81])
    #ax.set_xlim([1900,2100])
    ax.tick_params('x',labelsize=FS)
    
    ax.plot(dlatcdtN_profile,90-180*theta[0,:-1]/np.pi,color='tab:blue',label=r'North')
    
    ax.plot(dlatcdtS_profile,90-180*theta[0,:-1]/np.pi,color='tab:orange',label=r'South')
    
    ax.axvline(x=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')

    ax.legend(fontsize=15,loc='upper right')
    
    ax.tick_params('y',labelsize=FS)
    ax.set_xlim((-0.2, 0.2))
    ax.yaxis.get_offset_text().set_fontsize(12)
    #ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    #ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    plt.title(r'Centroid latitude yearly change (CMB)',fontsize=FS)
    plt.savefig(folder_quantity + 'dlatcdt_lats.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()

    # surf
    fig,ax = plt.subplots(figsize=(5,8))
    ax.set_ylabel(r'lat [deg]',fontsize=FS)
    ax.set_xlabel('deg / year ', fontsize=FS)
    #ax.set_ylim([77,81])
    #ax.set_xlim([1900,2100])
    ax.tick_params('x',labelsize=FS)
    
    ax.plot(dlatcdtN_surf_profile,90-180*theta[0,:-1]/np.pi,color='tab:blue',label=r'North')
    
    ax.plot(dlatcdtS_surf_profile,90-180*theta[0,:-1]/np.pi,color='tab:orange',label=r'South')
    
    ax.axvline(x=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')

    ax.legend(fontsize=15,loc='upper right')
    
    ax.tick_params('y',labelsize=FS)
    ax.set_xlim((-0.1, 0.1))
    ax.yaxis.get_offset_text().set_fontsize(12)
    #ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    #ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    plt.title(r"Centroid latitude yearly change (Earth's surface)",fontsize=FS)
    plt.savefig(folder_quantity + 'dlatcdt_surf_lats.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()
    

    #cmb
    fig,ax = plt.subplots(figsize=(5,8))
    ax.set_ylabel(r'lat [deg]',fontsize=FS)
    ax.set_xlabel('deg / year ', fontsize=FS)
    #ax.set_ylim([77,81])
    #ax.set_xlim([1900,2100])
    ax.tick_params('x',labelsize=FS)
    
    for loc in range(len(names_loc)):
        ax.plot(ddistdt_profile[:,loc],90-180*theta[0,:-1]/np.pi,label=names_loc[loc])
    
    ax.axvline(x=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')

    ax.legend(fontsize=15,loc='lower right')
    
    ax.tick_params('y',labelsize=FS)
    ax.set_xlim((-0.07, 0.07))
    ax.yaxis.get_offset_text().set_fontsize(12)
    #ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    #ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    plt.title(r'City distances yearly change',fontsize=FS)
    plt.savefig(folder_quantity + 'ddistdt_lats.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()
    
    
    # surf
    fig,ax = plt.subplots(figsize=(5,8))
    ax.set_ylabel(r'lat [deg]',fontsize=FS)
    ax.set_xlabel('deg / year ', fontsize=FS)
    #ax.set_ylim([77,81])
    #ax.set_xlim([1900,2100])
    ax.tick_params('x',labelsize=FS)
    
    for loc in range(len(names_loc)):
        ax.plot(ddistdt_surf_profile[:,loc],90-180*theta[0,:-1]/np.pi,label=names_loc[loc])
    
    ax.axvline(x=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')

    ax.legend(fontsize=15,loc='lower right')
    
    ax.tick_params('y',labelsize=FS)
    ax.set_xlim((-0.045, 0.045))
    ax.yaxis.get_offset_text().set_fontsize(12)
    #ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    #ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    plt.title(r"City distances yearly change (Earth's surface)",fontsize=FS)
    plt.savefig(folder_quantity + 'ddistdt_surf_lats.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    

    
    ##########
    # cleanup    
    os.system('rm '+folder_out+'CGMlats_'+str(t2020)+'_'+filename3)
    os.system('rm '+folder_out+'CGMlons_'+str(t2020)+'_'+filename3)
    os.system('rm '+folder_out+'CGMmlt_'+str(t2020)+'_'+filename3)
    
    os.system('rm '+folder_quantity+'CGMlats_'+str(t2020)+'_'+'mymodel.txt')
    os.system('rm '+folder_quantity+'CGMlons_'+str(t2020)+'_'+'mymodel.txt')
    os.system('rm '+folder_quantity+'CGMmlt_'+str(t2020)+'_'+'mymodel.txt')

    os.system('rm '+folder_quantity+'mymodel.txt')

    
else:
   print("Something went wrong")