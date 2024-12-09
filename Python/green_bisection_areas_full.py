#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 16:40:01 2021

@author: stefanomaffei

Green's function approach to analyse the sensitivity of the auroral zone shape
to different Gauss coefficient's time variations

Uses a bisection algorithm to find the edges of the auroral zones

Uses an algorithm from Bevis and Cambareri, 1987 (and used in Zossi, Fagre, Amit and Elias, 2020)
to compute areas of the zones

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

from scipy.special import roots_legendre, eval_legendre

import SH_library


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

# main function    
if __name__ == "__main__":   
    
    plt.close('all')

    # for the dipolar magmodel as background
    folder_magmodel   = '../models/magmodel/'
    filename_magmodel = 'magmodel_1590-2020.txt'
    #filename_magmodel = 'magmodel_1590-2020_dipole.txt'
    #filename3 = 'magmodel_1590-2020_dipole.txt'
    #filename3 = 'magmodel_1590-2020_axial_dipole_smallH.txt'    

    # ifor IGRF model: the target year is the year itself
    t_ref = 2020
    t1t = 2020
    
    
    
    folder_base = './'
    folder_out = folder_base+'aacgmv2_coords/'
    folder_models = folder_base+'models/'
    
    folder_green = folder_base+'Green_time_magmodel'+str(t_ref)+'_bisection/'
    #folder_green = folder_base+'Green_time_magmodel2020axial_dipole_bisection/'
    #folder_green = folder_base+'Green_time_magmodel2020dipole_bisection/'
    
    folder_quantity = folder_green + 'auroral_zones/'
    if not os.path.exists(folder_quantity):
        os.makedirs(folder_quantity) 
    
   
    
    folder3 = folder_magmodel
    filename3 = filename_magmodel
    
    dtime_ref = dt.datetime(t1t, 1, 1)        

    # prepare bisecting longitude. If the reference zones are already calculated, load them        
    filename_latsN = folder_quantity+"latsN_"+str(t_ref)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3
    filename_latsS = folder_quantity+"latsS_"+str(t_ref)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3

    if os.path.exists(filename_latsN) and os.path.exists(filename_latsS):
        lons_b = np.loadtxt(filename_latsN)[:,0]*np.pi/180
    else:
        lons_b = np.linspace(0,2*np.pi,num=int(2*np.pi/lon_res_b))
    
    # lats and lons for plotting the coordinates
    '''
    # uniform theta grid
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
    # reference model 
    ################################        
    header1, header2, header3, header4, MFcoeffs, SVcoeffs = aacgm_functions.read_magmodel_files(folder3+filename3)
    
    # find temporal index in array
    idx_time = np.where(MFcoeffs[:,0]==t1t)
    idx_time = idx_time[0][0]
    
    # for the way I set up things the SVcoeffs might not make sense. Let's prepare explicit SV arrays:
    # magmodel 2020-2025 forecast
    _,_,_,_,_,SVcoeffs_igrf_20_25 = aacgm_functions.read_magmodel_files(folder_magmodel+filename_magmodel)
    SVcoeffs_igrf_20_25_mat = SH_library.lin2matCoeffs(SVcoeffs_igrf_20_25)

    
    # current model SV
    SVcoeffs3 =  (MFcoeffs[idx_time,1:] - MFcoeffs[idx_time-1,1:])/(MFcoeffs[idx_time,0] - MFcoeffs[idx_time-1,0])
    
    MFcoeffs_mat = SH_library.lin2matCoeffs(MFcoeffs[idx_time,1:])
    
    ########################################################
    # Green's function approach for the temporal variations   
    #######################################################    

    ##############################
    # csv file name with results
    # name for fixed rel_increase and dg = rel_increase*rms_g
    #csv_name = folder_quantity + 'full_green_results_target_res_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_rel_increase_'+str(rel_increase[1])+'.csv'
    # name for exponential increase 
    csv_name = folder_quantity + 'full_green_results_target_res_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_.csv'
    
    ##############################
    # initialise stuff
    
    coeff_cs = [] # 'g' or 'h'
    coeff_l  = np.array([])
    coeff_m  = np.array([])    
    
    areaN = np.array([])
    areaNplus = np.array([])
    areaNdiff_plus  = np.array([]) # areas covered by modified zone
    areaNdiff_minus = np.array([])  # areas not covered by modified zone
    areaS = np.array([])
    areaSplus = np.array([])
    areaSdiff_plus  = np.array([]) # areas covered by modified zone
    areaSdiff_minus = np.array([])  # areas not covered by modified zone
    
    # green's functions for total area change
    G_areaN_dg = np.array([])
    G_areaS_dg = np.array([])
        
    # green's functions for surface coverage change
    G_N_plus_dg  =   np.array([])
    G_N_minus_dg =   np.array([])

    G_S_plus_dg  =   np.array([])
    G_S_minus_dg =   np.array([])
    
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
        print('calculating reference auroral zone (IGRF)')
        
        # initialise lats arrays: n x 2 where n is range(lons), and the other dimension contains the polward and equatorward boundaries
        latsN = np.zeros((lons_b.shape[0],2))
        latsS = np.zeros((lons_b.shape[0],2))
        
        # northern zone
        # polar edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t1t, 1, 1, np.rad2deg(lons_b), latNpol_target, target_res, latNpol_i, latNeq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        # equatorial edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t1t, 1, 1, np.rad2deg(lons_b), latNeq_target, target_res, latNpol_i, latNeq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        
        # southern zone
        # polar edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t1t, 1, 1, np.rad2deg(lons_b), latSpol_target, target_res, latSpol_i, latSeq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        # equatorial edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t1t, 1, 1, np.rad2deg(lons_b), latSeq_target, target_res, latSpol_i, latSeq_i, max_iterations,folder_out))
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
    
    
    # plot for checking purposes
    
    
    ## plot zones 
    fig = plt.figure(figsize=(15,8))
    
    # north pole
    ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.Orthographic(0, 90))
    
    ax1.set_global() # important (apparently)
    ax1.coastlines()
    ax1.gridlines(crs=ccrs.PlateCarree(), 
                      linewidth=1, color='black', linestyle=':')    
        
    ax1._autoscaleXon = False
    ax1._autoscaleYon = False
    ax1.fill_between(np.rad2deg(lons_b), latsN[:,0], latsN[:,1], 
         alpha=0.5,color='red', linewidth =0,
         transform=ccrs.PlateCarree())
    
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

    #ax2.text(45,-85, '2020+', color='blue',fontsize=18,transform=ccrs.PlateCarree())
    ax2.text(-60,-70, '2020', color='red',fontsize=18,transform=ccrs.PlateCarree())
    
    ax2.set_title('South pole',fontsize=20)
    
    #plt.savefig(folder_case+'polar_zones_g='+str(MFcoeffs_plus[-1,idx])+'_res_lon='+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)    
    
    plt.show(block=False)             
            
    
    
    
    
    for l in range(1,Lmax+1,1):
    #for l in range(1,2,1):
        #for m in range(0,1,1):
        for m in range(0,l+1,1):
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
                if os.path.exists(csv_name):
                    df_in = pd.read_csv(csv_name, index_col=0)
            
                for ir in range(len(rel_increase)):
                    # define increment
                    ## use the rms of the coefficients at this l to define the increments
                    #rms_g = np.sqrt( sum( MFcoeffs[-1,idx:int( sum(coeffs_in_l(np.linspace(1,l,l))) +1)]**2 ) ) 
                    #increase = rms_g * rel_increase[ir]
                    ## use the experimental increase law I found in the Mathematica script:
                    increase = math.exp(3.7206+0.296276*l) * rel_increase[ir]
                    
                    # less sophisticated than the single case script: if file exists, load it and do not calculate anything
                    if os.path.exists(csv_name):
                        idx_E = np.where( (df_in['g/h']== cs) 
                                           & (df_in['l'].to_numpy().astype(int)== int(l)) 
                                           & (df_in['m'].to_numpy().astype(int)== int(m)) 
                                           #& (df_in['dg[nT]'].to_numpy().round(4) == increase.round(4))  # not for exponential-increase, since I manually corrected g_8^8
                                           & (df_in['res_lon[deg]'].to_numpy().round(4) == np.rad2deg(lon_res_b).round(4)) 
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
                    areaN = np.append(areaN,0) # reference area from igrf
                    areaNplus = np.append(areaNplus,0) # updated area from igrf + perturbation
                    areaNdiff_plus = np.append(areaNdiff_plus,0) # areas covered by modified zone
                    areaNdiff_minus = np.append(areaNdiff_minus,0)  # areas not covered by modified zone
                    areaS = np.append(areaS,0)
                    areaSplus = np.append(areaSplus,0)
                    areaSdiff_plus  = np.append(areaSdiff_plus,0) # areas covered by modified zone
                    areaSdiff_minus = np.append(areaSdiff_minus,0)  # areas not covered by modified zone
                    
                    G_areaN_dg = np.append(G_areaN_dg,0)
                    G_areaS_dg = np.append(G_areaS_dg,0)
                
                    G_N_plus_dg  = np.append(G_N_plus_dg,0)
                    G_N_minus_dg = np.append(G_N_minus_dg,0)
                
                    G_S_plus_dg   = np.append(G_S_plus_dg,0)
                    G_S_minus_dg  = np.append(G_S_minus_dg,0)
                
                    dg = np.append(dg,increase)
                    g  = np.append(g,MFcoeffs[idx_time,idx]+dg[inotE])
                                        
                    ####################################################
                    # create mymodel where I only modify one coefficient
                    # modify one coefficient
                    print('perturbating coefficient')

                    MFcoeffs_plus = 1*MFcoeffs
                    MFcoeffs_plus[idx_time,idx] = MFcoeffs_plus[idx_time,idx]+dg[inotE]
                    # save mymodel
                    aacgm_functions.write_magmodel_files(header1, header2, header3, header4, MFcoeffs_plus, SVcoeffs, folder_quantity+'mymodel.txt')
                    
                    
                    #############################
                    # mymodel auroral zones
                    filename_latsNplus = folder_coeff+"latsNplus_"+str(t_ref)+"_dg_"+str(dg[inotE])+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3
                    filename_latsSplus = folder_coeff+"latsSplus_"+str(t_ref)+"_dg_"+str(dg[inotE])+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3
                
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
                        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder_quantity,'mymodel.txt', t1t, 1, 1, np.rad2deg(lons_b), latNpol_target, target_res, latNpol_i, latNeq_i, max_iterations,folder_out))
                        my_process3.start()
                        my_process3.join() 
                        # equatorial edge
                        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder_quantity,'mymodel.txt', t1t, 1, 1, np.rad2deg(lons_b), latNeq_target, target_res, latNpol_i, latNeq_i, max_iterations,folder_out))
                        my_process3.start()
                        my_process3.join() 
                        
                        # southern zone
                        # polar edge
                        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder_quantity,'mymodel.txt', t1t, 1, 1, np.rad2deg(lons_b), latSpol_target, target_res, latSpol_i, latSeq_i, max_iterations,folder_out))
                        my_process3.start()
                        my_process3.join() 
                        # equatorial edge
                        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder_quantity,'mymodel.txt', t1t, 1, 1, np.rad2deg(lons_b), latSeq_target, target_res, latSpol_i, latSeq_i, max_iterations,folder_out))
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
                    
                    
                    
                    #########################################################################
                    ## find intersections of zones' edges: 
                    # assumes small variations (the polar edges never intersects the equatorial ones)
                    print('calculating intersections between the zones')
                    
                    # find intersections of the polward contours of the northern zone
                    lons_x_Npol, lats_x_Npol = aacgm_functions.intersection_bisection(folder3, filename3, 
                                                                                                    folder_quantity, 'mymodel.txt', 
                                                                                                    np.rad2deg(lons_b), latsN[:,0], latsNplus[:,0], 
                                                                                                    latNpol_target, latNpol_target,
                                                                                                    t1t, 1, 1, 
                                                                                                    t1t, 1, 1, 
                                                                                                    target_res, target_res_zero, max_iterations, folder_out)
                    # find intersections of the equatorward contours of the northern zone
                    lons_x_Neq, lats_x_Neq = aacgm_functions.intersection_bisection(folder3, filename3, 
                                                                                                    folder_quantity, 'mymodel.txt', 
                                                                                                    np.rad2deg(lons_b), latsN[:,1], latsNplus[:,1], 
                                                                                                    latNeq_target, latNeq_target,
                                                                                                    t1t, 1, 1, 
                                                                                                    t1t, 1, 1, 
                                                                                                    target_res, target_res_zero, max_iterations, folder_out)        
                    # find intersections of the polward contours of the southern zone
                    lons_x_Spol, lats_x_Spol = aacgm_functions.intersection_bisection(folder3, filename3, 
                                                                                                    folder_quantity, 'mymodel.txt', 
                                                                                                    np.rad2deg(lons_b), latsS[:,0], latsSplus[:,0], 
                                                                                                    latSpol_target, latSpol_target,
                                                                                                    t1t, 1, 1, 
                                                                                                    t1t, 1, 1, 
                                                                                                    target_res, target_res_zero, max_iterations, folder_out)      
                    # find intersections of the equatorward contours of the northern zone
                    lons_x_Seq, lats_x_Seq = aacgm_functions.intersection_bisection(folder3, filename3, 
                                                                                                    folder_quantity, 'mymodel.txt', 
                                                                                                    np.rad2deg(lons_b), latsS[:,1], latsSplus[:,1], 
                                                                                                    latSeq_target, latSeq_target,
                                                                                                    t1t, 1, 1, 
                                                                                                    t1t, 1, 1, 
                                                                                                    target_res, target_res_zero, max_iterations, folder_out)      
                    
            
            
                    ## form the new polygons and calculate areas differences
                    # area will be in km^2
                    print('calculating area and surface coverage differences')
                    
                    areas_diff_Npol, poly_list_Npol = aacgm_functions.area_diff_bisection(np.rad2deg(lons_b), lons_x_Npol, lats_x_Npol, latsN, latsNplus, 0, r_a/1000)
            
                    areas_diff_Neq, poly_list_Neq = aacgm_functions.area_diff_bisection(np.rad2deg(lons_b), lons_x_Neq, lats_x_Neq, latsN, latsNplus, 1, r_a/1000)
            
                    areas_diff_Spol, poly_list_Spol = aacgm_functions.area_diff_bisection(np.rad2deg(lons_b), lons_x_Spol, lats_x_Spol, latsS, latsSplus, 0, r_a/1000)
            
                    areas_diff_Seq, poly_list_Seq = aacgm_functions.area_diff_bisection(np.rad2deg(lons_b), lons_x_Seq, lats_x_Seq, latsS, latsSplus, 1, r_a/1000)

            
                    
                    ## calculate areas differences
            
                    areaNdiff_plus[inotE]  = np.sum(np.concatenate((areas_diff_Npol[areas_diff_Npol>0],areas_diff_Neq[areas_diff_Neq>0])))
            
                    areaNdiff_minus[inotE] = abs(np.sum(np.concatenate((areas_diff_Npol[areas_diff_Npol<0],areas_diff_Neq[areas_diff_Neq<0]))))
            
                    areaSdiff_plus[inotE] = np.sum(np.concatenate((areas_diff_Spol[areas_diff_Spol>0],areas_diff_Seq[areas_diff_Seq>0])))
            
                    areaSdiff_minus[inotE] = abs(np.sum(np.concatenate((areas_diff_Spol[areas_diff_Spol<0],areas_diff_Seq[areas_diff_Seq<0]))))
            
                    ###########################
                    # calculate reference areas
                    areaNpol_b = aacgm_functions.spherical_polygon_area(np.flip(latsN[:,0]),np.flip(np.rad2deg(lons_b)),r_a/1000)
                    areaNeq_b  = aacgm_functions.spherical_polygon_area(np.flip(latsN[:,1]),np.flip(np.rad2deg(lons_b)),r_a/1000)
                    areaN[inotE] = areaNeq_b - areaNpol_b
                    
                    areaSpol_b = aacgm_functions.spherical_polygon_area(latsS[:,0],np.rad2deg(lons_b),r_a/1000)
                    areaSeq_b  = aacgm_functions.spherical_polygon_area(latsS[:,1],np.rad2deg(lons_b),r_a/1000)
                    areaS[inotE] = areaSeq_b - areaSpol_b          
                    
                    # calculate perturbed areas
                    areaNpol_b = aacgm_functions.spherical_polygon_area(np.flip(latsNplus[:,0]),np.flip(np.rad2deg(lons_b)),r_a/1000)
                    areaNeq_b  = aacgm_functions.spherical_polygon_area(np.flip(latsNplus[:,1]),np.flip(np.rad2deg(lons_b)),r_a/1000)
                    areaNplus[inotE] = areaNeq_b - areaNpol_b
                    
                    areaSpol_b = aacgm_functions.spherical_polygon_area(latsSplus[:,0],np.rad2deg(lons_b),r_a/1000)
                    areaSeq_b  = aacgm_functions.spherical_polygon_area(latsSplus[:,1],np.rad2deg(lons_b),r_a/1000)
                    areaSplus[inotE] = areaSeq_b - areaSpol_b  
                    
                    #############################################
                    # area size difference green's function calculation
                    
                    G_areaN_dg[inotE] = (areaNplus[inotE] - areaN[inotE]) / dg[inotE]
                    G_areaS_dg[inotE] = (areaSplus[inotE] - areaS[inotE]) / dg[inotE]
                                            
                    #############################################
                    # surface coverage difference green's function calculation
                    
                    # absolute changes                
                    G_N_plus_dg[inotE]  =   areaNdiff_plus[inotE]  / dg[inotE] 
                    G_N_minus_dg[inotE] =   areaNdiff_minus[inotE]  / dg[inotE] 
                
                    G_S_plus_dg[inotE]  =   areaSdiff_plus[inotE]  / dg[inotE] 
                    G_S_minus_dg[inotE] =   areaSdiff_minus[inotE]  / dg[inotE]  

                    toc = time.time()
                    elapsed_time = toc-tic
                    print('elapsed time for current Gauss coefficient: '+str(elapsed_time)+' seconds')
                # end of rel_increase

                    
                
    ##################
    # save some output
    print('')
    print('* saving output of calculation * ')
    if areaN.size != 0: # check that not all calculations have been skipped
        data = {'g/h': coeff_cs,
                'l': coeff_l,
                'm': coeff_m,
                'res_lon[deg]':  np.rad2deg(lon_res_b),
                'dg[nT]': dg,
                'g[nT]': g,
                'areaN[km^2]':areaN,
                'areaS[km^2]':areaS,
                'G_areaN_dg[km^2/nT]':G_areaN_dg,
                'G_areaS_dg[km^2/nT]':G_areaS_dg,
                'G_N_plus_dg[km^2/nT]':G_N_plus_dg,
                'G_N_minus_dg[km^2/nT]': G_N_minus_dg,
                'G_S_plus_dg[km^2/nT]': G_S_plus_dg,
                'G_S_minus_dg[km^2/nT]': G_S_minus_dg,
               }
        
        df = pd.DataFrame (data, columns = ['g/h',
                                            'l',
                                            'm',
                                            'res_lon[deg]',
                                            'dg[nT]',
                                            'g[nT]',
                                            'areaN[km^2]',
                                            'areaS[km^2]',
                                            'G_areaN_dg[km^2/nT]',
                                            'G_areaS_dg[km^2/nT]',
                                            'G_N_plus_dg[km^2/nT]',
                                            'G_N_minus_dg[km^2/nT]',
                                            'G_S_plus_dg[km^2/nT]',
                                            'G_S_minus_dg[km^2/nT]'])
        
        # probably not the smartes way to go about it. Shold try and implement the list approach for each iteration of the for loop
        # https://stackoverflow.com/questions/10715965/create-pandas-dataframe-by-appending-one-row-at-a-time
        if os.path.exists(csv_name):
            df_out = df_in.append(df, ignore_index=True)
        else:
            df_out = df
            
        df_out.reset_index() # because index may be jumbled when appending
        df_out.to_csv(csv_name) 
        
        #df_out = pd.DataFrame(rows_data) 
        #df_out.to_csv(folder_case + 'green_results_'+cs+'_'+str(l)+'_'+str(m)+'.csv') 
                    
                    
                    
                    
                    
                    
    ###########
    # some plots
    print('* plotting * ')

    #reload csv for plotting
    df_in = pd.read_csv(csv_name, index_col=0)
    
    df_plot = df_in[ (df_in["res_lon[deg]"].to_numpy().round(8)==np.rad2deg(lon_res_b).round(8)) ]
    df_plot = df_plot.sort_values(by='dg[nT]')
    
    # reference areas  (they should all be the same)
    area_zone_N = df_plot['areaN[km^2]'][0]
    area_zone_S = df_plot['areaS[km^2]'][0]
    
    
    # area coverage green's functions: matrix with positive values of dg
    G_N_surf_dgplus_matrix = np.zeros((Lmax,2*Lmax+1))
    G_N_surf_dgplus_matrix[:,:] = math.nan

    G_S_surf_dgplus_matrix = np.zeros((Lmax,2*Lmax+1))
    G_S_surf_dgplus_matrix[:,:] = math.nan
    
    # area coverage, positive-negative
    G_N_surf_signed_dgplus_matrix = np.zeros((Lmax,2*Lmax+1))
    G_N_surf_signed_dgplus_matrix[:,:] = math.nan

    G_S_surf_signed_dgplus_matrix = np.zeros((Lmax,2*Lmax+1))
    G_S_surf_signed_dgplus_matrix[:,:] = math.nan
    
    # surface area change, for positive values of dg
    G_areaN_dgplus_matrix = np.zeros((Lmax,2*Lmax+1))
    G_areaN_dgplus_matrix[:,:] = math.nan

    G_areaS_dgplus_matrix = np.zeros((Lmax,2*Lmax+1))
    G_areaS_dgplus_matrix[:,:] = math.nan
    
    ## Green's function in Br (CMB)
    dAdBc_N = np.zeros(theta.shape)
    dAdBc_S = np.zeros(theta.shape)
    dCdBr_N = np.zeros(theta.shape)
    dCdBr_S = np.zeros(theta.shape)
    
    ## Green's function in Br (surface)
    dAdBr_N = np.zeros(theta.shape)
    dAdBr_S = np.zeros(theta.shape)

    ## to obtain the actual variations:
    # I am using the current SV coeffs
    SVcoeffs_mat = SH_library.lin2matCoeffs(SVcoeffs3)
    # d(surf)/dg * SV = d(surf)/dt
    dAdt_N = np.copy(G_areaN_dgplus_matrix)
    dAdt_S = np.copy(G_areaS_dgplus_matrix)
    # mix up the SV coefficients
    
    # d(coverage)/dg * SV = d(coverage)/dt
    dCdt_N = np.copy(G_N_surf_dgplus_matrix)
    dCdt_S = np.copy(G_S_surf_dgplus_matrix)    
    
    
    
    for l in range(1,Lmax+1,1):
        #df_l = df_plus[df_plus['l']==l] # older version
        df_l = df_plot[df_plot['l']==l]
        
        SVcoeffs_l = SVcoeffs_mat[SVcoeffs_mat[:,0]==l]
        
       
        
        for m in range(0,l+1,1):
            df_m = df_l[df_l['m']==m]
            SVcoeffs_m = SVcoeffs_l[SVcoeffs_l[:,1]==m]
        
            
           
            for cs in ['g','h']:
                if cs =='g':
                    SV = SVcoeffs_m[0,2]
                    
                    
                    df_cs = df_m[df_m['g/h']=='g']
                    # now pick positive or negative:
                    
                    df_P = df_cs[df_cs['dg[nT]']>0]
                    df_M = df_cs[df_cs['dg[nT]']<0]

                    G_areaN_dgplus_matrix[l-1,Lmax+m] = df_cs['G_areaN_dg[km^2/nT]'].mean()
                    G_areaS_dgplus_matrix[l-1,Lmax+m] = df_cs['G_areaS_dg[km^2/nT]'].mean()
#                    else:
#                        if SV<0:
#                            G_areaN_dgplus_matrix[l-1,Lmax+m] = df_p['G_areaN_dg[km^2/nT]'].to_numpy()
#                            G_areaS_dgplus_matrix[l-1,Lmax+m] = df_p['G_areaS_dg[km^2/nT]'].to_numpy()
#                        else:
#                            G_areaN_dgplus_matrix[l-1,Lmax+m] = df_m['G_areaN_dg[km^2/nT]'].to_numpy()
#                            G_areaS_dgplus_matrix[l-1,Lmax+m] = df_m['G_areaS_dg[km^2/nT]'].to_numpy()                        
                                     
                    # will hve to alter these onese too
                    G_N_surf_dgplus_matrix[l-1,Lmax+m] = df_P['G_N_minus_dg[km^2/nT]'].to_numpy() + df_P['G_N_plus_dg[km^2/nT]'].to_numpy()
                    G_S_surf_dgplus_matrix[l-1,Lmax+m] = df_P['G_S_minus_dg[km^2/nT]'].to_numpy() + df_P['G_S_plus_dg[km^2/nT]'].to_numpy()
                    G_N_surf_signed_dgplus_matrix[l-1,Lmax+m] = -df_P['G_N_minus_dg[km^2/nT]'].to_numpy() + df_P['G_N_plus_dg[km^2/nT]'].to_numpy()
                    G_S_surf_signed_dgplus_matrix[l-1,Lmax+m] = -df_P['G_S_minus_dg[km^2/nT]'].to_numpy() + df_P['G_S_plus_dg[km^2/nT]'].to_numpy()
                            

                    
                                        

                    
                    # d(surf)/dg * SV = d(surf)/dt
                    dAdt_N[l-1,Lmax+m] = G_areaN_dgplus_matrix[l-1,Lmax+m]*SV
                    dAdt_S[l-1,Lmax+m] = G_areaS_dgplus_matrix[l-1,Lmax+m]*SV

                    

                    # d(coverage)/dg * SV = d(coverage)/dt
                    dCdt_N[l-1,Lmax+m] = G_N_surf_dgplus_matrix[l-1,Lmax+m]*SV
                    dCdt_S[l-1,Lmax+m] = G_S_surf_dgplus_matrix[l-1,Lmax+m]*SV  
                    
                    # Transform to Green's function in terms of Br at the CMB
                    SH = SH_library.SchmidtSH(l,m,theta,lons,'c')
                    norm = 4*np.pi/(2*l+1)
                    dgdBc_map =  (1/(l+1)) *  (r_cmb/r_a)**(l+2) * SH / norm
                    dAdBc_N = dAdBc_N + G_areaN_dgplus_matrix[l-1,Lmax+m]*dgdBc_map
                    dAdBc_S = dAdBc_S + G_areaS_dgplus_matrix[l-1,Lmax+m]*dgdBc_map
                    dCdBr_N = dCdBr_N + G_N_surf_dgplus_matrix[l-1,Lmax+m]*dgdBc_map
                    dCdBr_S = dCdBr_S + G_S_surf_dgplus_matrix[l-1,Lmax+m]*dgdBc_map                    

                    dgdBr_map =  (1/(l+1)) *  (r_a/r_a)**(l+2) * SH / norm
                    dAdBr_N = dAdBr_N + G_areaN_dgplus_matrix[l-1,Lmax+m]*dgdBr_map
                    dAdBr_S = dAdBr_S + G_areaS_dgplus_matrix[l-1,Lmax+m]*dgdBr_map
                    
                   
                        
                else:
                    if m==0:
                        continue
                    else:
                        #df_cs = df_m[df_m['g/h']=='h']
                        SV = SVcoeffs_m[0,3]
                       
                        #G_N_surf_dgplus_matrix[l-1,Lmax-m] = df_cs['G_N_minus_dg[km^2/nT]'].to_numpy() + df_cs['G_N_plus_dg[km^2/nT]'].to_numpy()
                        #G_S_surf_dgplus_matrix[l-1,Lmax-m] = df_cs['G_S_minus_dg[km^2/nT]'].to_numpy() + df_cs['G_S_plus_dg[km^2/nT]'].to_numpy()
                        #G_N_surf_signed_dgplus_matrix[l-1,Lmax-m] = -df_cs['G_N_minus_dg[km^2/nT]'].to_numpy() + df_cs['G_N_plus_dg[km^2/nT]'].to_numpy()
                        #G_S_surf_signed_dgplus_matrix[l-1,Lmax-m] = -df_cs['G_S_minus_dg[km^2/nT]'].to_numpy() + df_cs['G_S_plus_dg[km^2/nT]'].to_numpy()
                        #G_areaN_dgplus_matrix[l-1,Lmax-m] = df_cs['G_areaN_dg[km^2/nT]'].to_numpy()
                        #G_areaS_dgplus_matrix[l-1,Lmax-m] = df_cs['G_areaS_dg[km^2/nT]'].to_numpy()
                        
                        df_cs = df_m[df_m['g/h']=='h']
                        # now pick positive or negative:
                        df_P = df_cs[df_cs['dg[nT]']>0]
                        df_M = df_cs[df_cs['dg[nT]']<0]

                        #if df_p['G_areaN_dg[km^2/nT]']*df_m['G_areaS_dg[km^2/nT]']>0: # same sign
                        G_areaN_dgplus_matrix[l-1,Lmax-m] = df_cs['G_areaN_dg[km^2/nT]'].mean()
                        G_areaS_dgplus_matrix[l-1,Lmax-m] = df_cs['G_areaS_dg[km^2/nT]'].mean()
                        #else:
                        #    if SV<0:
                        #        G_areaN_dgplus_matrix[l-1,Lmax-m] = df_p['G_areaN_dg[km^2/nT]'].to_numpy()
                        #        G_areaS_dgplus_matrix[l-1,Lmax-m] = df_p['G_areaS_dg[km^2/nT]'].to_numpy()
                        #    else:
                        #        G_areaN_dgplus_matrix[l-1,Lmax-m] = df_m['G_areaN_dg[km^2/nT]'].to_numpy()
                        #        G_areaS_dgplus_matrix[l-1,Lmax-m] = df_m['G_areaS_dg[km^2/nT]'].to_numpy()                        
                                         
                        # will hve to alter these onese too
                        G_N_surf_dgplus_matrix[l-1,Lmax-m] = df_P['G_N_minus_dg[km^2/nT]'].to_numpy() + df_P['G_N_plus_dg[km^2/nT]'].to_numpy()
                        G_S_surf_dgplus_matrix[l-1,Lmax-m] = df_P['G_S_minus_dg[km^2/nT]'].to_numpy() + df_P['G_S_plus_dg[km^2/nT]'].to_numpy()
                        G_N_surf_signed_dgplus_matrix[l-1,Lmax-m] = -df_P['G_N_minus_dg[km^2/nT]'].to_numpy() + df_P['G_N_plus_dg[km^2/nT]'].to_numpy()
                        G_S_surf_signed_dgplus_matrix[l-1,Lmax-m] = -df_P['G_S_minus_dg[km^2/nT]'].to_numpy() + df_P['G_S_plus_dg[km^2/nT]'].to_numpy()
                          
                        # d(surf)/dg * SV = d(surf)/dt
                        dAdt_N[l-1,Lmax-m] = G_areaN_dgplus_matrix[l-1,Lmax-m]*SV
                        dAdt_S[l-1,Lmax-m] = G_areaS_dgplus_matrix[l-1,Lmax-m]*SV
                        
                        
                        
                        # d(coverage)/dg * SV = d(coverage)/dt
                        dCdt_N[l-1,Lmax-m] = G_N_surf_dgplus_matrix[l-1,Lmax-m]*SV
                        dCdt_S[l-1,Lmax-m] = G_S_surf_dgplus_matrix[l-1,Lmax-m]*SV  
                        
                        # Transform to Green's function in terms of Br at the CMB
                        SH = SH_library.SchmidtSH(l,m,theta,lons,'s')
                        norm = 4*np.pi/(2*l+1)
                        dgdBc_map =  (1/(l+1)) *  (r_cmb/r_a)**(l+2) * SH / norm
                        dAdBc_N = dAdBc_N + G_areaN_dgplus_matrix[l-1,Lmax-m]*dgdBc_map
                        dAdBc_S = dAdBc_S + G_areaS_dgplus_matrix[l-1,Lmax-m]*dgdBc_map
                        dCdBr_N = dCdBr_N + G_N_surf_dgplus_matrix[l-1,Lmax-m]*dgdBc_map
                        dCdBr_S = dCdBr_S + G_S_surf_dgplus_matrix[l-1,Lmax-m]*dgdBc_map   

                        dgdBr_map =  (1/(l+1)) *  (r_a/r_a)**(l+2) * SH / norm
                        dAdBr_N = dAdBr_N + G_areaN_dgplus_matrix[l-1,Lmax-m]*dgdBr_map
                        dAdBr_S = dAdBr_S + G_areaS_dgplus_matrix[l-1,Lmax-m]*dgdBr_map
                        
        
    ## surface coverage change                
    fig,axs = plt.subplots(2,1,figsize=(8,8))
    FS = 12
    FT = 15
    FST = 18
    
    # find min and max of colorbar:
    Zmax = np.max([np.nanmax(G_N_surf_dgplus_matrix),np.nanmax(G_S_surf_dgplus_matrix)])    
    Zmin = np.min([np.nanmin(G_N_surf_dgplus_matrix),np.nanmin(G_S_surf_dgplus_matrix)])    
    
    # NORTH
    ax = axs[0]
    Z = G_N_surf_dgplus_matrix 
    cf = ax.matshow(Z
                    #,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale = 1, linthresh = 1)
                    ,norm=colors.LogNorm(vmin=Zmin, vmax=Zmax)
                    ,cmap = 'plasma'
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
    
    
    ax.set_title(r'Northern Hemisphere, $\mathcal{C}_N(\delta\beta_l^m)$', fontsize=FT, y=1.13)    

    # SOUTH
    ax = axs[1]
    Z = G_S_surf_dgplus_matrix 
    cf = ax.matshow(Z
                    #,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale = 1, linthresh = 1)
                    ,norm=colors.LogNorm(vmin=Zmin, vmax=Zmax)
                    ,cmap = 'plasma'
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
    
    
    ax.set_title(r'Southern Hemisphere, $\mathcal{C}_S(\delta\beta_l^m)$', fontsize=FT, y=1.15)    

    #clb = fig.colorbar(cf, ax=axs.ravel().tolist(), shrink=0.95, pad=0.1)
    cb_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
    clb = fig.colorbar(cf, cax=cb_ax)

    clb.ax.tick_params(labelsize=FS)
    clb.set_label(r'[km$^2$/nT]',fontsize=FS)
    
    plt.suptitle('Sensitivity of auroral zones surface coverage change', fontsize=FST)
    plt.tight_layout()
                        
    plt.savefig(folder_quantity + 'G_surf_unsigned_dgplus_matrix_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show(block=False)     
    




    ## surface coverage (signed) change                
    fig,axs = plt.subplots(2,1,figsize=(8,8))
    FS = 12
    FT = 15
    FST = 18
    
    # find min and max of colorbar:
    Zmax = np.max([np.nanmax(G_N_surf_signed_dgplus_matrix),np.nanmax(G_S_surf_signed_dgplus_matrix)])    
    Zmin = np.min([np.nanmin(G_N_surf_signed_dgplus_matrix),np.nanmin(G_S_surf_signed_dgplus_matrix)])    
    
    # NORTH
    ax = axs[0]
    Z = G_N_surf_signed_dgplus_matrix 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale = 1, linthresh = 10)
                    #,norm=colors.LogNorm(vmin=Zmin, vmax=Zmax)
                    ,cmap = 'RdYlBu_r'
                    )

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
    
    
    ax.set_title(r'Northern Hemisphere, $\mathcal{C}_N(\delta\beta_l^m)$', fontsize=FT, y=1.13)    

    # SOUTH
    ax = axs[1]
    Z = G_S_surf_signed_dgplus_matrix 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale =1, linthresh = 10)
                    #,norm=colors.LogNorm(vmin=Zmin, vmax=Zmax)
                    ,cmap = 'RdYlBu_r'
                    )


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
    
    
    ax.set_title(r'Southern Hemisphere, $\mathcal{C}_S(\delta\beta_l^m)$', fontsize=FT, y=1.15)    

    #clb = fig.colorbar(cf, ax=axs.ravel().tolist(), shrink=0.95, pad=0.1)
    cb_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
    clb = fig.colorbar(cf, cax=cb_ax)

    clb.ax.tick_params(labelsize=FS)
    clb.set_label(r'[km$^2$/nT]',fontsize=FS)
    
    plt.suptitle('Sensitivity of auroral zones surface area change (signed)', fontsize=FST)
    plt.tight_layout()
                        
    plt.savefig(folder_quantity + 'G_surf_signed_dgplus_matrix_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show(block=False)     
    
    
    
    

    

    ## surface area change                
    fig,axs = plt.subplots(2,1,figsize=(8,8))
    FS = 12
    FT = 15
    FST = 18
    
    # find min and max of colorbar:
    Zmax = np.max([np.nanmax(G_areaN_dgplus_matrix),np.nanmax(G_areaS_dgplus_matrix)])    
    Zmin = np.min([np.nanmin(G_areaN_dgplus_matrix),np.nanmin(G_areaS_dgplus_matrix)])    
    
    # NORTH
    ax = axs[0]
    Z = G_areaN_dgplus_matrix 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale = 1, linthresh = 10)
                    #,norm=colors.LogNorm(vmin=Zmin, vmax=Zmax)
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
    
    
    #ax.set_title(r'Northern Hemisphere, $\mathcal{A}_N(\delta\beta_l^m)$', fontsize=FT, y=1.13)    
    ax.set_title(r'Northern Hemisphere, $\frac{\partial {A}_N}{\partial\beta_l^m}$', fontsize=FT, y=1.15)    

    # SOUTH
    ax = axs[1]
    Z = G_areaS_dgplus_matrix 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale =1, linthresh = 10)
                    #,norm=colors.LogNorm(vmin=Zmin, vmax=Zmax)
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
    
    
    #ax.set_title(r'Southern Hemisphere, $\mathcal{A}_S(\delta\beta_l^m)$', fontsize=FT, y=1.15)    
    ax.set_title(r'Southern Hemisphere, $\frac{\partial {A}_S}{\partial\beta_l^m}$', fontsize=FT, y=1.15)    

    #clb = fig.colorbar(cf, ax=axs.ravel().tolist(), shrink=0.95, pad=0.1)
    cb_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
    clb = fig.colorbar(cf, cax=cb_ax)

    clb.ax.tick_params(labelsize=FS)
    clb.set_label(r'[km$^2$/nT]',fontsize=FS)
    
    plt.suptitle('Sensitivity of auroral zones surface area change', fontsize=FST)
    plt.tight_layout()
                        
    plt.savefig(folder_quantity + 'G_area_dgplus_matrix_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show(block=False)     
    
    
    LmaxSV = 13
    ## dA/dt               
    fig,axs = plt.subplots(2,1,figsize=(8,8))
    FS = 12
    FT = 15
    FST = 18
    
    # find min and max of colorbar:
    Zmax = np.max([np.nanmax(dAdt_N),np.nanmax(dAdt_S)])    
    Zmin = np.min([np.nanmin(dAdt_N),np.nanmin(dAdt_S)])    
    
    # NORTH
    ax = axs[0]
    Z = dAdt_N[:LmaxSV,Lmax-LmaxSV:2*Lmax+1-(Lmax-LmaxSV)] 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale = 1, linthresh = 10)
                    ,cmap = 'RdYlBu_r'
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
    
    
    ax.set_title(r'Northern Hemisphere, $\frac{\delta A_N(\delta\beta_l^m)}{\delta t}$', fontsize=FT, y=1.13)    

    # SOUTH
    ax = axs[1]
    Z = dAdt_S[:LmaxSV,Lmax-LmaxSV:2*Lmax+1-(Lmax-LmaxSV)] 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale =1, linthresh = 10)
                    ,cmap = 'RdYlBu_r'
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
    
    
    ax.set_title(r'Southern Hemisphere, $\frac{\delta A_S(\delta\beta_l^m)}{\delta t}$', fontsize=FT, y=1.15)    

    #clb = fig.colorbar(cf, ax=axs.ravel().tolist(), shrink=0.95, pad=0.1)
    cb_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
    clb = fig.colorbar(cf, cax=cb_ax)

    clb.ax.tick_params(labelsize=FS)
    clb.set_label(r'[km$^2$/year]',fontsize=FS)
    
    plt.suptitle('Auroral zones surface area yearly change', fontsize=FST)
    plt.tight_layout()
                        
    plt.savefig(folder_quantity + 'dAdt_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show(block=False)  
    
    
    ## dC/dt               
    fig,axs = plt.subplots(2,1,figsize=(8,8))
    FS = 12
    FT = 15
    FST = 18
    
    # find min and max of colorbar:
    Zmax = np.max([np.nanmax(abs(dCdt_N)),np.nanmax(abs(dCdt_S))])    
    Zmin = np.min([np.nanmin(abs(dCdt_N)),np.nanmin(abs(dCdt_S))])    
    
    # NORTH
    ax = axs[0]
    Z = abs(dCdt_N[:LmaxSV,Lmax-LmaxSV:2*Lmax+1-(Lmax-LmaxSV)]) 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale =1, linthresh = 10)
                    ,cmap = 'plasma'
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
    
    
    ax.set_title(r'Northern Hemisphere, $\frac{\delta C_N(\delta\beta_l^m)}{\delta t}$', fontsize=FT, y=1.13)    

    # SOUTH
    ax = axs[1]
    Z = abs(dCdt_S[:LmaxSV,Lmax-LmaxSV:2*Lmax+1-(Lmax-LmaxSV)]) 
    cf = ax.matshow(Z
                    ,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale =1, linthresh = 10)
                    ,cmap = 'plasma'
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
    
    
    ax.set_title(r'Southern Hemisphere, $\frac{\delta C_S(\delta\beta_l^m)}{\delta t}$', fontsize=FT, y=1.15)    

    #clb = fig.colorbar(cf, ax=axs.ravel().tolist(), shrink=0.95, pad=0.1)
    cb_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
    clb = fig.colorbar(cf, cax=cb_ax)

    clb.ax.tick_params(labelsize=FS)
    clb.set_label(r'[km$^2$/year]',fontsize=FS)
    
    plt.suptitle('Auroral zones surface coverage yearly change', fontsize=FST)
    plt.tight_layout()
                        
    plt.savefig(folder_quantity + 'dCdt_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show(block=False)  
    
    ## same plot as the last two, but as a function of l only
    
    # dAdt
    
    FS = 20
    
    fig,ax = plt.subplots(figsize=(8,5))
    ax.set_xlabel(r'$l$',fontsize=FS)
    ax.set_ylabel('km$^2$ / year ', fontsize=FS)
    #ax.set_ylim([77,81])
    #ax.set_xlim([1900,2100])
    ax.tick_params('x',labelsize=FS)
    ax.axhline(y=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    
    ax.plot(np.arange(1,LmaxSV+1,1),np.nansum(dAdt_N[:LmaxSV,Lmax-LmaxSV:2*Lmax+1-(Lmax-LmaxSV)],axis=1),'o-',color='tab:blue',label=r'North')
    
    ax.plot(np.arange(1,LmaxSV+1,1),np.nansum(dAdt_S[:LmaxSV,Lmax-LmaxSV:2*Lmax+1-(Lmax-LmaxSV)],axis=1),'o-',color='tab:orange',label=r'South')
    
    ax.legend(fontsize=15,loc='upper right')
    
    ax.tick_params('y',labelsize=FS)
    #ax.yaxis.get_offset_text().set_fontsize(12)
    #ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    #ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    plt.title(r'Auroral zone surface area yearly change',fontsize=FS)
    plt.savefig(folder_quantity + 'dAdt_l_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()

    
    # dCdt
    fig,ax = plt.subplots(figsize=(8,5))
    ax.set_xlabel(r'$l$',fontsize=FS)
    ax.set_ylabel('km$^2$ / year ', fontsize=FS)
    #ax.set_ylim([77,81])
    #ax.set_xlim([1900,2100])
    ax.tick_params('x',labelsize=FS)
    ax.axhline(y=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    
    ax.plot(np.arange(1,LmaxSV+1,1),np.nansum(abs(dCdt_N[:LmaxSV,Lmax-LmaxSV:2*Lmax+1-(Lmax-LmaxSV)]),axis=1),'o-',color='tab:blue',label=r'North')
    
    ax.plot(np.arange(1,LmaxSV+1,1),np.nansum(abs(dCdt_S[:LmaxSV,Lmax-LmaxSV:2*Lmax+1-(Lmax-LmaxSV)]),axis=1),'o-',color='tab:orange',label=r'South')
    
    ax.legend(fontsize=15,loc='upper right')
    
    ax.tick_params('y',labelsize=FS)
    ax.yaxis.get_offset_text().set_fontsize(12)
    #ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    #ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    plt.title(r'Auroral zone surface coverage yearly change',fontsize=FS)
    plt.savefig(folder_quantity + 'dCdt_l_'+str(target_res)+'_target_res_zero_'+str(target_res_zero)+'_exponential_increase_res_lon_'+str(np.rad2deg(lon_res_b))+'deg_bisection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()

    
    #####################################
    # Sensitivity to dipole tilt and dipole magnitude
    
    g10 = MFcoeffs[idx_time,1]
    g11 = MFcoeffs[idx_time,2]
    h11 = MFcoeffs[idx_time,3]
    m1 = np.sqrt(g10**2 + g11**2 + h11**2)
    theta_d = np.arccos(g10/m1)
    phi_d = np.arctan(h11/g11)
    
    g10dot = SVcoeffs3[0]
    g11dot = SVcoeffs3[1]
    h11dot = SVcoeffs3[2]
    
    m1_dot = np.sqrt((g10+g10dot)**2 + (g11+g11dot)**2 + (h11+h11dot)**2) - m1
    theta_d_dot = np.arccos((g10+g10dot)/(m1+m1_dot)) - theta_d
    phi_d_dot = np.arctan((h11+h11dot)/(g11+g11dot)) - phi_d


    # north
    dCdg10_N = G_N_surf_dgplus_matrix[0,13]
    dCdg11_N = G_N_surf_dgplus_matrix[0,14]
    dCdh11_N = G_N_surf_dgplus_matrix[0,12]
    
    dAdg10_N = G_areaN_dgplus_matrix[0,13]
    dAdg11_N = G_areaN_dgplus_matrix[0,14]
    dAdh11_N = G_areaN_dgplus_matrix[0,12]  
    
    # usefuk?
    dCdm_N = m1 * ( dCdg10_N/g10 + dCdg11_N/g11 + dCdh11_N/h11 )
    dCdthetad_N = m1**2 * ( -dCdg10_N + dCdg11_N*(g11**2+h11**2)/(g11*g10) + dCdh11_N*(g11**2+h11**2)/(h11*g10) ) / np.sqrt(g11**2+h11**2)  
    
    dAdm_N = m1 * ( dAdg10_N/g10 + dAdg11_N/g11 + dAdh11_N/h11 )
    dAdthetad_N = m1**2 * ( -dAdg10_N + dAdg11_N*(g11**2+h11**2)/(g11*g10) + dAdh11_N*(g11**2+h11**2)/(h11*g10) ) / np.sqrt(g11**2+h11**2)
    dAdphid_N = (g11**2+h11**2) * (dAdg11_N/g11 - dAdh11_N/h11 )


    # south
    dCdg10_S = G_S_surf_dgplus_matrix[0,13]
    dCdg11_S = G_S_surf_dgplus_matrix[0,14]
    dCdh11_S = G_S_surf_dgplus_matrix[0,12]
    
    dAdg10_S = G_areaS_dgplus_matrix[0,13]
    dAdg11_S = G_areaS_dgplus_matrix[0,14]
    dAdh11_S = G_areaS_dgplus_matrix[0,12]  
    
    #usefuk?
    dCdm_S = m1 * ( dCdg10_S/g10 + dCdg11_S/g11 + dCdh11_S/h11 )    
    dCdthetad_S = m1**2 * ( -dCdg10_S + dCdg11_S*(g11**2+h11**2)/(g11*g10) + dCdh11_S*(g11**2+h11**2)/(h11*g10) ) / np.sqrt(g11**2+h11**2)
 
    dAdm_S = m1 * ( dAdg10_S/g10 + dAdg11_S/g11 + dAdh11_S/h11 )
    dAdthetad_S = m1**2 * ( -dAdg10_S + dAdg11_S*(g11**2+h11**2)/(g11*g10) + dAdh11_S*(g11**2+h11**2)/(h11*g10) ) / np.sqrt(g11**2+h11**2)
    dAdphid_S = (g11**2+h11**2) * (dAdg11_S/g11 - dAdh11_S/h11 )
      
    
    
    #####################################
    # plot of STD from MPG forecast
    
    
    
    ## Green's functions for A_i as a function of Br(CMB)

    fig = plt.figure(figsize=(10,6))
    
    FS = 20 
    Zs = dAdBc_N
    #Zmax = max(np.nanmax(abs(dAdBc_N)),np.nanmax(abs(dAdBc_S)))
    Zmax = 10.5
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
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsN[:,i], '--',
                 color='r',
                 transform=ccrs.PlateCarree())
        
        
    clb.set_label(r'$\partial A_N/\partial B_c$ [km$^2$/nT]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dAdBc_N.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    
    
    
    fig = plt.figure(figsize=(10,6))
    Zs = dAdBc_S
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
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsS[:,i], '--',
                 color='r',
                 transform=ccrs.PlateCarree())
        
    clb.set_label(r'$\partial A_S/\partial B_c$ [km$^2$/nT]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dAdBc_S.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    

    ## Green's functions for A_i as a function of Br(CMB)

    fig = plt.figure(figsize=(10,6))
    Zs = dAdBr_N
    #Zmax = max(np.nanmax(abs(dAdBr_N)),np.nanmax(abs(dAdBr_S)))
    Zmax = 450
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
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsN[:,i], '--',
                 color='r',
                 transform=ccrs.PlateCarree())
        
        
    clb.set_label(r'$\partial A_N/\partial B_r$ [km$^2$/nT]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dAdBr_N.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    
    fig = plt.figure(figsize=(10,6))
    Zs = dAdBr_S
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
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsS[:,i], '--',
                 color='r',
                 transform=ccrs.PlateCarree())
    
    clb.set_label(r'$\partial A_S/\partial B_r$ [km$^2$/nT]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dAdBr_S.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    

    
    #multiply the above by the actual SV from current model
    
    SVcoeffsmat = SH_library.lin2matCoeffs(SVcoeffs3)
    #SVcoeffsmat = SH_library.lin2matCoeffs(SVcoeffs_igrf_20_25_mat)
    
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
    
    dAdtN_profile =np.zeros(theta[0,:-1].shape)
    dAdtS_profile =np.zeros(theta[0,:-1].shape)

    dAdtN_surf_profile =np.zeros(theta[0,:-1].shape)
    dAdtS_surf_profile =np.zeros(theta[0,:-1].shape)

    
    for ith in range(len(theta[0,:])-1):
        dth = abs(theta[iph,ith+1]-theta[iph,ith])
        for iph in range(len(lons[:,0])-1):
            dph = abs(lons[iph+1,ith]-lons[iph,ith])
            thav = (theta[iph,ith+1]+theta[iph,ith])/2  
            # cmb
            dAdtN_profile[ith] = dAdtN_profile[ith] + dBcdt[iph,ith]*dAdBc_N[iph,ith]*dph*np.sin(theta[iph,ith])
            dAdtS_profile[ith] = dAdtS_profile[ith] + dBcdt[iph,ith]*dAdBc_S[iph,ith]*dph*np.sin(theta[iph,ith])
            #surf
            dAdtN_surf_profile[ith] = dAdtN_surf_profile[ith] + dBrdt[iph,ith]*dAdBr_N[iph,ith]*dph*np.sin(theta[iph,ith])
            dAdtS_surf_profile[ith] = dAdtS_surf_profile[ith] + dBrdt[iph,ith]*dAdBr_S[iph,ith]*dph*np.sin(theta[iph,ith])
    
    dAdtN_tot=0
    dAdtS_tot=0
    dAdtN_surf_tot=0
    dAdtS_surf_tot=0

    # assume uniform phi grid
    dph = abs(lons_lin[1]-lons_lin[0])
    for ilon in range(len(lons_lin)-1):
        for ilat in range(len(theta_lin)):
            # check the actual time diff
            dAdtN_tot      = dAdtN_tot      + dBcdt[ilon,ilat]*dAdBc_N[ilon,ilat]* dph * weights[ilat]
            dAdtN_surf_tot = dAdtN_surf_tot + dBrdt[ilon,ilat]*dAdBr_N[ilon,ilat]* dph * weights[ilat]
            dAdtS_tot      = dAdtS_tot      + dBcdt[ilon,ilat]*dAdBc_S[ilon,ilat]* dph * weights[ilat]
            dAdtS_surf_tot = dAdtS_surf_tot + dBrdt[ilon,ilat]*dAdBr_S[ilon,ilat]* dph * weights[ilat]
            
    
    
    
    fig = plt.figure(figsize=(10,6))
    expmax=5
    
    Zs = dBcdt*dAdBc_N/10**expmax
    Zmax = 1.5
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
    
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsN[:,i], '--',
                 color='r',
                 transform=ccrs.PlateCarree())
        
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r'$\partial A_N/\partial t  \times 10^{5}$ (CMB) [km$^2$/yr]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dAdt_N.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()   
    
    
    
    fig = plt.figure(figsize=(10,6))
    Zs = dBcdt*dAdBc_S/10**expmax
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
    
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsS[:,i], '--',
                 color='r',
                 transform=ccrs.PlateCarree())
        
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r'$\partial A_S/\partial t  \times 10^{5}$ (CMB) [km$^2$/yr]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dAdt_S.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()   
    
    
    # integrated in longitude
    fig,ax = plt.subplots(figsize=(5,8))
    ax.set_ylabel(r'lat [deg]',fontsize=FS)
    ax.set_xlabel(r'km$^2$ / year $\times 10^{5}$ ', fontsize=FS)
    #ax.set_ylim([77,81])
    #ax.set_xlim([1900,2100])
    ax.tick_params('x',labelsize=FS)
    
    ax.plot(dAdtN_profile/10**expmax,90-180*theta[0,:-1]/np.pi,color='tab:blue',label=r'North')
    
    ax.plot(dAdtS_profile/10**expmax,90-180*theta[0,:-1]/np.pi,color='tab:orange',label=r'South')
    
    ax.axvline(x=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')

    ax.legend(fontsize=15,loc='upper right')
    
    ax.tick_params('y',labelsize=FS)
    ax.yaxis.get_offset_text().set_fontsize(12)
    #ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    #ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    plt.title(r'Auroral zone surface area yearly change (CMB)',fontsize=FS)
    plt.savefig(folder_quantity + 'dAdt_lats.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()



    #Earth surface

    fig = plt.figure(figsize=(10,6))
    expmax=4
    Zs = dBrdt*dAdBr_N/10**expmax
    Zmax=2.8
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
    
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsN[:,i], '--',
                 color='r',
                 transform=ccrs.PlateCarree())
        
        
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r"$\partial A_N/\partial t  \times 10^{4}$ (Earth's surface) [km$^2$/yr]", fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dAdt_N_surf.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()   
    
    
    
    fig = plt.figure(figsize=(10,6))
    Zs = dBrdt*dAdBr_S/10**expmax
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
    
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsS[:,i], '--',
                 color='r',
                 transform=ccrs.PlateCarree())
        
        
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r"$\partial A_S/\partial t  \times 10^{4}$  (Earth's surface) [km$^2$/yr]", fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dAdt_S_surf.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()   
    
 
    # integrated in longitude
    fig,ax = plt.subplots(figsize=(5,8))
    ax.set_ylabel(r'lat  [deg]',fontsize=FS)
    ax.set_xlabel(r'km$^2$ / year $\times 10^{4}$', fontsize=FS)
    #ax.set_ylim([77,81])
    #ax.set_xlim([1900,2100])
    ax.tick_params('x',labelsize=FS)
    
    ax.plot(dAdtN_surf_profile/10**expmax,90-180*theta[0,:-1]/np.pi,color='tab:blue',label=r'North')
    
    ax.plot(dAdtS_surf_profile/10**expmax,90-180*theta[0,:-1]/np.pi,color='tab:orange',label=r'South')
    
    ax.axvline(x=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')

    ax.legend(fontsize=15,loc='upper right')
    
    ax.tick_params('y',labelsize=FS)
    ax.yaxis.get_offset_text().set_fontsize(12)
    #ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    #ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
    plt.title(r"Auroral zone surface area yearly change (Earth's surface)",fontsize=FS)
    plt.savefig(folder_quantity + 'dAdt_surf_lats.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()
    
    ## Green's functions for C_i as a function of Br(CMB)

    fig = plt.figure(figsize=(10,6))
    Zs = dCdBr_N
    #Zmax = max(np.nanmax(abs(dCdBr_N)),np.nanmax(abs(dCdBr_S)))
    Zmax = 165
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
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r'$\partial C_N/\partial Br$ [km$^2$/nT]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dCdBr_N.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    
    fig = plt.figure(figsize=(10,6))
    Zs = dCdBr_S
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
    clb = plt.colorbar(cf
                       #,ticks=[-90,-45,0,45,90]
                       ,fraction=0.05
                       ,orientation='horizontal'
                       ,aspect = 40
                       ,pad = 0.05
                       #,format='%.1e'
                       )
    clb.set_label(r'$\partial C_S/\partial Br$ [km$^2$/nT]', fontsize = FS)
    #clb.ax.set_xticklabels(['-90','-45','0','45','90'])
    clb.ax.tick_params(labelsize=FS)
    #plt.title('AACGM latitudes',fontsize=20)
    plt.savefig(folder_quantity + '/dCdBr_S.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()    
    
    
    
    
    ##########
    # cleanup    
    os.system('rm '+folder_out+'CGMlats_'+str(t_ref)+'_'+filename3)
    os.system('rm '+folder_out+'CGMlons_'+str(t_ref)+'_'+filename3)
    os.system('rm '+folder_out+'CGMmlt_'+str(t_ref)+'_'+filename3)
    
    os.system('rm '+folder_quantity+'CGMlats_'+str(t_ref)+'_'+'mymodel.txt')
    os.system('rm '+folder_quantity+'CGMlons_'+str(t_ref)+'_'+'mymodel.txt')
    os.system('rm '+folder_quantity+'CGMmlt_'+str(t_ref)+'_'+'mymodel.txt')

    os.system('rm '+folder_quantity+'mymodel.txt')

    
else:
   print("Something went wrong")