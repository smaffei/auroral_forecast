#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 16:40:01 2021

@author: stefanomaffei



"""

import aacgmv2
import aacgm_functions # my module to implement aacgmv2
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import math
import os
from multiprocessing import Process, Pool, Queue
import pandas as pd
import time
import cartopy.feature as cfeature


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



# define linestiles:
igrf_LS = '--'
mpg_LS  = '-'
ipgp_LS = ':'

# parameters for ensemble plots
Nens_plot = 64 # how many members to plot
alpha_ens = 0.1
lw_ens = 0.5


cities_file = "../datasets/simplemaps_worldcities_basicv1.73/worldcities.csv"
cities = pd.read_csv(cities_file)

# for auroral zones
city_list_aurorae_N = ["Tromso","Reykjavik","Yellowknife","Nuuk","Fairbanks","Salekhard","Yakutsk","Kodiak","Edmonton","Stockholm","Toronto"] # "Rovaniemi",
off_lon_aurorae_N = [    1,         3,             -2,       2,      21,             2,      0,         15,      2,          0,          -3]
off_lat_aurorae_N = [    -4,         -4,            1,       2,      -15,             -1,      -1,        -10,      -12,          -4,          -12]
city_list_aurorae_S = ["Hobart","Dunedin","Perth"]
off_lon_aurorae_S = [-2,-2,-2]
off_lat_aurorae_S = [0,0,0]

# for danger zones
city_list_America = ["Wood Buffalo","Yellowknife","Minneapolis","Winnipeg","Edmonton","Kodiak","Fairbanks","Iqaluit","Juneau","San Francisco","Los Angeles","Denver","Charlotte","New York","Quebec City","Seattle","Labrador City"]
city_list_Europe =  [ "Inverness","Oulu","Berlin","Leeds","Oslo","Reykjavik","Moscow","Warsaw","Paris","Copenhagen","London","Kyiv","Trondheim","Stockholm","Minsk","Vienna","Omsk","Perm","Nur-Sultan","Ulaanbaatar"]
city_list_Oceania =  ["Camberra","Perth","Albany","Melbourne","Camberra","Sydney","Hobart","Dunedin","Wellington","Auckland"]


FSL = 10
    
    
######################################
# define some functions
######################################
def CGMlat_bisection(folder_quantity,t1t,t_ref,folder3,filename3):
    if not os.path.exists(folder_quantity):
        os.makedirs(folder_quantity) 
    # prepare bisecting longitude. If the reference zones are already calculated, load them        
    filename_latsN = folder_quantity+"latsN_"+str(t_ref)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3
    filename_latsS = folder_quantity+"latsS_"+str(t_ref)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3

    if os.path.exists(filename_latsN) and os.path.exists(filename_latsS):
        lons_b = np.loadtxt(filename_latsN)[:,0]*np.pi/180
    else:
        lons_b = np.linspace(0,2*np.pi,num=int(2*np.pi/lon_res_b))
    
    # lats and lons for plotting the coordinates
    lats_lin = np.linspace(-np.pi/2,np.pi/2,num=180)
    lons_lin = np.linspace(0,2*np.pi,num=360)
    
    lats, lons = np.meshgrid(lats_lin,lons_lin)
    theta = -lats+np.pi/2
    
    #############################
    # reference model auroral zones
    if os.path.exists(filename_latsN) and os.path.exists(filename_latsS):
        print('loading auroral zone boundaries')
        latsN = np.loadtxt(filename_latsN)[:,1:3]
        latsS = np.loadtxt(filename_latsS)[:,1:3]
    else:
        print('calculating auroral zone boundaries')
        
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
    return lons_b, latsN, latsS
        


# same function, but for danger zone (it's all very messy)
# nominal latitudinal bounds for danger zones:
latNdanger_pol_target = 60
latNdanger_eq_target  = 50
latSdanger_pol_target = -60
latSdanger_eq_target  = -50
# first guess for the search for the bisection algorithm
latNdanger_pol_i = 40
latNdanger_eq_i = 79
latSdanger_pol_i = -79
latSdanger_eq_i  = -40

def CGMdanger_lat_bisection(folder_quantity,t1t,t_ref,folder3,filename3):
    if not os.path.exists(folder_quantity):
        os.makedirs(folder_quantity) 
    # prepare bisecting longitude. If the reference zones are already calculated, load them        
    filename_latsN = folder_quantity+"latsNdanger_"+str(t_ref)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3
    filename_latsS = folder_quantity+"latsSdanger_"+str(t_ref)+"_target_res_"+str(target_res)+"_target_res_zero_"+str(target_res_zero)+"_res_lon_"+str(np.rad2deg(lon_res_b))+"_"+filename3

    if os.path.exists(filename_latsN) and os.path.exists(filename_latsS):
        lons_b = np.loadtxt(filename_latsN)[:,0]*np.pi/180
    else:
        lons_b = np.linspace(0,2*np.pi,num=int(2*np.pi/lon_res_b))
    
    # lats and lons for plotting the coordinates
    lats_lin = np.linspace(-np.pi/2,np.pi/2,num=180)
    lons_lin = np.linspace(0,2*np.pi,num=360)
    
    lats, lons = np.meshgrid(lats_lin,lons_lin)
    theta = -lats+np.pi/2
    
    #############################
    # reference model auroral zones
    if os.path.exists(filename_latsN) and os.path.exists(filename_latsS):
        print('loading danger zone boundaries')
        latsN = np.loadtxt(filename_latsN)[:,1:3]
        latsS = np.loadtxt(filename_latsS)[:,1:3]
    else:
        print('calculating danger zone boundaries')
        
        # initialise lats arrays: n x 2 where n is range(lons), and the other dimension contains the polward and equatorward boundaries
        latsN = np.zeros((lons_b.shape[0],2))
        latsS = np.zeros((lons_b.shape[0],2))
        
        # northern zone
        # polar edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t1t, 1, 1, np.rad2deg(lons_b), latNdanger_pol_target, target_res, latNdanger_pol_i, latNdanger_eq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        # equatorial edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t1t, 1, 1, np.rad2deg(lons_b), latNdanger_eq_target, target_res, latNdanger_pol_i, latNdanger_eq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        
        # southern zone
        # polar edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t1t, 1, 1, np.rad2deg(lons_b), latSdanger_pol_target, target_res, latSdanger_pol_i, latSdanger_eq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
        # equatorial edge
        my_process3 = Process(target=aacgm_functions.lat_bisection, args=(folder3,filename3, t1t, 1, 1, np.rad2deg(lons_b), latSdanger_eq_target, target_res, latSdanger_pol_i, latSdanger_eq_i, max_iterations,folder_out))
        my_process3.start()
        my_process3.join() 
                     
        # read the calculated values
        latsN[:,0] = np.loadtxt(folder_out+'bisection_'+str(latNdanger_pol_target)+'_lats_'+filename3)
        latsN[:,1] = np.loadtxt(folder_out+'bisection_'+str(latNdanger_eq_target)+'_lats_'+filename3)
    
        np.savetxt(filename_latsN,np.hstack([lons_b[:,None]*180/np.pi, latsN]),header='lons lats(polar edge) lats(equatorward edge)')
        os.system('rm '+folder_out+'bisection_'+str(latNdanger_pol_target)+'_lats_'+filename3)
        os.system('rm '+folder_out+'bisection_'+str(latNdanger_eq_target)+'_lats_'+filename3)
        
        
        latsS[:,0] = np.loadtxt(folder_out+'bisection_'+str(latSdanger_pol_target)+'_lats_'+filename3)
        latsS[:,1] = np.loadtxt(folder_out+'bisection_'+str(latSdanger_eq_target)+'_lats_'+filename3)
        
        np.savetxt(filename_latsS,np.hstack([lons_b[:,None]*180/np.pi, latsS]),header='lons lats(polar edge) lats(equatorward edge)')
        os.system('rm '+folder_out+'bisection_'+str(latSdanger_pol_target)+'_lats_'+filename3)
        os.system('rm '+folder_out+'bisection_'+str(latSdanger_eq_target)+'_lats_'+filename3)
    return lons_b, latsN, latsS

###############
# main function    
if __name__ == "__main__":   
    
    plt.close('all')

    # magmodel file
    # provided in the aacgmv2 package. 
    # relocated in ../models/ for convenience  
    folder_magmodel   = '../models/'
    filename_magmodel = 'magmodel_1590-2020.txt'

    # igrf-13 forecast
    igrf_folder = '../models/'
    igrf_file   = 'magmodel_2020-2140.txt'
    igrf_name    = 'igrf'
    igrf_t_offset = 2020-1590
    igrf_t_in = 2020
    igrf_t_fin = 2100    

    # Aubert's 2015 forecast
    ipgp_folder = '../models/'    
    ipgp_file   = 'ipgpMF4aacgmv2_from1590.txt' 
    ipgp_name    = 'ipgp15'
    ipgp_t_offset = 2015-1590
    ipgp_t_in = 2015
    ipgp_t_fin = 2100    
    
    # mpg forecast 
    mpg_folder = '../models/'    
    mpg_file   = 'mpgMF4aacgmv2_from1590_L13.txt' 
    mpg_name    = 'mpg20'
    mpg_t_offset = 2015-1590
    mpg_t_in = 2015
    mpg_t_fin = 2100

    
    
    folder_base = './'
    folder_out = folder_base+'aacgmv2_coords/'
    folder_models = folder_base+'models/'
    
    folder_figures = folder_base+'figures/'
    
    os.makedirs(folder_figures)
    
    # magmodel, 2020
    folder_coords = folder_base+'coords'+str(2020)+'_bisection/'
    folder_quantity = folder_coords + 'auroral_zones/' 
    lons_b, latsN2020, latsS2020 = CGMlat_bisection(folder_quantity,2020,2020,folder_magmodel,filename_magmodel)

    folder_coords = folder_base+'coords'+str(2020)+'_bisection/'
    folder_quantity = folder_coords + 'danger_zones/' 
    lons_b, latsNdanger_2020, latsSdanger_2020 = CGMdanger_lat_bisection(folder_quantity,2020,2020,folder_magmodel,filename_magmodel)

    # igrf forecast
    folder_coords = folder_base+'coords'+str(2070)+'_bisection/'
    folder_quantity = folder_coords + 'auroral_zones/' 
    lons_b, latsN_igrf, latsS_igrf = CGMlat_bisection(folder_quantity,2070 -igrf_t_offset,2070 ,igrf_folder,igrf_file)    

    folder_coords = folder_base+'coords'+str(2070)+'_bisection/'
    folder_quantity = folder_coords + 'danger_zones/'     
    _, latsNdanger_igrf, latsSdanger_igrf = CGMdanger_lat_bisection(folder_quantity,2070 -igrf_t_offset,2070 ,igrf_folder,igrf_file)
    
    
    # ipgp
    folder_coords = folder_base+'coords'+str(2070)+'_bisection/'
    folder_quantity = folder_coords + 'auroral_zones/' 
    lons_b, latsN_ipgp, latsS_ipgp = CGMlat_bisection(folder_quantity,2070 -ipgp_t_offset,2070 ,ipgp_folder,ipgp_file)    

    folder_coords = folder_base+'coords'+str(2070)+'_bisection/'
    folder_quantity = folder_coords + 'danger_zones/'     
    _, latsNdanger_ipgp, latsSdanger_ipgp = CGMdanger_lat_bisection(folder_quantity,2070 -ipgp_t_offset,2070 ,ipgp_folder,ipgp_file)
    
    
    # mpg
    folder_coords = folder_base+'coords'+str(2070)+'_bisection/'
    folder_quantity = folder_coords + 'auroral_zones/' 
    lons_b, latsN_mpg, latsS_mpg = CGMlat_bisection(folder_quantity,2070 -mpg_t_offset,2070 ,mpg_folder,mpg_file)

    folder_coords = folder_base+'coords'+str(2070)+'_bisection/'
    folder_quantity = folder_coords + 'danger_zones/'     
    _, latsNdanger_mpg, latsSdanger_mpg = CGMdanger_lat_bisection(folder_quantity,2070 -mpg_t_offset,2070 ,mpg_folder,mpg_file)
    
    
    # mpg, ensemble
    latsN_mpg_ens = np.zeros((latsN_mpg.shape[0],latsN_mpg.shape[1],Nens_plot))
    latsS_mpg_ens = np.zeros((latsS_mpg.shape[0],latsS_mpg.shape[1],Nens_plot))

    latsNdanger_mpg_ens = np.zeros((latsNdanger_mpg.shape[0],latsNdanger_mpg.shape[1],Nens_plot))
    latsSdanger_mpg_ens = np.zeros((latsSdanger_mpg.shape[0],latsSdanger_mpg.shape[1],Nens_plot))
    
    for n_ens in range(Nens_plot):
        num_ens = str(n_ens+1).zfill(3)
        ens_file = 'mpgMF4aacgmv2_from1590_e'+num_ens+'_L13.txt'
        
        folder_coords = folder_base+'coords'+str(2070)+'_bisection/'
        folder_quantity = folder_coords + 'auroral_zones/mpg_ensemble/'
        _, latsN_mpg_ens[:,:,n_ens], latsS_mpg_ens[:,:,n_ens] = CGMlat_bisection(folder_quantity,2070 -mpg_t_offset,2070 ,mpg_folder+'Ensemble/',ens_file)
        
        folder_coords = folder_base+'coords'+str(2070)+'_bisection/'
        folder_quantity = folder_coords + 'danger_zones/mpg_ensemble/'
        _, latsNdanger_mpg_ens[:,:,n_ens], latsSdanger_mpg_ens[:,:,n_ens] = CGMdanger_lat_bisection(folder_quantity,2070 -mpg_t_offset,2070 ,mpg_folder+'Ensemble/',ens_file)
             
        
        
        
        
    #########
    # plots
    #########
    
    
    ##############
    ## plot zones 
    fig = plt.figure(figsize=(10,6))
    
    # north pole
    ax1 = fig.add_subplot(1, 2, 1, projection=ccrs.Orthographic(0, 90))
    
    ax1.set_global() # important (apparently)
    ax1.coastlines()
    ax1.gridlines()    
    

    for i in [0,1]:
        ax1.plot(np.rad2deg(lons_b), latsN2020[:,i], 
             igrf_LS,color='red',
             transform=ccrs.PlateCarree())

        ax1.plot(np.rad2deg(lons_b), latsN_igrf[:,i], 
             igrf_LS,color='blue',
             transform=ccrs.PlateCarree())
        
        ax1.plot(np.rad2deg(lons_b), latsN_ipgp[:,i], 
             ipgp_LS,color='blue',
             transform=ccrs.PlateCarree())
    
        ax1.plot(np.rad2deg(lons_b), latsN_mpg[:,i], 
             mpg_LS,color='blue',
             transform=ccrs.PlateCarree())
    
        for n_ens in range(Nens_plot):
            ax1.plot(np.rad2deg(lons_b), latsN_mpg_ens[:,i,n_ens], 
            mpg_LS,color='blue',
            linewidth=lw_ens,alpha =alpha_ens,
            transform=ccrs.PlateCarree() )

    for city_idx in range(len(city_list_aurorae_N)):
        city_name = city_list_aurorae_N[city_idx]
        city = cities[cities["city_ascii"]==city_name]
        city = city[city["lat"]>0]
        if city.empty:
                continue
        else:
            if city.shape[0]>1:
                city.sort_values('population')
                city = city.iloc[0]
            
            ax1.plot(city["lng"],city["lat"],'o',color='teal',ms =4,transform=ccrs.PlateCarree())
            ax1.text(city["lng"]+off_lon_aurorae_N[city_idx],city["lat"]+off_lat_aurorae_N[city_idx],city_name,fontsize=FSL,color='teal',transform=ccrs.PlateCarree())
              
    ax1.text(-50, 40, str(2020), color='red',fontsize=18,transform=ccrs.PlateCarree())
    ax1.text(95, 63, str(2070), color='blue',fontsize=18,transform=ccrs.PlateCarree())  
    
    ax1.set_title('Northern auroral zone',fontsize=20)
    
    
    ax2 = fig.add_subplot(1, 2, 2, projection=ccrs.Orthographic(0, -90))
    
    ax2.set_global() # important (apparently)
    ax2.coastlines()
    ax2.gridlines(crs=ccrs.PlateCarree(), 
                      linewidth=1, color='black', linestyle=':')    
            
    for i in [0,1]:
        ax2.plot(np.rad2deg(lons_b), latsS2020[:,i], 
             igrf_LS,color='red',
             transform=ccrs.PlateCarree())

        ax2.plot(np.rad2deg(lons_b), latsS_igrf[:,i], 
             igrf_LS,color='blue',
             transform=ccrs.PlateCarree())
        
        ax2.plot(np.rad2deg(lons_b), latsS_ipgp[:,i], 
             ipgp_LS,color='blue',
             transform=ccrs.PlateCarree())
    
        ax2.plot(np.rad2deg(lons_b), latsS_mpg[:,i], 
             mpg_LS,color='blue',
             transform=ccrs.PlateCarree())
    
        for n_ens in range(Nens_plot):
            ax2.plot(np.rad2deg(lons_b), latsS_mpg_ens[:,i,n_ens], 
            mpg_LS,color='blue',
            linewidth=lw_ens,alpha =alpha_ens,
            transform=ccrs.PlateCarree() )

    for city_idx in range(len(city_list_aurorae_S)):
        city_name = city_list_aurorae_S[city_idx]
        city = cities[cities["city_ascii"]==city_name]
        city = city[city["lat"]<0]
        if city.empty:
                continue
        else:
            if city.shape[0]>1:
                city.sort_values('population')
                city = city.iloc[0]
            
            ax2.plot(city["lng"],city["lat"],'o',color='teal',ms =4,transform=ccrs.PlateCarree())
            ax2.text(city["lng"]+off_lon_aurorae_S[city_idx],city["lat"]+off_lat_aurorae_S[city_idx],city_name,fontsize=FSL,color='teal',transform=ccrs.PlateCarree())
     
    ax2.text(50,-85, str(2020), color='red',fontsize=18,transform=ccrs.PlateCarree())
    ax2.text(-10,-68, str(2070), color='blue',fontsize=18,transform=ccrs.PlateCarree())       
    
    ax2.set_title('Southern auroral zone',fontsize=20)
    
    plt.savefig(folder_figures +'polar_zones_ensemble_cities_cartopy.pdf',bbox_inches='tight',pad_inches=0.1)
    
    plt.show(block=False)             
            
    #############
    # danger zones
    
    # America
    fig = plt.figure(figsize=(15,6)) # bigger figure-> smaller writings
    lon0 = -180
    lon1 = -40
    lat0 = 20
    lat1 = 75
    nsquares = 12
    Dlat = abs(lat1 - lat0)/nsquares
    Dlon = abs(lon1 - lon0)/nsquares

    ax = fig.add_subplot(1,3,1,projection=ccrs.NearsidePerspective(central_longitude=-100, central_latitude=40, satellite_height=35785831/8, false_easting=0, false_northing=0))

    ax.set_global() # important (apparently)
    ax.coastlines()
    ax.gridlines()
    ax.add_feature(cfeature.BORDERS, edgecolor='gray')
    
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
    ax.add_feature(states_provinces, edgecolor='gray', linewidth=0.3)
    

        
    for i in [0,1]:
        ax.plot(np.rad2deg(lons_b), latsNdanger_2020[:,i], 
             igrf_LS,color='orangered',
             transform=ccrs.PlateCarree())

        ax.plot(np.rad2deg(lons_b), latsNdanger_igrf[:,i], 
             igrf_LS,color='deepskyblue',
             transform=ccrs.PlateCarree())

        ax.plot(np.rad2deg(lons_b), latsNdanger_ipgp[:,i], 
             ipgp_LS,color='deepskyblue',
             transform=ccrs.PlateCarree())
            
        ax.plot(np.rad2deg(lons_b), latsNdanger_mpg[:,i], 
             mpg_LS,color='deepskyblue',
             transform=ccrs.PlateCarree())

        for n_ens in range(Nens_plot):
            ax.plot(np.rad2deg(lons_b), latsNdanger_mpg_ens[:,i,n_ens], 
            mpg_LS,color='deepskyblue',
            linewidth=lw_ens,alpha =alpha_ens,
            transform=ccrs.PlateCarree() )

    for city_name in city_list_America:
        city = cities[cities["city_ascii"]==city_name]
        city = city[city["lat"]>lat0]
        city = city[city["lat"]<lat1]
        city = city[city["lng"]>lon0]
        city = city[city["lng"]<lon1]
        if city.empty:
                continue
        else:
            if city.shape[0]>1:
                city.sort_values('population')
                city = city.iloc[0]
            
            ax.plot(city["lng"],city["lat"],'o',color='brown',ms =4,transform=ccrs.PlateCarree())
            ax.text(city["lng"]+abs(lon1-lon0)/200,city["lat"]+abs(lat1-lat0)/200,city_name,fontsize=FSL,color='brown',transform=ccrs.PlateCarree())
            
    ax.text(-89,38, str(2020), color='orangered',fontsize=15,transform=ccrs.PlateCarree())
    ax.text(-88,58, str(2070), color='deepskyblue',fontsize=15,transform=ccrs.PlateCarree())    

    
    # Europe
    lon0 = -25
    lon1 = 180
    lat0 = 35
    lat1 = 75
    nsquares = 7
    Dlat = abs(lat1 - lat0)/nsquares
    Dlon = abs(lon1 - lon0)/nsquares
    
    ax1 = fig.add_subplot(1,3,2,projection=ccrs.NearsidePerspective(central_longitude=20, central_latitude=60, satellite_height=35785831/8, false_easting=0, false_northing=0))


    ax1.set_global() # important (apparently)
    ax1.coastlines()
    ax1.gridlines()
    ax1.add_feature(cfeature.BORDERS, edgecolor='gray')
    
    

    for i in [0,1]:
        ax1.plot(np.rad2deg(lons_b), latsNdanger_2020[:,i], 
             igrf_LS,color='orangered',
             transform=ccrs.PlateCarree())

        ax1.plot(np.rad2deg(lons_b), latsNdanger_igrf[:,i], 
             igrf_LS,color='deepskyblue',
             transform=ccrs.PlateCarree())

        ax1.plot(np.rad2deg(lons_b), latsNdanger_ipgp[:,i], 
             ipgp_LS,color='deepskyblue',
             transform=ccrs.PlateCarree())
            
        ax1.plot(np.rad2deg(lons_b), latsNdanger_mpg[:,i], 
             mpg_LS,color='deepskyblue',
             transform=ccrs.PlateCarree())

        for n_ens in range(Nens_plot):
            ax1.plot(np.rad2deg(lons_b), latsNdanger_mpg_ens[:,i,n_ens], 
            mpg_LS,color='deepskyblue',
            linewidth=lw_ens,alpha =alpha_ens,
            transform=ccrs.PlateCarree() )

    for city_name in city_list_Europe:
        city = cities[cities["city_ascii"]==city_name]
        city = city[city["lat"]>lat0]
        city = city[city["lat"]<lat1]
        city = city[city["lng"]>lon0]
        city = city[city["lng"]<lon1]
        if city.empty:
                continue
        else:
            ax1.plot(city["lng"],city["lat"],'o',color='brown',ms =4,transform=ccrs.PlateCarree())
            ax1.text(city["lng"]+abs(lon1-lon0)/200,city["lat"]+abs(lat1-lat0)/200,city_name,fontsize=FSL,color='brown',transform=ccrs.PlateCarree())

    ax1.text(47,67, str(2020), color='orangered',fontsize=15,transform=ccrs.PlateCarree())
    ax1.text(52,48, str(2070), color='deepskyblue',fontsize=15,transform=ccrs.PlateCarree())          
        
    # New Zealand
    lon0 = 100
    lon1 = 180
    lat0 = -60
    lat1 = -10
    nsquares = 10
    Dlat = abs(lat1 - lat0)/nsquares
    Dlon = abs(lon1 - lon0)/nsquares
    
    ax2 = fig.add_subplot(1,3,3,projection=ccrs.NearsidePerspective(central_longitude=140, central_latitude=-50, satellite_height=35785831/2, false_easting=0, false_northing=0))

    
    ax2.set_global() # important (apparently)
    ax2.coastlines()
    ax2.gridlines()
    ax2.add_feature(cfeature.BORDERS, linestyle=':')
    ax2.add_feature(cfeature.BORDERS, edgecolor='gray')
    ax2.add_feature(states_provinces, edgecolor='gray', linewidth=0.3)
    
           

    for i in [0,1]:
        ax2.plot(np.rad2deg(lons_b), latsSdanger_2020[:,i], 
             igrf_LS,color='orangered',
             transform=ccrs.PlateCarree())

        ax2.plot(np.rad2deg(lons_b), latsSdanger_igrf[:,i], 
             igrf_LS,color='deepskyblue',
             transform=ccrs.PlateCarree())

        ax2.plot(np.rad2deg(lons_b), latsSdanger_ipgp[:,i], 
             ipgp_LS,color='deepskyblue',
             transform=ccrs.PlateCarree())
            
        ax2.plot(np.rad2deg(lons_b), latsSdanger_mpg[:,i], 
             mpg_LS,color='deepskyblue',
             transform=ccrs.PlateCarree())

        for n_ens in range(Nens_plot):
            ax2.plot(np.rad2deg(lons_b), latsSdanger_mpg_ens[:,i,n_ens], 
            mpg_LS,color='deepskyblue',
            linewidth=lw_ens,alpha =alpha_ens,
            transform=ccrs.PlateCarree() )
        
    for city_name in city_list_Oceania:
        city = cities[cities["city_ascii"]==city_name]
        city = city[city["lat"]>lat0]
        city = city[city["lat"]<lat1]
        city = city[city["lng"]>lon0]
        city = city[city["lng"]<lon1]
        if city.empty:
                continue
        else:
            ax2.plot(city["lng"],city["lat"],'o',color='brown',ms =4,transform=ccrs.PlateCarree())
            ax2.text(city["lng"]+abs(lon1-lon0)/200,city["lat"]+abs(lat1-lat0)/200,city_name,fontsize=FSL,color='brown',transform=ccrs.PlateCarree())
    
    ax2.text(100,-32, str(2020), color='orangered',fontsize=15,transform=ccrs.PlateCarree())
    ax2.text(120,-52, str(2070), color='deepskyblue',fontsize=15,transform=ccrs.PlateCarree())
    
    plt.suptitle('Vulnerability to extreme space weather events',fontsize=20)
    
    plt.tight_layout()
    plt.savefig(folder_figures + 'Zones_ensemble_NearSideProjection.pdf',bbox_inches='tight',pad_inches=0.1)
    plt.show()        
    

    
else:
   print("Something went wrong")