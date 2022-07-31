#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 20:54:10 2021

@author: stefanomaffei
"""

# modules
import aacgmv2
import numpy as np
import datetime as dt
import os
import cartopy.crs as ccrs
import multiprocessing
from multiprocessing import Process, Pool, Queue
from joblib import Parallel, delayed
import math
from shapely import geometry

import importlib
import sys
from importlib import reload  

# parameters
r_cmb = 3485.0e3
r_a   = 6371.0e3


def angular_distance(lon1,lat1,lon2,lat2):
    # calculate angular distance between two points 1 and 2 the coordinates of which are in deg
    # returns answer in deg
    lon1 = lon1*np.pi/180
    lat1 = lat1*np.pi/180
    lon2 = lon2*np.pi/180
    lat2 = lat2*np.pi/180
    
    dist = math.acos(math.sin(lon1)*math.sin(lon2) + math.cos(lon1)*math.cos(lon2)*math.cos(lat2-lat1))*180/np.pi
    return dist
            
def lat_bisection(magfolder, magfile, year, month, day, lons, target_lat, target_res, lat1, lat2, max_iterations, folder_out):
    '''
    Bisection alogorithm to find the geographic latitude at which the CGM latitude is target_lat within the resolution res_target
    
    INPUT:
        magfile: input geomagnetic field model for aacgmv2
        year, month, day: the temporal instant in which to calculate the cgm coordinates
        lons: the longitudes (in degrees) along which the bisection needs to be performed: a 1d numpy array
        target_lat: target latitude, in degrees
        target_res: desired resolution, in degrees
        lat1, lat2 : initial interval for the bisection method, in degrees
        folder_out :  where to print the final result
    OUTPUT:
        a file containing the latitudes, for each longitude, calculated from the bisection algorithm
        
    '''
    lats = np.zeros(lons.shape) # collecting the final result for each longitude    
    
    # point to a specific magnetic field coefficient model file (magmodel is the default)
    aacgmv2.IGRF_COEFFS = magfolder+magfile
    aacgmv2.wrapper.set_coeff_path(aacgmv2.IGRF_COEFFS,True)      
    
    # will I need to define month and year as well?
    dtime = dt.datetime(year, month, day)
    
    for i in range(lons.shape[0]):
        # northern zone , polar border
        lat_pol_i = lat1
        lat_eq_i  = lat2
        res = 10*target_res
        iteration =0
        while res>=target_res:
            lat_mid_i = 0.5*(lat_pol_i + lat_eq_i)

            if iteration == 0:
                CGMlat_pol_i, CGMlon_pol_i, CGMmlt_pol_i = aacgmv2.get_aacgm_coord(lat_pol_i, lons[i], 0, dtime, method='TRACE')
                CGMlat_eq_i, CGMlon_eq_i, CGMmlt_eq_i = aacgmv2.get_aacgm_coord(lat_eq_i, lons[i], 0, dtime, method='TRACE')
                if np.isnan(CGMlat_eq_i):
                    # in case I calculate undefined values of CGM coordinates
                    CGMlat_eq_i = 0
                if np.isnan(CGMlat_pol_i):
                    CGMlat_pol_i = 0
                
            CGMlat_mid_i, CGMlon_mid_i, CGMmlt_mid_i = aacgmv2.get_aacgm_coord(lat_mid_i, lons[i], 0, dtime, method='TRACE')
            if np.isnan(CGMlat_mid_i):
                    CGMlat_mid_i = 0
                    
            #print('longitude = '+str(lons[i])+'; lat_pol_i = '+ str(lat_pol_i)+'; lat_eq_i = '+ str(lat_eq_i)+'; lat_mid_i = '+ str(lat_mid_i))
            #print('CGMlat_pol_i = '+str(CGMlat_pol_i)+'; CGMlat_eq_i = '+str(CGMlat_eq_i)+'; CGMlat_mid_i = '+str(CGMlat_mid_i))
            
            if (CGMlat_pol_i -target_lat)*(CGMlat_mid_i -target_lat) < 0:
                lat_eq_i = lat_mid_i
                CGMlat_eq_i = CGMlat_mid_i
            elif (CGMlat_eq_i -target_lat)*(CGMlat_mid_i -target_lat) < 0:
                lat_pol_i = lat_mid_i
                CGMlat_pol_i = CGMlat_mid_i
            elif (CGMlat_mid_i -target_lat) == 0:
                lats[i] = lat_mid_i
        
            res = abs(CGMlat_mid_i-target_lat)
            #print('error = '+str(res))
            #print(' ')
            
            iteration = iteration +1
            if iteration >= max_iterations:
                lats[i] = lat_mid_i
                print('lat_bisection: max iterations exceeded')
                break
            if res<target_res:
                lats[i] = lat_mid_i
                
        # use previous calculation as first guess for next iteration
        lat1 = np.amin([lats[i]+10,90])
        lat2 = np.amax([lats[i]-10,-90])
        
    
    np.savetxt(folder_out+'bisection_'+str(target_lat)+'_lats_'+magfile,lats)

    return lats

def intersection_bisection(folder_model1, magfile1, 
                           folder_model2, magfile2, 
                           lons_in, lats1, lats2, 
                           target_CGMlat1, target_CGMlat2,
                           y1, m1, d1, 
                           y2, m2, d2, 
                           target_res, target_res_zero, max_iterations, folder_out):
    '''
    

    Parameters
    ----------
    folder_model1 : TYPE
        folder for magmodel1.
    magfile1 : TYPE
        DESCRIPTION.
    folder_model2 : TYPE
        folder for magmodel2.
    magfile2 : TYPE
        DESCRIPTION.
    lons_in : TYPE
        array of longitudes, in degrees.
    lats1 : TYPE
        array of latitudes for the first curve, in degrees.
    lats2 : TYPE
        array of latitude for the second curve, in degrees.    
    target_CGMlat1 : TYPE
        target_CGM latitude on the first curve, in degrees.
    target_CGMlat2 : TYPE
        target_CGM latitude on the second curve, in degrees.
    y1, m1, d1: TYPE
        year, month, day for CGM calculations on the first curve
    y2, m2, d2: TYPE
        year, month, day for CGM calculations on the second curve
    target_res : TYPE
        target_resolution, in degrees, for the intersection longitude.
    target_res : TYPE
        maximum iterations for latitudinal bisection algorithm.
    target_res_zero : TYPE
        maximum iterations for intersection bisection algorithm.
    folder_out : TYPE
        where to write temporary files
    Returns
    -------
    lons_intersections, lats_intersections : TYPE
        longitude and latitude of the intersection, in degrees

    '''
    
    # function of which to find the zeros is lats1 - lats2
    func = lats1 - lats2
    # prepare an array with entries = -1 where func changes sign, =1 otherwise
    sgn_change = np.sign(func[1:]*func[:-1])
    # find indeces where func changes sign
    idxs = np.where(sgn_change==-1)
    
    if len(idxs[0])==0:
        print('no intersections')
        lons_intersections=math.nan
        lats_intersections=math.nan
        
    else:
        # initialise returned arrays
        lons_intersections = np.zeros(idxs[0].shape)
        lats_intersections = np.zeros(idxs[0].shape)
        
        # bisection between idx and idx+1
        for ilons in range(len(idxs[0])):
            idx = idxs[0][ilons]
            res = 10*target_res
            # initialise
            lon_a = lons_in[idx]
            lon_b = lons_in[idx+1]
            lat1_a = lats1[idx]
            lat2_a = lats2[idx]
            lat1_b = lats1[idx+1]
            lat2_b = lats2[idx+1]
            func_a = func[idx]
            func_b = func[idx+1]
            
            # start with simple linear interpolation
            a1 = (lat1_b - lat1_a)/(lon_b-lon_a)
            b1 = lat1_a - a1*lon_a
            
            a2 = (lat2_b - lat2_a)/(lon_b-lon_a)
            b2 = lat2_a - a2*lon_a
            
            lonx = (b2-b1)/(a1-a2)
            lat1x = a1*lonx + b1
            lat2x = a2*lonx + b2
            
            funcx = lat1x - lat2x
            if funcx<target_res_zero:
                print('solution found by linear interpolation')
                lons_intersections[ilons] = lonx
                lats_intersections[ilons] = 0.5*(lat1x+lat2x)
            else:
                # bisection needed
                if func_a*funcx < 0:
                    lon_b = lonx
                    lat1_b = lat1x
                    lat2_b = lat2x
                elif func_b*funcx < 0:
                    lon_a = lonx
                    lat1_a = lat1x
                    lat2_a = lat2x
                iteration = 0
                
                while res>target_res:
                    lon_mid = 0.5*(lon_a+lon_b)                            
                        
                    my_process3 = Process(target=lat_bisection, args=(folder_model1,magfile1, y1, m1, d1, np.array([lon_mid]), target_CGMlat1, target_res, max(lat1_a,lat2_a), min(lat1_a,lat2_a), max_iterations,folder_out))
                    my_process3.start()
                    my_process3.join() 
                    lat1_mid = np.loadtxt(folder_out+'bisection_'+str(target_CGMlat1)+'_lats_'+magfile1)
                    os.system('rm '+folder_out+'bisection_'+str(target_CGMlat1)+'_lats_'+magfile1)               
        
                    my_process3 = Process(target=lat_bisection, args=(folder_model2,magfile2, y2, m2, d2, np.array([lon_mid]), target_CGMlat2, target_res, max(lat1_a,lat2_a), min(lat1_a,lat2_a), max_iterations,folder_out))
                    my_process3.start()
                    my_process3.join() 
                    lat2_mid = np.loadtxt(folder_out+'bisection_'+str(target_CGMlat1)+'_lats_'+magfile2)
                    os.system('rm '+folder_out+'bisection_'+str(target_CGMlat1)+'_lats_'+magfile2)   
        
                    func_mid = lat1_mid - lat2_mid
        
                    #print('longitude = '+str(lons[i])+'; lat_pol_i = '+ str(lat_pol_i)+'; lat_eq_i = '+ str(lat_eq_i)+'; lat_mid_i = '+ str(lat_mid_i))
                    #print('CGMlat_pol_i = '+str(CGMlat_pol_i)+'; CGMlat_eq_i = '+str(CGMlat_eq_i)+'; CGMlat_mid_i = '+str(CGMlat_mid_i))
                    
                    if func_a*func_mid < 0:
                        lon_b = lon_mid
                        lat1_b = lat1_mid
                        lat2_b = lat2_mid
                    elif func_b*func_mid < 0:
                        lon_a = lon_mid
                        lat1_a = lat1_mid
                        lat2_a = lat2_mid
                    elif func_mid == 0:
                        print('found exact solution')
                        lons_intersections[ilons] = lon_mid
                        lats_intersections[ilons] = 0.5*(lat1_mid+lat2_mid)
                
                    res = func_mid
                    #print('error = '+str(res))
                    #print(' ')
                    
                    iteration = iteration +1
                    if iteration >= max_iterations:
                        lons_intersections[ilons] = lon_mid
                        lats_intersections[ilons] = 0.5*(lat1_mid+lat2_mid)
                        print('intersection_bisection: max iterations exceeded')
                        break
                    if res<target_res_zero:
                        lons_intersections[ilons] = lon_mid
                        lats_intersections[ilons] = 0.5*(lat1_mid+lat2_mid)
                            
    return lons_intersections, lats_intersections




def area_diff_bisection(lons_deg, lons_x, lats_x, latsBminus, latsBplus, pol_eq, r):
    
    poly_list = []
    
    
    if isinstance(lons_x, float) and math.isnan(lons_x): # no intersections: need to check overlaps in a different way
        areas_diff = np.zeros((1,))
        polygon_lons = np.concatenate((lons_deg,  # counterclockwise in the interior rim
                               np.flip(lons_deg), # clockwise on the exterior rim
                               np.array([lons_deg[0]]) # close back to the interior rim
                               ))
        '''
        # this did not work so well...
        if pol_eq ==0:
            lats1 = np.maximum.reduce([abs(latsBplus[:,0]),abs(latsBminus[:,0])])
            if (lats1 == abs(latsBplus[:,0])).all():
                lats2 = np.maximum.reduce([abs(latsBplus[:,1]),abs(latsBminus[:,0])])
            else:
                lats2 = np.maximum.reduce([abs(latsBplus[:,0]),abs(latsBminus[:,1])])
        else:
            lats1 = np.minimum.reduce([abs(latsBplus[:,1]),abs(latsBminus[:,1])])
            if (lats1 == abs(latsBplus[:,1])).all():
                lats2 = np.minimum.reduce([abs(latsBplus[:,0]),abs(latsBminus[:,1])])
            else:
                lats2 = np.minimum.reduce([abs(latsBplus[:,1]),abs(latsBminus[:,0])])            

        lats1 = np.sign(latsBplus[0,0])*lats1
        lats2 = np.sign(latsBplus[0,0])*lats2
        
        polygon_lats = np.concatenate((lats1,  
                       np.flip(lats2), 
                       np.array([lats1[0]]) 
                       ))
        idx = 0 # only one polygon
        areas_diff[idx] = spherical_polygon_area(polygon_lats,polygon_lons,r)
        
        if abs(areas_diff[idx]) > 2*np.pi*r**2:
        # let's not be fancy: if the area is greater than half the surface area of the sphere, flip the arrays to get the correct area
            polygon_lats = np.flip(polygon_lats)
            polygon_lons = np.flip(polygon_lons)
            areas_diff[idx] = spherical_polygon_area(polygon_lats,polygon_lons,r)        
        # need to adjust the sign based on which oval is covering the area calculated above
        if abs(np.sum(latsBplus[:,pol_eq])) > abs(np.sum(latsBminus[:,pol_eq])):
            sgn_area = (-1)**pol_eq
        else:
            sgn_area = -(-1)**pol_eq
        areas_diff[idx] = sgn_area*areas_diff[idx]
        
        # create polygon object (for plotting and checking use)
        poly1 = geometry.Polygon([[polygon_lons[p], polygon_lats[p]] for p in range(len(polygon_lons))])
        poly_list.append(poly1)
        
        '''
        
        # assuming small variations (the polar/equatorial edge
        # of the minus oval never crosses the equatorial/polar edge of the plus oval):
            
        # does not matter if the polygons are circled in the wrong direction, only the difference counts
        area1 = spherical_polygon_area(latsBminus[:,pol_eq],lons_deg,r_a/1000)
        area2 = spherical_polygon_area(latsBplus[:,pol_eq],lons_deg,r_a/1000)

        idx = 0 # only one polygon
        areas_diff[idx] = abs(area1-area2)
        # adjust the sign
        if abs(np.sum(latsBplus[:,pol_eq])) > abs(np.sum(latsBminus[:,pol_eq])):
            sgn_area = (-1)**pol_eq
        else:
            sgn_area = -(-1)**pol_eq
        areas_diff[idx] = sgn_area*areas_diff[idx]

        ##################################################
        # ** I am leaving the poly_list empty in this case
        # ** cartopy is not plotting it anyway, 
        # ** and I was probably doing something wrong
        ##################################################
        
    else: # if there are intersections
        areas_diff = np.zeros(lons_x.shape)
        #northern auroral zones: polward polygons
        for idx in range(len(lons_x)):
            # find left
            ilon = (np.abs(lons_deg-lons_x[idx])).argmin()
            if lons_deg[ilon]<lons_x[idx]:
                ilon = ilon +1
            # right edge
            if idx == len(lons_x) -1:
                idx_next = 0
            else:
                idx_next = idx +1
            ilon_next = (np.abs(lons_deg-lons_x[idx_next])).argmin()
            if lons_deg[ilon_next]>lons_x[idx_next]:
                ilon_next = ilon_next - 1
                
            # check if the polygon crosses the zero meridian
            if lons_x[idx_next] < lons_x[idx]:
                # form the new polygon
                lons1 = np.concatenate((lons_deg[ilon:]-360,lons_deg[:ilon_next]))
                lons2 = np.flip(np.concatenate((lons_deg[ilon:]-360,lons_deg[:ilon_next])))
                lats1 = np.concatenate((latsBminus[ilon:,pol_eq],latsBminus[:ilon_next,pol_eq]))
                lats2 = np.flip(np.concatenate((latsBplus[ilon:,pol_eq],latsBplus[:ilon_next,pol_eq])))
                                
                polygon_lons = np.concatenate((np.array([lons_x[idx]]), 
                                               lons1,
                                               np.array([lons_x[idx_next]]),
                                               lons2
                                              ))
                polygon_lats = np.concatenate((np.array([lats_x[idx]]), 
                                               lats1,
                                               np.array([lats_x[idx_next]]),
                                               lats2
                                              ))
            
            else:    
                # form the new polygon                
                lons1 = lons_deg[ilon:ilon_next]
                lons2 = np.flip(lons_deg[ilon:ilon_next])
                lats1 = latsBminus[ilon:ilon_next,pol_eq]
                lats2 = np.flip(latsBplus[ilon:ilon_next,pol_eq])
                
                polygon_lons = np.concatenate((np.array([lons_x[idx]]), 
                                               lons1,
                                               np.array([lons_x[idx_next]]),
                                               lons2
                                              ))
                polygon_lats = np.concatenate((np.array([lats_x[idx]]), 
                                               lats1,
                                               np.array([lats_x[idx_next]]),
                                               lats2
                                              ))
                
            areas_diff[idx] = spherical_polygon_area(polygon_lats,polygon_lons,r)
            if abs(areas_diff[idx]) > 2*np.pi*r**2:
            # let's not be fancy: if the area is greater than half the surface area of the sphere, flip the arrays to get the correct area
                polygon_lats = np.flip(polygon_lats)
                polygon_lons = np.flip(polygon_lons)
                areas_diff[idx] = spherical_polygon_area(polygon_lats,polygon_lons,r)        
            # need to adjust the sign based on which oval is covering the area calculated above
            if abs(np.sum(lats2)) > abs(np.sum(lats1)):
                sgn_area = (-1)**pol_eq
            else:
                sgn_area = -(-1)**pol_eq
            areas_diff[idx] = sgn_area*areas_diff[idx]
            
            # create polygon object (for plotting and checking use)
            ##################################################
            # ** Actually fails for very small variations.
            # ** I am just going to skip this
            # ** Maybe solve this better some other time
            ##################################################
            #poly1 = geometry.Polygon([[polygon_lons[p], polygon_lats[p]] for p in range(len(polygon_lons))])
            #poly_list.append(poly1)
        
    return areas_diff, poly_list

def calc_aacgm_coord(magfile,year,month,day,lats,lons):
    '''
    
    INPUT:
        magfile = full path of magnetic field coefficient file (see the magmodel for an example)
        t = year of interest
        lats,lons = latitude and longitude (in meshgrid format)
    OUTPUT:
        CGMlats
        CGMlons
        CGMmlt
    '''

    # point to a specific magnetic field coefficient model file (magmodel is the default)
    aacgmv2.IGRF_COEFFS = magfile
    aacgmv2.wrapper.set_coeff_path(aacgmv2.IGRF_COEFFS,True)      
    
    # will I need to define month and year as well?
    dtime = dt.datetime(year, month, day)


    CGMlats = np.zeros(lats.shape)
    CGMlons = np.zeros(lats.shape)
    CGMmlt   = np.zeros(lats.shape)
    
    for i in range(lats.shape[0]):
        for j in range(lats.shape[1]):
            CGMlats[i,j], CGMlons[i,j], CGMmlt[i,j] = np.array(aacgmv2.get_aacgm_coord(np.rad2deg(lats[i,j]), np.rad2deg(lons[i,j]), 0, dtime, method='TRACE'))

    return CGMlats, CGMlons, CGMmlt

def print_aacgm_coords(folder,filename, folder_out,t,lats,lons):
    '''
    wrapper for aacgm_functions.calc_aacgm_coord: it executes the function 
    and prints the results to files.
    
    I wanted to use multiprocessing.Queue to avoid printing things to file, 
    but I could not figure out how to use it, so...
    
    INPUT:
        folder+filename : path and name of the geomagnetic field model file
        t:  year of interest
        lats, lons : meshgrid lat and lon, as needed for aacgm_functions.calc_aacgm_coord
        folder_out: where to put the output files
    '''
    # run the function that calculated aacgm coordinates from lats and lons
    CGMlats, CGMlons, CGMmlt = calc_aacgm_coord(folder+filename,t,1,1,lats,lons)
    # print output to files
    np.savetxt(folder_out+'CGMlats_'+str(t)+'_'+filename,CGMlats)
    np.savetxt(folder_out+'CGMlons_'+str(t)+'_'+filename,CGMlons)
    np.savetxt(folder_out+'CGMmlt_'+str(t)+'_'+filename,CGMmlt)
    

def print_aacgm_coord_parallel(lats,lons,height,year,month,day,model_folder,model_file,folder_out):
    '''
    Parallelized code to calculate AACGM coordinates
    
    INPUT:
        lats,lons = latitude and longitude, in meshgrid form
        height = height above Earth surface, 0 being the radius of the Earth
        dtime = datetime format of the instant of interest
        model_folder = folder where the model file is located (must end with a slash)
        model_file   = model coefficients
        folder_out   = a folder where temporary files can be written (and deleted by the end of the function)
    OUTPUT: 
        CGMlats, CGMlons, CGMmlt
    '''
    #import aacgmv2 as aacgmv2_loc
    #aacgmv2 = reload(aacgmv2) -> error
    #importlib.reload(sys.modules['aacgmv2'])
    
    aacgmv2.IGRF_COEFFS = model_folder+model_file
    aacgmv2.wrapper.set_coeff_path(aacgmv2.IGRF_COEFFS,True)  
        
    CGMlats = np.zeros(lats.shape)
    CGMlons = np.zeros(lats.shape)
    CGMmlts   = np.zeros(lats.shape)
    
    dtime = dt.datetime(year, month, day)
    
    num_cores = multiprocessing.cpu_count()
    
    def calc_aacgm_coords_latwise(lats,lons,height,i,dtime,model_folder,model_file,folder_out):    
        CGMlat = np.zeros(lats.shape[1])
        CGMlon = np.zeros(lats.shape[1])
        CGMmlt   = np.zeros(lats.shape[1])
        for j in range(lats.shape[1]):
            CGMlat[j], CGMlon[j], CGMmlt[j] = np.array(aacgmv2.get_aacgm_coord(np.rad2deg(lats[i,j]), np.rad2deg(lons[i,j]), height, dtime, method='TRACE'))
        # cheap trick to access output: print each longitudinal slice to file    
        np.savetxt(folder_out+'temp_CGMlat_i='+str(i)+'.txt',CGMlat)
        np.savetxt(folder_out+'temp_CGMlon_i='+str(i)+'.txt',CGMlon)
        np.savetxt(folder_out+'temp_CGMmlt_i='+str(i)+'.txt',CGMmlt)
        return 0 #CGMlats, CGMlon, CGMmlt
    
    Parallel(n_jobs=num_cores)(delayed(calc_aacgm_coords_latwise)(lats,lons,0,i,dtime,model_folder,model_file,folder_out) for i in range(lats.shape[0]))
    
    # load the output and rebuild the full matrices
    CGMlats  = np.zeros(lats.shape)
    CGMlons  = np.zeros(lats.shape)
    CGMmlts  = np.zeros(lats.shape)
    
    for i in range(lats.shape[0]):
        CGMlats[i,:] = np.loadtxt(folder_out+'temp_CGMlat_i='+str(i)+'.txt')
        CGMlons[i,:] = np.loadtxt(folder_out+'temp_CGMlon_i='+str(i)+'.txt')
        CGMmlts[i,:] = np.loadtxt(folder_out+'temp_CGMmlt_i='+str(i)+'.txt')
    os.system('rm '+folder_out+'temp_CGM*_i=*txt')  
    
    # print final output to files
    np.savetxt(folder_out+'CGMlats_'+str(year)+'_'+model_file,CGMlats)
    np.savetxt(folder_out+'CGMlons_'+str(year)+'_'+model_file,CGMlons)
    np.savetxt(folder_out+'CGMmlt_'+str(year)+'_'+model_file,CGMmlts)
    
    return 0 #CGMlats, CGMlons, CGMmlts
    

def print_aacgm_coord_parallel_v2(lats,lons,height,year,month,day,model_folder,model_file,folder_out):
    '''
    CAUSES ERRORS
    Parallelized code to calculate AACGM coordinates
    
    INPUT:
        lats,lons = latitude and longitude, in meshgrid form
        height = height above Earth surface, 0 being the radius of the Earth
        dtime = datetime format of the instant of interest
        model_folder = folder where the model file is located (must end with a slash)
        model_file   = model coefficients
        folder_out   = a folder where temporary files can be written (and deleted by the end of the function)
    OUTPUT: 
        CGMlats, CGMlons, CGMmlt
    '''
    CGMlats = np.zeros(lats.shape)
    CGMlons = np.zeros(lats.shape)
    CGMmlts   = np.zeros(lats.shape)
    
    dtime = dt.datetime(year, month, day)
    
    num_cores = multiprocessing.cpu_count()
    
    def calc_aacgm_coords_latwise(lats,lons,height,i,dtime,model_folder,model_file,folder_out):
        aacgmv2.IGRF_COEFFS = model_folder+model_file
        aacgmv2.wrapper.set_coeff_path(aacgmv2.IGRF_COEFFS,True)      
        CGMlat = np.zeros(lats.shape[1])
        CGMlon = np.zeros(lats.shape[1])
        CGMmlt   = np.zeros(lats.shape[1])
        for j in range(lats.shape[1]):
            CGMlat[j], CGMlon[j], CGMmlt[j] = np.array(aacgmv2.get_aacgm_coord(np.rad2deg(lats[i,j]), np.rad2deg(lons[i,j]), height, dtime, method='TRACE'))
        # cheap trick to access output: print each longitudinal slice to file    
        np.savetxt(folder_out+'temp_CGMlat_i='+str(i)+'.txt',CGMlat)
        np.savetxt(folder_out+'temp_CGMlon_i='+str(i)+'.txt',CGMlon)
        np.savetxt(folder_out+'temp_CGMmlt_i='+str(i)+'.txt',CGMmlt)
        return 0 #CGMlats, CGMlon, CGMmlt
        
    # I need to fire up a separate process in a new kernel to be able to change the coefficient file
    def aacgm_process(lats,lons,height,i,dtime,model_folder,model_file,folder_out):
        my_process3 = Process(target=calc_aacgm_coords_latwise, args=(lats,lons,height,i,dtime,model_folder,model_file,folder_out))
        my_process3.start()
        my_process3.join() 
        return 0
    
    Parallel(n_jobs=num_cores)(delayed(aacgm_process)(lats,lons,0,i,dtime,model_folder,model_file,folder_out) for i in range(lats.shape[0]))
    
    # load the output and rebuild the full matrices
    CGMlats  = np.zeros(lats.shape)
    CGMlons  = np.zeros(lats.shape)
    CGMmlts  = np.zeros(lats.shape)
    
    for i in range(lats.shape[0]):
        CGMlats[i,:] = np.loadtxt(folder_out+'temp_CGMlat_i='+str(i)+'.txt')
        CGMlons[i,:] = np.loadtxt(folder_out+'temp_CGMlon_i='+str(i)+'.txt')
        CGMmlts[i,:] = np.loadtxt(folder_out+'temp_CGMmlt_i='+str(i)+'.txt')
    os.system('rm '+folder_out+'temp_CGM*_i=*txt')  
    
    # print final output to files
    np.savetxt(folder_out+'CGMlats_'+str(year)+'_'+model_file,CGMlats)
    np.savetxt(folder_out+'CGMlons_'+str(year)+'_'+model_file,CGMlons)
    np.savetxt(folder_out+'CGMmlt_'+str(year)+'_'+model_file,CGMmlts)
    
    return 0 #CGMlats, CGMlons, CGMmlts
    
    
def read_magmodel_files(magfile):
    '''
    Reads the coefficients from the magmodel files, they should all be in the same format,
    which is the format of the magmodel_1590-2020.txt file in the aacgmv2 directory
    INPUT:
        the magmodel file location
    OUTPUT: 
        header1
        header2
        header3
        header4
        coeffs_MF
        coeffs_SV
    '''
    
    coeffs_MF = []
    coeffs_SV = []

    f = open(magfile,"r")
    header1 = f.readline()
    header2 = f.readline()
    header3 = f.readline()
    line = f.readline()
    header4 = line+''
    x = [list(map(float, line.split()[3:-1]))]
    years = np.squeeze(np.array(x))
    coeffs_MF = np.append(coeffs_MF,years)
    while line:
        line = f.readline()
        if bool(line):
            x = [list(map(float, line.split()[1:]))] 
            coeffs_MF = np.append(coeffs_MF,np.squeeze(np.array(x))[2:-1])
            coeffs_SV = np.append(coeffs_SV,np.squeeze(np.array(x))[-1])
    f.close()
    
    cols = len(years)
    rows = len(coeffs_MF)
    # the transpose below allows continuity with previously written functions
    coeffs_MF = np.transpose(np.reshape(coeffs_MF,(int(rows/cols),cols)))
    coeffs_SV = np.transpose(coeffs_SV)
    
    # because aacgm_v2 is somehow picky about the headers, I need to export them too
    return header1, header2, header3, header4, coeffs_MF, coeffs_SV



def write_magmodel_files(header1, header2, header3, header4, coeffs_MF, coeffs_SV, magfile):
    '''
    Writes magmodel files from given coefficients, in the format required by aacgm_v2
    To be used in conjuction with read_magmodel_files, at least to read the 
    headers from the the magmodel_1590-2020.txt file in the aacgmv2 directory
    
    ******************************************************************
    No control is made on the number of columns and rows: 
        the input files need to have the sizes described below 
        and describe a dataset that spans 1950-2020 every 5 years
    ******************************************************************
    
    INPUT:
        header1
        header2
        header3
        header4
        coeffs_MF: assumed sizes from magmodel_1590-2020.txt is (87,196)
        coeffs_SV: assumed sizes from magmodel_1590-2020.txt is (195,)
        magfile: location and name of the file to be written
    OUTPUT:
        
    '''
    
    f = open(magfile,'w')
    # write the header
    f.write(header1)
    f.write(header2)
    f.write(header3)
    f.write(header4)
    # write the coefficients
    
    # write the coefficients
    Nd = coeffs_MF.shape[1] # number of columns in coeffs_MF
    l=1
    m=0
    cs='g'
    for j in range(1,Nd): # cycle over the columns degree and order of SH
        #print(l,m,cs)
        f.write( cs+'  '+str(l)+'  '+str(m)+'  ' )
        for i in range(coeffs_MF.shape[0]): # cicle over the years
            f.write(str(coeffs_MF[i,j]) + '  ')       
        f.write(str(coeffs_SV[j-1]) + '\n')    # keep the last value of the forecast
        # update indexes
        if l==m:
            if l==0:
                l=1
                m=0
                cs='g'
            elif l>0:
                if cs=='h':
                    m=0
                    l=l+1
                    cs='g'
                elif cs=='g':
                    cs='h'
        elif m<l:
            if m==0:
                m=1
            elif m>0:
                if cs=='g':
                    cs='h'
                elif cs=='h':
                    cs='g'
                    m=m+1
                
    f.close()
    
    return

def spherical_polygon_area(vlat,vlon,r):
    '''
    Routine to calculate the surface area of a spherical polygon on a sphere of radius r.
    Based on Bevis and Cambareri, 1987 : "Computing the area of a spherical polygon of arbitrary shape"
    
    IMPORTANT:
    The  vertices  are  numbered  sequentially  around  the  border  of  the 
    spherical  polygon.  Vertex  1  lies  between  vertex  nv  and  vertex  2. 
    The  user  must  follow  the  convention  whereby  in  moving  around  the 
    polygon  border  in  the  direction  of  increasing  vertex  number  clockwise 
    bends  occur  at  salient  vertices.  A vertex  is  salient  if  the  interior 
    angle'is  less  than  180  degrees. (In  the  case  of  a  convex  polygon 
    this  convention  implies  that  the  vertices  are  numbered  in  clockwise 
    sequence). 
    Two  adjacent  vertices  may  never  be  exactly  180  degrees  apart 
    because  they  could  be  connected  by  infinitely  many  different 
    great  circle  arcs,  and  thus  the  border  of  the  spherical 
    polygon  would  not  be  uniquely  defined. 
    
    IN PRACTICE:
    typically one needs to flip the input vectors vlon and vlat, since vlon is
    commonly 0<vlon<360 in increasing order.

    Parameters
    ----------
    vlat : 1-D array, in degrees
        the list of latitudes of the vertices. These are positive in the North and negative in the South.
    vlon : 1-D array, in degrees
        the list of longitudes of the vertices. These are positive to the East and negative to the West.
    r : the radius of the sphere 

    Returns
    -------
    area : float
        the surface area of the polygon, in units of r^2.

    '''
    sum_ang = 0
    nv = vlat.shape[0]
    for i in range(nv):
        if i == 0:
            flat = vlat[1]
            flon = vlon[1]
            blat = vlat[-1]
            blon = vlon[-1]
        elif i == nv -1:
            flat = vlat[0]
            flon = vlon[0]
            blat = vlat[i-1]
            blon = vlon[i-1]
        else:
            flat = vlat[i+1]
            flon = vlon[i+1]
            blat = vlat[i-1]
            blon = vlon[i-1]
        
        fang = trnsfrm_lon(vlat[i],vlon[i],flat,flon)
        bang = trnsfrm_lon(vlat[i],vlon[i],blat,blon)
        
        fvb=bang-fang 
        if fvb<0:
            fvb = fvb+2*np.pi
            
        sum_ang = sum_ang + fvb
    
    area = (sum_ang - np.pi*(nv-2))* r**2
        
    
    return area

def trnsfrm_lon(plat,plon,qlat,qlon):
    '''
    subroutine to find the 'longitude' of point Q in a geographic coordinate system
    for which point P acts as a 'North Pole'. 

    Parameters
    ----------
    plat : float
        latitude of point P (in degrees).
    plon : float
        longitude of point P (in degress).
    qlat : float
        latitude of point Q (in degrees).
    qlon : float
        longitude of point P (in degress).

    Returns
    -------
    tranlon : float
        in RADIANS.

    '''
    t = np.sin(np.deg2rad(qlon-plon)) * np.cos(np.deg2rad(qlat)) 
    b = np.sin(np.deg2rad(qlat)) * np.cos(np.deg2rad(plat)) - np.cos(np.deg2rad(qlat)) * np.sin(np.deg2rad(plat)) * np.cos(np.deg2rad(qlon-plon)) 
    tranlon = math.atan2(t,b)
    
    return tranlon
    