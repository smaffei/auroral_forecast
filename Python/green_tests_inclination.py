#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 19:10:03 2021

@author: stefanomaffei
"""

# test the Green's function approach on known quantities (Inclination)

import sys
import aacgmv2
import aacgm_functions # my module to implement aacgmv2
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import math
import os
from multiprocessing import Process, Pool, Queue
from scipy.special import roots_legendre, eval_legendre

#from skimage import measure
import matplotlib.ticker as ticker
import matplotlib.colors as colors

from decimal import *

import SH_library

plt.close('all')


# change plot fonts globally   
#plt.rcParams["font.family"] = "Times"
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
plt.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
plt.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

# define coordinates of Sulmona basin location
#colat_SUL = 47.8488
#lon_SUL   = 13.8229 

#Leeds
colat_SUL = 90-53.810380928219175
lon_SUL =-1.5494718924623057 

# define Earth surface and CMB radii
r_a = 6371.0     # Earth's surface
r_c = 3485.0     # CMB


'''
# uniform theta grid:
# define latitude and longitude
lats_lin = np.linspace(-np.pi/2,np.pi/2,num=100)
#lats_lin = np.linspace(50,80,num=70)*np.pi/180.

lons_lin = np.linspace(0,2*np.pi,num=100)

lats, lons = np.meshgrid(lats_lin,lons_lin)
theta = -lats+np.pi/2
'''

# for quadrature-based integrals
nleg = 100
costh, weights = roots_legendre(nleg)
theta_lin = np.arccos(costh)
lats_lin = -theta_lin + np.pi/2

lons_lin = np.linspace(0,2*np.pi,num=100)

theta, lons = np.meshgrid(theta_lin,lons_lin)
lats = -theta + np.pi/2



    
    
folder3   = '../models/magmodel/'
filename3 = 'magmodel_1590-2020.txt'

folder_base = './'
folder_out = folder_base+'aacgmv2_coords/'
folder_models = folder_base+'models/'
folder_figures = folder_base+'figures/'

folder_green = folder_base+'Green_time_magmodel2020_brute_force/'
folder_quantity = folder_green + 'Inclination/'


# for magmodel
t2020 = 2020
dtime2020 = dt.datetime(t2020, 1, 1)

# read ref model
header1, header2, header3, header4, MFcoeffs, SVcoeffs = aacgm_functions.read_magmodel_files(folder3+filename3)

# isolate dipole coeffs
g10_ref = MFcoeffs[-1,1]
g11_ref = MFcoeffs[-1,2]
h11_ref = MFcoeffs[-1,3]

# SV over 2015-2020:
SVcoeffs3 =  (MFcoeffs[-1,1:] - MFcoeffs[-2,1:])/(MFcoeffs[-1,0] - MFcoeffs[-2,0])

# plot reference model intensity in 2020
MFref_mat = SH_library.lin2matCoeffs(MFcoeffs[-1,1:])
MFm_mat = SH_library.lin2matCoeffs(MFcoeffs[-2,1:])


Br_a, Bt_a, Bp_a = SH_library.calcB(MFref_mat,theta,lons,r_a,r_a)
Br_c, Bt_c, Bp_c = SH_library.calcB(MFref_mat,theta,lons,r_a,r_c)

Brm_a, Btm_a, Bpm_a = SH_library.calcB(MFm_mat,theta,lons,r_a,r_a)
Brm_c, Btm_c, Bpm_c = SH_library.calcB(MFm_mat,theta,lons,r_a,r_c)

dBrdt = (Br_c - Brm_c) /5

dBradt = (Br_a - Brm_a) /5
    
F_a = np.sqrt(Br_a**2 + Bt_a**2 + Bp_a**2)
Incl = np.arctan(-Br_a/np.sqrt(Bt_a**2 + Bp_a**2))*180/np.pi 



fig = plt.figure(figsize=(10,6))
Zs = F_a/1000 
Zmax = 70
Zmin = 20
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cs = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs,cmap='coolwarm' 
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )
clb = plt.colorbar(cs
                   ,ticks=[20,30,40,50,60,70]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   )

CS = ax.contour(lons*180/np.pi,lats*180/np.pi,Zs, [30], colors='white')

clb.set_label(r'$F  [\mu T]$', fontsize = 20)
clb.ax.set_xticklabels(['20','30','40','50','60','70'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/F_surf.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()
    
fig = plt.figure(figsize=(10,6))
Zs = dBrdt/1000
Zmax = 30
Zmin = -30
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
clb.set_label(r'$\dot{B_c}$  [$\mu$T/year]', fontsize = 20)
#clb.ax.set_xticklabels(['20','30','40','50','60','70'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/CMB_SV.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()

fig = plt.figure(figsize=(10,6))
Zs = dBradt/1000
Zmax = 220/1000
Zmin = -220/1000
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
clb.set_label(r'$\dot{B_e}$  [$\mu$T/year]', fontsize = 20)
#clb.ax.set_xticklabels(['20','30','40','50','60','70'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/surf_SV.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()



fig = plt.figure(figsize=(10,6))
Zs = Incl
Zmax = 90
Zmin = -90
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cs = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               ,cmap='PiYG' 
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )
clb = plt.colorbar(cs
                   ,ticks=[-90,-45,0,45,90]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   )
clb.set_label(r'$Inclination  [deg]$', fontsize = 20)
clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
#plt.savefig(folder_case + '/CGMlats.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()

#################################
# get reference value of F at SUL
Br_SUL, Bt_SUL, Bp_SUL = SH_library.calcB(MFref_mat,np.array([colat_SUL*np.pi/180]),np.array([lon_SUL*np.pi/180]),r_a,r_a)
F_SUL = np.sqrt(Br_SUL**2 + Bt_SUL**2 + Bp_SUL**2)
H_SUL = np.sqrt(Bt_SUL**2 + Bp_SUL**2)
Incl_SUL = np.arctan(-Br_SUL/H_SUL)*180/np.pi 


Brm_SUL, Btm_SUL, Bpm_SUL = SH_library.calcB(MFm_mat,np.array([colat_SUL*np.pi/180]),np.array([lon_SUL*np.pi/180]),r_a,r_a)
Fm_SUL = np.sqrt(Brm_SUL**2 + Btm_SUL**2 + Bpm_SUL**2)
Hm_SUL = np.sqrt(Btm_SUL**2 + Bpm_SUL**2)
Inclm_SUL = np.arctan(-Brm_SUL/Hm_SUL)*180/np.pi 

# reference
dIdt_SUL = (Incl_SUL - Inclm_SUL)/5

##################################
##################################
# Green Function calculation
##################################
##################################

##############################
# run a few cases explicitly
###############################
l=13
m=1
cs = 'g' # g for the cosine harmonic, h for the sine

if cs=='g':
    cs_idx = 0
else:
    cs_idx = 1    

# name the folder for this case
folder_case = folder_quantity+cs+'_'+str(l)+'_'+str(m)+'/'

if not os.path.exists(folder_case):
    os.makedirs(folder_case)
    
# calculate the index position for the given l,m,cs
# not necessary in a for loop
coeffs_in_l = lambda x : 2*x+1
if m==0:
    m_pos = 1
else:
    m_pos = 2*m    
idx = int( sum(coeffs_in_l(np.linspace(1,l-1,l-1))) + m_pos + cs_idx )

# variations to be imposed:
# this was found to experimentally give a 10% error from the linear answer for
# all l values    
rel_increase = np.linspace(-0.05,0.05,20)*(math.factorial(l)**(1/3))
#rel_increase = np.linspace(-0.01,0.01,20)

Incl_SUL_plus     = np.zeros(rel_increase.shape)
G_Incl_dg = np.zeros(rel_increase.shape)
dg            = np.zeros(rel_increase.shape)

for ir in range(len(rel_increase)):
    # define increment, relative to the 2020 value of the coefficient itself
    #dg[ir] = MFcoeffs[-1,idx] * rel_increase[ir]
    # but some coefficients are zero!!!!!!
    # so let's define it in terms of the rms value of the coefficients with the same l
    rms_g = np.sqrt( sum( MFcoeffs[-1,idx:int( sum(coeffs_in_l(np.linspace(1,l,l))) +1)]**2 ) ) 
    dg[ir] = rms_g * rel_increase[ir]
    
    # modify the coefficient for this case
    MFcoeffs_plus = 1*MFcoeffs
    MFcoeffs_plus[-1,idx] = MFcoeffs[-1,idx]+dg[ir]
    # prepare modified coefficient list for magnetic field calculation    
    MFplus_mat = SH_library.lin2matCoeffs(MFcoeffs_plus[-1,1:])
    
    # recalculate target quantity
    Br_SUL_plus, Bt_SUL_plus, Bp_SUL_plus = SH_library.calcB(MFplus_mat,np.array([colat_SUL*np.pi/180]),np.array([lon_SUL*np.pi/180]),r_a,r_a)
    H_SUL_plus = np.sqrt(Bt_SUL_plus**2 + Bp_SUL_plus**2)
    Incl_SUL_plus[ir] = np.arctan(-Br_SUL_plus/H_SUL_plus)*180/np.pi 
    
    # calculate Green's Function component
    G_Incl_dg[ir] = (Incl_SUL_plus[ir] -Incl_SUL)/dg[ir]

# theoretical expectiation
# calculate the kernels for each B components
MF_g = 0*MFcoeffs[-1,1:int(sum(coeffs_in_l(np.linspace(1,l,l)))+1)]
MF_g[idx-1] = 1 
MF_g_mat = SH_library.lin2matCoeffs(MF_g)

dBrdg, dBtdg, dBpdg = SH_library.calcB(MF_g_mat,np.array([colat_SUL*np.pi/180]),np.array([lon_SUL*np.pi/180]),r_a,r_a)

# calculate theoretical linear green's functions
dIdg = (180./np.pi)*(Bt_SUL*Br_SUL*dBtdg + Bp_SUL*Br_SUL*dBpdg - dBrdg*H_SUL**2 ) / (H_SUL*F_SUL**2)

# theroetical, linearised 2015-2020 variation
dIt_lin = (180./np.pi)*(Bt_SUL*Br_SUL*(Bt_SUL-Btm_SUL)/5 + Bp_SUL*Br_SUL*(Bp_SUL-Bpm_SUL)/5 - H_SUL**2 * (Br_SUL-Brm_SUL)/5) / (H_SUL*F_SUL**2)


# calculate theoretical linear variation for given increments
#dg_1_0 = MFcoeffs[-1,idx] * rel_increase
dI_lin = dIdg * dg


# plot it

fig,axs = plt.subplots(1,2,figsize=(10,6))

FS = 12

coeff_label = cs+'_{'+str(l)+'}^{'+str(m)+'}'

ax = axs[0]
g = (MFcoeffs[-1,idx]+dg)/1000
ax.plot(g,dI_lin + Incl_SUL,'ko-', label='linearised')
ax.plot(g,Incl_SUL_plus,'ro-', label='finite differences')

ax.axvline(x = MFcoeffs[-1,idx]/1000, ymin = 0, ymax =90
           , color='b',linestyle='--', linewidth=2)

ax.tick_params(axis='x', which = 'major', labelsize=FS)
ax.tick_params(axis='y', which = 'major', labelsize=FS)

ax.set_xlabel(r'$'+coeff_label+' [\mu T]$',fontsize=FS) 
ax.set_ylabel(r'Inclination [deg]',fontsize=FS) 

ax.legend(loc='best'
          ,fontsize=12
          ,framealpha = 1)

ax = axs[1]
ax.plot(rel_increase,dIdg*np.ones(rel_increase.shape),'ko-')
ax.plot(rel_increase,G_Incl_dg,'ro-')

ax.axvline(x = 0, ymin = 0, ymax =90
           , color='b',linestyle='--', linewidth=2)

ax.axhline(y = dIdg + 0.1*dIdg
           , color='k', linestyle='--', linewidth=2)
ax.annotate('10 %'
            ,(0.9*np.max(rel_increase),dIdg + abs(0.11*dIdg))
            ,fontsize = 12)
ax.annotate('10 %'
            ,(0.9*np.max(rel_increase),dIdg - abs(0.15*dIdg))
            ,fontsize = 12)
ax.axhline(y = dIdg - 0.1*dIdg
           , color='k', linestyle='--', linewidth=2)

ax.tick_params(axis='x', which = 'major', labelsize=FS)
ax.tick_params(axis='y', which = 'major', labelsize=FS)

ax.set_xlabel(r'$\delta '+coeff_label+' / r_{'+str(l)+'}$',fontsize=FS) 
ax.set_ylabel(r'slope',fontsize=FS) 

plt.tight_layout()
plt.savefig(folder_case + 'slopes_'+cs+'_'+str(l)+'_'+str(m)+'.pdf',bbox_inches='tight',pad_inches=0.1)

plt.show()



################################
# for loop over all coefficients
#################################

# initialise indexes
cs = 'g'
l  = 0
m  = 0

# initialize output quantities
G_Incl_dg = np.zeros(MFcoeffs[-1,1:].shape)
Lt = 13
dIdg = np.zeros(MFcoeffs[-1,1:].shape)
dIdBr = np.zeros(theta.shape)
dIdBra = np.zeros(theta.shape)

dFdg = np.zeros((MFcoeffs[-1,1:].shape[0],theta.shape[0],theta.shape[1]))

#dBrdbc = np.zeros(theta.shape)

# loop 
for idx in range(1,MFcoeffs[-1,1:].shape[0]+1):    
    # update SH coeff label
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
                
    coeff_label = cs+'_'+str(l)+'^'+str(m)
    
    # relative increase of the current Gauss coefficient:
    # experimental value to have Green's functions less 1% away from the theoretical linear kernel
    # and a 1% variation on the l=1 coefficients
    #
    # initial version was  0.01*math.factorial(l), but the increase ended up being too much for higher coefficients
    # even np.sqrt(math.factorial(l)) was big
    rel_increase = 0.01*(math.factorial(l)**(1/3))

    # define increment
    #dg = MFcoeffs[-1,idx] * rel_increase
    # use the rms of the coefficients at this l to define the increments
    rms_g = np.sqrt( sum( MFcoeffs[-1,idx:int( sum(coeffs_in_l(np.linspace(1,l,l))) +1)]**2 ) ) 
    dg = rms_g * rel_increase
    # modify the coefficient for this case
    MFcoeffs_plus = 1*MFcoeffs
    MFcoeffs_plus[-1,idx] = MFcoeffs[-1,idx]+dg
    # prepare modified coefficient list for magnetic field calculation    
    MFplus_mat = SH_library.lin2matCoeffs(MFcoeffs_plus[-1,1:])
    
    # recalculate target quantity
    Br_SUL_plus, Bt_SUL_plus, Bp_SUL_plus = SH_library.calcB(MFplus_mat,np.array([colat_SUL*np.pi/180]),np.array([lon_SUL*np.pi/180]),r_a,r_a)
    H_SUL_plus = np.sqrt(Bt_SUL_plus**2 + Bp_SUL_plus**2)
    Incl_SUL_plus = np.arctan(-Br_SUL_plus/H_SUL_plus)*180/np.pi 
    
    # calculate Green's Function component
    G_Incl_dg[idx-1] = (Incl_SUL_plus -Incl_SUL)/dg

    # theoretical expectiation
    # calculate the kernels for each B components
    MF_g = 0*MFcoeffs[-1,1:int(sum(coeffs_in_l(np.linspace(1,l,l)))+1)]
    MF_g[idx-1] = 1 
    MF_g_mat = SH_library.lin2matCoeffs(MF_g)
    
    dBrdg, dBtdg, dBpdg = SH_library.calcB(MF_g_mat,np.array([colat_SUL*np.pi/180]),np.array([lon_SUL*np.pi/180]),r_a,r_a)
    # calculate theoretical linear green's functions
    dIdg[idx-1] = (180./np.pi)*(Bt_SUL*Br_SUL*dBtdg + Bp_SUL*Br_SUL*dBpdg - dBrdg*H_SUL**2 ) / (H_SUL*F_SUL**2)
    # calculate green's function of I with respect to Br
    # dI/ dBr = sum( dI/dglm * dglm/dBr ) = sum( dI/dglm * ( Br(glm=1) )^-1 )    
    if cs == 'g':
        SH = SH_library.SchmidtSH(l,m,theta,lons,'c')
        dSHdt = SH_library.DthSchmidtSH(l,m,theta,lons,'c')
        dSHdp = SH_library.DphSchmidtSH(l,m,theta,lons,'c')
    else:
        SH = SH_library.SchmidtSH(l,m,theta,lons,'s')
        dSHdt = SH_library.DthSchmidtSH(l,m,theta,lons,'s')
        dSHdp = SH_library.DphSchmidtSH(l,m,theta,lons,'s')

    norm = 4*np.pi/(2*l+1)
    dgdBr_map =  (1/(l+1)) *  (r_c/r_a)**(l+2) * SH / norm
    dIdBr = dIdBr + dIdg[idx-1] * dgdBr_map

    dgdBra_map =  (1/(l+1)) *  (r_a/r_a)**(l+2) * SH / norm
    dIdBra = dIdBra + dIdg[idx-1] * dgdBra_map
        
    dFdg[idx-1,:,:] = np.divide( np.multiply(Br_a,(l+1)*SH) - np.multiply(Bt_a,dSHdt)+np.multiply(Bp_a,np.divide(dSHdp,np.sin(theta))),F_a)



# average over longitudes


dIdt_profile =np.zeros(theta[0,:-1].shape)
dIdt_surf_profile =np.zeros(theta[0,:-1].shape)



check = 0
dIdt_tot=0
dIdt_surf_tot=0

for ith in range(len(theta[0,:])-1):
    for iph in range(len(lons[:,0])-1):
        dth = abs(theta[iph,ith+1]-theta[iph,ith])
        dph = abs(lons[iph+1,ith]-lons[iph,ith])
        thav = (theta[iph,ith+1]+theta[iph,ith])/2  
        phav = (lons[iph+1,ith]+lons[iph,ith])/2  
        
        dIdt_profile[ith] = dIdt_profile[ith] + dBrdt[iph,ith]*dIdBr[iph,ith]*dph*np.sin(thav)
        dIdt_surf_profile[ith] = dIdt_surf_profile[ith] + dBradt[iph,ith]*dIdBra[iph,ith]*dph*np.sin(thav)

        #check= check + dth*dph*np.sin(thav) # surfce area
        check= check + dth*dph*np.sin(thav)*np.cos(thav)**13*np.sin(phav)**13 # should be zero
        #check= check + dth*dph*np.sin(thav)*np.cos(thav)**2*np.sin(phav)**2 # should be 2*np.pi/3


#dIdt_tot = np.sum(dIdtdth_profile)    
#dIdt_surf_tot = np.sum(dIdtdth_surf_profile)    
#####################################
# some plots
######################################



fig = plt.figure(figsize=(10,6))
Zs = np.divide(1,dFdg[0,:,:])
Zmax = 1
Zmin = -1
#Zmax = 10
#Zmin = -10
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               ,cmap='PiYG' 
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )
plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)
             
clb = plt.colorbar(cf
                   #,ticks=[-90,-45,0,45,90]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   ,format='%.1e'
                   )
clb.set_label(r'$\partial F/\partial g$ [nT/nT]', fontsize = 20)
#clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
#plt.savefig(folder_quantity + '/dIdBr.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()




fig = plt.figure(figsize=(10,6))

Zs = dIdBr
Zmax = np.nanmax(abs(Zs))
Zmin = -Zmax
#Zmax = 10
#Zmin = -10
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               #,cmap='PiYG' 
               ,cmap='PuOr_r'
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )
plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)
             
clb = plt.colorbar(cf
                   #,ticks=[-90,-45,0,45,90]
                   ,ticks=[-1e-4, 0, 1e-4]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   ,format='%.1e'
                   )
clb.set_label(r'$\partial I/\partial B_c$ (numerical solution) [deg/nT]', fontsize = 20)
#clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/dIdBr.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()



fig = plt.figure(figsize=(10,6))
Zs = dIdBra
Zmax = np.nanmax(abs(Zs))
Zmin = -Zmax
#Zmax = 10
#Zmin = -10
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               ,cmap='PiYG' 
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )
plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)
             
clb = plt.colorbar(cf
                   #,ticks=[-90,-45,0,45,90]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   ,format='%.1e'
                   )
clb.set_label(r'$\partial I/\partial B_r$ [deg/nT]', fontsize = 20)
#clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/dIdBra.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()
    


# cmb
fig = plt.figure(figsize=(10,6))
Zs = dIdBr*dBrdt 
Zmax = np.nanmax(abs(Zs))
Zmin = -Zmax
#Zmax = 10
#Zmin = -10
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               ,cmap='PiYG' 
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )
plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)
             
clb = plt.colorbar(cf
                   #,ticks=[-90,-45,0,45,90]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   ,format='%.1e'
                   )
clb.set_label(r'$\partial I/\partial t$ (CMB) [deg/year]', fontsize = 20)
#clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/dIdt.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()


# earth's surface
fig = plt.figure(figsize=(10,6))
Zs = dIdBra*dBradt 
Zmax = np.nanmax(abs(Zs))
Zmin = -Zmax
#Zmax = 10
#Zmin = -10
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               ,cmap='PiYG' 
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )
plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)
             
clb = plt.colorbar(cf
                   #,ticks=[-90,-45,0,45,90]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   ,format='%.1e'
                   )
clb.set_label(r'$\partial I/\partial t$ (Earth''s surface) [deg/year]', fontsize = 20)
#clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/dIdt_surf.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()


# integrated in longitude
# cmb
fig,ax = plt.subplots(figsize=(5,8))
ax.set_ylabel(r'lat [deg]',fontsize=FS)
ax.set_xlabel('deg / year ', fontsize=FS)
#ax.set_ylim([77,81])
#ax.set_xlim([1900,2100])
ax.tick_params('x',labelsize=FS)

ax.plot(dIdt_profile,90-180*theta[0,:-1]/np.pi,color='tab:blue')

ax.axvline(x=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')

#ax.legend(fontsize=12,loc='upper right')

ax.tick_params('y',labelsize=FS)
ax.yaxis.get_offset_text().set_fontsize(12)
#ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
#ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
plt.title(r'dIdt (from CMB variations)',fontsize=FS)
plt.savefig(folder_quantity + '/dIdt_lats.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()



# integrated in longitude
# Earth's surface
fig,ax = plt.subplots(figsize=(5,8))
ax.set_ylabel(r'lat [deg]',fontsize=FS)
ax.set_xlabel('deg / year ', fontsize=FS)
#ax.set_ylim([77,81])
#ax.set_xlim([1900,2100])
ax.tick_params('x',labelsize=FS)

ax.plot(dIdt_surf_profile,90-180*theta[0,:-1]/np.pi,color='tab:blue')

ax.axvline(x=0, color='k',alpha=0.6, linewidth=0.8,linestyle='--')

#ax.legend(fontsize=12,loc='upper right')

ax.tick_params('y',labelsize=FS)
ax.yaxis.get_offset_text().set_fontsize(12)
#ax.axvline(x=1859, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
#ax.axvline(x=2019, color='k',alpha=0.6, linewidth=0.8,linestyle='--')
plt.title(r'dIdt (from surface variations)',fontsize=FS)
plt.savefig(folder_quantity + '/dIdt_surf_lats.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()


##########
# plot it
fig,ax = plt.subplots(figsize=(10,6))

FS = 20

coeff_label = cs+'_'+str(l)+'^'+str(m)

g = (MFcoeffs[-1,idx]+dg)/1000
ax.plot(G_Incl_dg,'ro-', label='numerical solution')
ax.plot(dIdg,'k+',markersize=15, label='analytical solution')

ax.tick_params(axis='x', which = 'major', labelsize=FS)
ax.tick_params(axis='y', which = 'major', labelsize=FS)

ax.set_xlabel(r'Spherical Harmonic degree',fontsize=FS) 
ax.set_ylabel(r'$\partial I / \partial \beta_l^m$ [deg/ nT]',fontsize=FS) 

ax.legend(loc='best'
          ,fontsize=15
          ,framealpha = 1)

plt.grid(True)

# redefine the x axis ticks and labels
l_ticks        = np.zeros((13,))
l_ticks_labels = []
for il in range(13):
    l_ticks[il] = int(sum(coeffs_in_l(np.linspace(1,il,il))))
    l_ticks_labels.append(str(il+1))
    
ax.set_xticks(l_ticks)
ax.set_xticklabels(l_ticks_labels)

plt.savefig(folder_quantity + 'Inclination_Green.pdf',bbox_inches='tight',pad_inches=0.1)

plt.show()



# create and plot it in matrix form

GImatrix  = np.zeros((Lt,2*Lt+1))
GImatrix[:,:] = math.nan

SVmatrix  = np.zeros((Lt,2*Lt+1))
SVmatrix[:,:] = math.nan

idx = 0
# fill it with positive variations
for l in range(1,Lt+1,1):
    for m in range(0,l+1,1):
        for cs in ['g','h']:
            if cs =='g':
                GImatrix[l-1,Lt+m] = G_Incl_dg[idx]
                SVmatrix[l-1,Lt+m] = SVcoeffs3[idx]
                idx=idx+1
            else:
                if m==0:
                    continue
                else:
                    GImatrix[l-1,Lt-m] = G_Incl_dg[idx] 
                    SVmatrix[l-1,Lt-m] = SVcoeffs3[idx]
                    idx=idx+1
                    
dIdtmatrix = GImatrix * SVmatrix

fig,ax = plt.subplots(1,1,figsize=(8,4.5))
FS = 12
FT = 15
FST = 18


# find min and max of colorbar:
    
# we manually adjust for the scientific notation 
   
Zmax1 = np.nanmax(abs(GImatrix))
expmax1=math.log10(Zmax1)
if expmax1>0:
    expmax = math.ceil(expmax1)
else :
    expmax = math.floor(expmax1)

Z = GImatrix/10**expmax

Zmax = np.nanmax(abs(Z))
Zmin = -Zmax

cf = ax.matshow(Z
                #,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale = 1, linthresh = 2*10**(-4))
                #,cmap = 'RdYlBu_r'
                ,cmap = 'seismic'
                )

ax.xaxis.set_ticks_position('bottom')

m_ticks = np.arange(0,2*Lt+1,1) # need the +1...simply obsciene
m_labels = list(map(str, abs(np.arange(-Lt,Lt+1,1))))
ax.set_xticks(m_ticks)
ax.set_xticklabels(m_labels)

l_ticks = np.arange(0,Lt,1) # need the +1...simply obsciene
l_labels = list(map(str, np.arange(1,Lt+1,1)))
ax.set_yticks(l_ticks)
ax.set_yticklabels(l_labels)

ax.tick_params(axis='x', which = 'major', labelsize=FS)
ax.tick_params(axis='y', which = 'major', labelsize=FS)

ax.set_xlabel(r'$\leftarrow \ h_l^m\qquad$ order, $m$ $\qquad g_l^m \ \rightarrow$',fontsize=FS) 
ax.set_ylabel(r'degree, $l$',fontsize=FS)

ax.tick_params(labeltop=True, labelright=True, top=True, right=True)


#ax.set_title(r'$\frac{\partial I}{\partial \beta_l^m}$', fontsize=FT, y=1.13)    

clb = fig.colorbar(cf, shrink=0.95, pad=0.1
                   #, ticks = [Zmin, Zmin/2, 0, Zmax/2, Zmax]
                   #, ticks = [-5/1000, -2.5/1000, 0, 2.5/1000, 5/1000]
                   )

#cb_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
#clb = fig.colorbar(cf, cax=cb_ax)
#clb.formatter.set_powerlimits((0, 0))

clb.ax.tick_params(labelsize=FS)
clb.set_label(r'deg/ nT $\times 10^{'+str(expmax)+'} $',fontsize=FS)
clb.formatter.set_powerlimits((0, 0))
clb.update_ticks()

plt.suptitle(r'Inclination sensitivity to Gauss coefficients, $\frac{\partial I}{\partial \beta_l^m}$', fontsize=FST)
plt.tight_layout()
                    
plt.savefig(folder_quantity + 'Inclination_Green_matrix.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show(block=False)  


##### SV coeffs (2020-2015) in the same format:

    
fig,ax = plt.subplots(1,1,figsize=(8,4.5))
FS = 12
FT = 15
FST = 18


# find min and max of colorbar:
    
# we manually adjust for the scientific notation 
   
Zmax1 = np.nanmax(abs(SVmatrix))
expmax1=math.log10(Zmax1)
if expmax1>0:
    expmax = math.ceil(expmax1)
else :
    expmax = math.floor(expmax1)

Z = SVmatrix/10**expmax
Zmax = np.nanmax(abs(Z))
Zmin = -Zmax

cf = ax.matshow(Z
                #,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale = 1, linthresh = 2*10**(-4))
                #,cmap = 'RdYlBu_r'
                ,cmap = 'seismic'
                )

ax.xaxis.set_ticks_position('bottom')

m_ticks = np.arange(0,2*Lt+1,1) # need the +1...simply obsciene
m_labels = list(map(str, abs(np.arange(-Lt,Lt+1,1))))
ax.set_xticks(m_ticks)
ax.set_xticklabels(m_labels)

l_ticks = np.arange(0,Lt,1) # need the +1...simply obsciene
l_labels = list(map(str, np.arange(1,Lt+1,1)))
ax.set_yticks(l_ticks)
ax.set_yticklabels(l_labels)

ax.tick_params(axis='x', which = 'major', labelsize=FS)
ax.tick_params(axis='y', which = 'major', labelsize=FS)

ax.set_xlabel(r'$\leftarrow \ h_l^m\qquad$ order, $m$ $\qquad g_l^m \ \rightarrow$',fontsize=FS) 
ax.set_ylabel(r'degree, $l$',fontsize=FS)

ax.tick_params(labeltop=True, labelright=True, top=True, right=True)


#ax.set_title(r'$\frac{\partial I}{\partial \beta_l^m}$', fontsize=FT, y=1.13)    

clb = fig.colorbar(cf, shrink=0.95, pad=0.1
                   #, ticks = [Zmin, Zmin/2, 0, Zmax/2, Zmax]
                   #, ticks = [-5/1000, -2.5/1000, 0, 2.5/1000, 5/1000]
                   )

#cb_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
#clb = fig.colorbar(cf, cax=cb_ax)
#clb.formatter.set_powerlimits((0, 0))

clb.ax.tick_params(labelsize=FS)
clb.set_label(r'nT/ yr $\times 10^{'+str(expmax)+'} $',fontsize=FS)
clb.formatter.set_powerlimits((0, 0))
clb.update_ticks()

plt.suptitle(r'SV during 2015-2020, $\frac{\partial \beta_l^m}{\partial t}$', fontsize=FST)
plt.tight_layout()
                    
plt.savefig(folder_quantity + 'SV_matrix.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show(block=False)  



##### contributions to dIdt
fig,ax = plt.subplots(1,1,figsize=(8,4.5))
FS = 12
FT = 15
FST = 18


# find min and max of colorbar:
    
# we manually adjust for the scientific notation 
   
Zmax1 = np.nanmax(abs(dIdtmatrix))
expmax1=math.log10(Zmax1)
if expmax1>0:
    expmax = math.ceil(expmax1)
else :
    expmax = math.floor(expmax1)

Z = dIdtmatrix/10**expmax
Zmax = np.nanmax(abs(Z))
Zmin = -Zmax

cf = ax.matshow(Z
                #,norm=colors.SymLogNorm(vmin=Zmin, vmax=Zmax, linscale = 1, linthresh = 2*10**(-4))
                #,cmap = 'RdYlBu_r'
                ,cmap = 'seismic'
                )

ax.xaxis.set_ticks_position('bottom')

m_ticks = np.arange(0,2*Lt+1,1) # need the +1...simply obsciene
m_labels = list(map(str, abs(np.arange(-Lt,Lt+1,1))))
ax.set_xticks(m_ticks)
ax.set_xticklabels(m_labels)

l_ticks = np.arange(0,Lt,1) # need the +1...simply obsciene
l_labels = list(map(str, np.arange(1,Lt+1,1)))
ax.set_yticks(l_ticks)
ax.set_yticklabels(l_labels)

ax.tick_params(axis='x', which = 'major', labelsize=FS)
ax.tick_params(axis='y', which = 'major', labelsize=FS)

ax.set_xlabel(r'$\leftarrow \ h_l^m\qquad$ order, $m$ $\qquad g_l^m \ \rightarrow$',fontsize=FS) 
ax.set_ylabel(r'degree, $l$',fontsize=FS)

ax.tick_params(labeltop=True, labelright=True, top=True, right=True)


#ax.set_title(r'$\frac{\partial I}{\partial \beta_l^m}$', fontsize=FT, y=1.13)    

clb = fig.colorbar(cf, shrink=0.95, pad=0.1
                   #, ticks = [Zmin, Zmin/2, 0, Zmax/2, Zmax]
                   #, ticks = [-5/1000, -2.5/1000, 0, 2.5/1000, 5/1000]
                   )

#cb_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
#clb = fig.colorbar(cf, cax=cb_ax)
#clb.formatter.set_powerlimits((0, 0))

clb.ax.tick_params(labelsize=FS)
clb.set_label(r'deg/ yr $\times 10^{'+str(expmax)+'} $',fontsize=FS)
clb.formatter.set_powerlimits((0, 0))
clb.update_ticks()

plt.suptitle(r'dIdt during 2015-2020, $\frac{\partial I}{\partial t}$', fontsize=FST)
plt.tight_layout()
                    
plt.savefig(folder_quantity + 'SV_matrix.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show(block=False)     

#########
# dI/dBr
########

G_incl_dBr_c = np.zeros(Br_c.shape)
GX = np.zeros(Br_c.shape)
GY = np.zeros(Br_c.shape)
GZ = np.zeros(Br_c.shape)

Gr = np.zeros(Br_c.shape)
Gth = np.zeros(Br_c.shape)
Gph = np.zeros(Br_c.shape)

rho = r_c/r_a



# this increment refers to the last gauss coefficient being modified in the above loop
Br_a_plus, Bt_a_plus, Bp_a_plus = SH_library.calcB(MFplus_mat,theta,lons,r_a,r_a)
Br_c_plus, Bt_c_plus, Bp_c_plus = SH_library.calcB(MFplus_mat,theta,lons,r_a,r_c)

### I thiknk someting is wrong in this calculation
for ilon in range(len(lons_lin)):
    for ilat in range(len(lats_lin)):
        # theoretical Green's function
        
        ###########################
        # johnson + Constable 1997
        #rdots = np.sin(theta[ilon,ilat])*np.sin(np.deg2rad(colat_SUL)) * np.cos(np.deg2rad(lon_SUL)-lons[ilon,ilat]) + np.cos(theta[ilon,ilat])*np.cos(np.deg2rad(colat_SUL))
        #sdotx = - np.cos(theta[ilon,ilat])*np.sin(np.deg2rad(colat_SUL)) * np.cos(np.deg2rad(lon_SUL)-lons[ilon,ilat]) + np.cos(theta[ilon,ilat])*np.sin(np.deg2rad(colat_SUL))
        #sdoty = np.sin(theta[ilon,ilat])*np.sin(np.deg2rad(lon_SUL)-lons[ilon,ilat])
        rdots = np.cos(theta[ilon,ilat])*np.cos(np.deg2rad(colat_SUL))+np.sin(theta[ilon,ilat])*np.sin(np.deg2rad(colat_SUL))*np.cos(np.deg2rad(lon_SUL)-lons[ilon,ilat])
        sdotx = -( np.cos(np.deg2rad(colat_SUL))*np.sin(theta[ilon,ilat])*np.cos(np.deg2rad(lon_SUL)-lons[ilon,ilat]) - np.sin(np.deg2rad(colat_SUL))*np.cos(theta[ilon,ilat]) )
        sdoty = -np.sin(theta[ilon,ilat])*np.sin(np.deg2rad(lon_SUL)-lons[ilon,ilat])
        
        R = np.sqrt(1-2*rdots*rho + rho**2)
        T = 1+R-rdots*rho
        GX[ilon,ilat] = -(1+2*R-rho**2)*rho**3 *sdotx / (R**3 * T) /(4*np.pi)
        GY[ilon,ilat] = -(1+2*R-rho**2)*rho**3 *sdoty / (R**3 * T) /(4*np.pi)
        GZ[ilon,ilat] = (rho**2 - rho**2*(1-rho**2)/R**3)/(4*np.pi)
        
        #################################
        # LEt's follow Hammer+Finlay 2019 (which is Gubbins+Roberts, 1983)
        mu = np.cos(theta[ilon,ilat])*np.cos(np.deg2rad(colat_SUL))+np.sin(theta[ilon,ilat])*np.sin(np.deg2rad(colat_SUL))*np.cos(np.deg2rad(lon_SUL)-lons[ilon,ilat])
        h = r_c/r_a
        R = np.sqrt(r_a**2+r_c**2-2*r_a*r_c*mu)
        f = R/r_a
        
        dNcdmu = (h/(4*np.pi)) * ( (1-2*h*mu+3*h**2)/f**3 + mu/(f*(f+h-mu)) -1/(1-mu) )

        Gr[ilon,ilat]  =  h**2*(1-h**2)/f**3  /(4*np.pi)
        Gth[ilon,ilat] = -dNcdmu * ( np.cos(np.deg2rad(colat_SUL))*np.sin(theta[ilon,ilat])*np.cos(np.deg2rad(lon_SUL)-lons[ilon,ilat]) - np.sin(np.deg2rad(colat_SUL))*np.cos(theta[ilon,ilat]) )
        Gph[ilon,ilat] = dNcdmu * ( np.sin(theta[ilon,ilat])*np.sin(np.deg2rad(lon_SUL)-lons[ilon,ilat]) )


    
# johnson + Constable 1997
#GIc = (180./np.pi)*( H_SUL*GZ - (-Br_SUL/H_SUL)*(-Bt_SUL*GX+Bp_SUL*GY) ) / (H_SUL**2 * Br_SUL**2)
GIc = (180./np.pi)*(-Bt_SUL*Br_SUL*GX + Bp_SUL*Br_SUL*GY + GZ*H_SUL**2 ) / (H_SUL*F_SUL**2)
GDc = (180./np.pi)*(-Bt_SUL*GY - Bp_SUL*GX)/H_SUL**2
GFc = (1/F_SUL)*(-Br_SUL*GZ - Bt_SUL*GX + Bp_SUL*GY)

# Gubbins+Roberts, 1983
GIs = (180./np.pi)*(Bt_SUL*Br_SUL*Gth + Bp_SUL*Br_SUL*Gph - Gr*H_SUL**2 ) / (H_SUL*F_SUL**2)
GDs = (180./np.pi)*(-Bt_SUL*Gph + Bp_SUL*Gth)/H_SUL**2
GFs = (1/F_SUL)*(Br_SUL*Gr + Bt_SUL*Gth + Bp_SUL*Gph)



###########################################################
# check the green's functions against the actual variation
check = 0

dI_dBr_check = 0
dI_dBr_GR1983 = 0
dI_dBr_JC1993 = 0

dBra_dBr_GR1983 = 0
dBta_dBr_GR1983 = 0
dBpa_dBr_GR1983 = 0

dX_dBr_JC1993 = 0
dY_dBr_JC1993 = 0
dZ_dBr_JC1993 = 0

dIdt_tot = 0
dIdt_surf_tot=0

g10_check=0
h11_check=0

g10dot_check =0

# surface integrals:
'''                   
for ilon in range(len(lons_lin)-1):
    for ilat in range(len(lats_lin)-1):
        dth = abs(theta[ilon,ilat+1]-theta[ilon,ilat])
        thav = (theta[ilon,ilat+1]+theta[ilon,ilat])/2
        phav = (lons[ilon+1,ilat]+lons[ilon,ilat])/2
        
        # check the result
        check= check + dth*dph*np.sin(thav) 
        
        # check the green's functions
        dI_dBr_check = dI_dBr_check + dIdBr[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])*dth*dph*np.sin(thav) 
        dI_dBr_GR1983 = dI_dBr_GR1983 + GIs[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])*dth*dph*np.sin(thav)
        dI_dBr_JC1993 = dI_dBr_JC1993 + GIc[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])*dth*dph*np.sin(thav)
        
        dBra_dBr_GR1983 = dBra_dBr_GR1983 + Gr[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])*dth*dph*np.sin(thav)
        dBta_dBr_GR1983 = dBta_dBr_GR1983 + Gth[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])*dth*dph*np.sin(thav)
        dBpa_dBr_GR1983 = dBpa_dBr_GR1983 + Gph[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])*dth*dph*np.sin(thav)

        dX_dBr_JC1993 = dX_dBr_JC1993 + GX[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])*dth*dph*np.sin(thav)
        dY_dBr_JC1993 = dY_dBr_JC1993 + GY[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])*dth*dph*np.sin(thav)
        dZ_dBr_JC1993 = dZ_dBr_JC1993 + GZ[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])*dth*dph*np.sin(thav)
        
        # check the actual time diff
        dIdt_tot = dIdt_tot + dBrdt[ilon,ilat]*dIdBr[ilon,ilat]*dth*dph*np.sin(thav)
        dIdt_surf_tot = dIdt_surf_tot + dBradt[ilon,ilat]*dIdBra[ilon,ilat]*dth*dph*np.sin(thav)
        
        # reproduce some gauss coeffs (accurate to 2nd significant figure, not great)
        g10_check = g10_check +  (1/2)*(r_c/r_a)**3*(3/(4*np.pi)) * SH_library.SchmidtSH(1,0,thav,phav,'c')  * Br_c[ilon,ilat]  *dth*dph*np.sin(thav)
        h11_check = h11_check +  (1/2)*(r_c/r_a)**3*(3/(4*np.pi)) * SH_library.SchmidtSH(1,1,thav,phav,'s')  * Br_c[ilon,ilat]  *dth*dph*np.sin(thav)
        
        norm = 4*np.pi/(2*1+1)
        dgdBr_map =  (1/(1+1)) *  (r_c/r_a)**(1+2) * SH_library.SchmidtSH(1,0,thav,phav,'c') / norm
        g10dot_check = g10dot_check + dgdBr_map * (Br_c[ilon,ilat] - Brm_c[ilon,ilat])/(MFcoeffs[-1,0]-MFcoeffs[-2,0]) *dth*dph*np.sin(thav)
'''

# assume uniform phi grid
dph = abs(lons_lin[1]-lons_lin[0])
for ilon in range(len(lons_lin)-1):
    for ilat in range(len(theta_lin)):
        dI_dBr_check = dI_dBr_check + dIdBr[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])* dph * weights[ilat]
        dI_dBr_GR1983 = dI_dBr_GR1983 + GIs[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])* dph * weights[ilat]
        dI_dBr_JC1993 = dI_dBr_JC1993 + GIc[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])* dph * weights[ilat]
       
        dBra_dBr_GR1983 = dBra_dBr_GR1983 + Gr[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])* dph * weights[ilat]
        dBta_dBr_GR1983 = dBta_dBr_GR1983 + Gth[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])* dph * weights[ilat]
        dBpa_dBr_GR1983 = dBpa_dBr_GR1983 + Gph[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])* dph * weights[ilat]

        dX_dBr_JC1993 = dX_dBr_JC1993 + GX[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])* dph * weights[ilat]
        dY_dBr_JC1993 = dY_dBr_JC1993 + GY[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])* dph * weights[ilat]
        dZ_dBr_JC1993 = dZ_dBr_JC1993 + GZ[ilon,ilat]*(Br_c_plus[ilon,ilat]-Br_c[ilon,ilat])* dph * weights[ilat]
        
        # check the actual time diff
        dIdt_tot = dIdt_tot + dBrdt[ilon,ilat]*dIdBr[ilon,ilat]* dph * weights[ilat]
        dIdt_surf_tot = dIdt_surf_tot + dBradt[ilon,ilat]*dIdBra[ilon,ilat]* dph * weights[ilat]
        
        # reproduce some gauss coeffs 
        g10_check = g10_check +  (1/2)*(r_c/r_a)**3*(3/(4*np.pi)) * SH_library.SchmidtSH(1,0,theta_lin[ilat],lons[ilon,ilat],'c')  * Br_c[ilon,ilat] * dph * weights[ilat]
        h11_check = h11_check +  (1/2)*(r_c/r_a)**3*(3/(4*np.pi)) * SH_library.SchmidtSH(1,1,theta_lin[ilat],lons[ilon,ilat],'s')  * Br_c[ilon,ilat] * dph * weights[ilat]
 

print('****************************************') 
print('Let''s check these green''s functions')   
print('The actual inclination change is:')   
print('Inlc_plus - Incl = '+str(Incl_SUL_plus[0] -Incl_SUL[0])+' deg')   
print('To be compared with:')    
print('')     
print('1. The green''s function calculation:')         
print('\int dIdBr * (B_c_plus - B_c) dS = '+str(dI_dBr_check)+ ' deg')         
print('where dIdBr = dIdg * dgdBr_map')      
print('')
print('2. Gubbins + Roberts, 1983 Green''s functions calculation:')
print('\int GIs * (B_c_plus - B_c) dS = '+str(dI_dBr_GR1983)+ ' deg')     
print('In this specific case')
print(' \int Gr * (B_c_plus - B_c) dS =' + str(dBra_dBr_GR1983) + ' (vs '+str(Br_SUL_plus[0] - Br_SUL[0])+')')
print(' \int Gth * (B_c_plus - B_c) dS =' + str(dBta_dBr_GR1983) + ' (vs '+str(Bt_SUL_plus[0] - Bt_SUL[0])+')')
print(' \int Gph * (B_c_plus - B_c) dS =' + str(dBpa_dBr_GR1983) + ' (vs '+str(Bp_SUL_plus[0] - Bp_SUL[0])+')')
print('')
print('3. Johnson + Constable, 1993 Green''s functions calculation:')
print('\int GIc * (B_c_plus - B_c) dS = '+str(dI_dBr_JC1993)+ ' deg')         
print('In this specific case')
print(' \int GX * (B_c_plus - B_c) dS =' + str(dX_dBr_JC1993) + ' (vs '+str(-Bt_SUL_plus[0] + Bt_SUL[0])+')')
print(' \int GY * (B_c_plus - B_c) dS =' + str(dY_dBr_JC1993) + ' (vs '+str(Bp_SUL_plus[0] - Bp_SUL[0])+')')
print(' \int GZ * (B_c_plus - B_c) dS =' + str(dZ_dBr_JC1993) + ' (vs '+str(-Br_SUL_plus[0] + Br_SUL[0])+')')
print('')
print('')
print('****************************************')            
print('Let''s check the 2020-2015 difference')
print('The actual inclination change is:')   
print('(Inlc(2020) - Incl(2015))/5 = '+str((Incl_SUL -Inclm_SUL)/5)+' deg/yr')   
print('linearised variation = '+str(dIt_lin)+' deg/yr')   
print('\int dIdBc dBcdt dS = '+str(dIdt_tot)+' deg/yr')
print('\int dIdBe dBedt dS = '+str(dIdt_surf_tot)+' deg/yr')
print('\sum dI/dbeta * dbeta/dt = '+str(np.nansum(dIdtmatrix))+' deg/yr')
print('****************************************')            


#####################
# plots
######################

fig = plt.figure(figsize=(10,6))
Zs = GIs
Zmax = np.nanmax(abs(Zs))
Zmin = -Zmax
#Zmax = 1
#Zmin = -1
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               ,cmap='PiYG' 
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )

plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)

clb = plt.colorbar(cf
                   #,ticks=[-90,-45,0,45,90]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   ,format='%.1e'
                   )
clb.set_label(r'$G_I$  [deg/nT]', fontsize = 20)
#clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/GI_Hammer2019.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()

fig = plt.figure(figsize=(10,6))
Zs = GDs
Zmax = np.nanmax(abs(Zs))
Zmin = -Zmax
#Zmax = 1
#Zmin = -1
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               ,cmap='PiYG' 
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )

plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)


clb = plt.colorbar(cf
                   #,ticks=[-90,-45,0,45,90]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   ,format='%.1e'
                   )
clb.set_label(r'$G_D$  [deg/nT]', fontsize = 20)
#clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/GD_Hammer2019.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()

fig = plt.figure(figsize=(10,6))
Zs = GIc
Zmax = np.nanmax(abs(Zs))
Zmin = -Zmax
#Zmax = 1
#Zmin = -1
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               #,cmap='PiYG' 
               ,cmap='PuOr_r'
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )

plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)

clb = plt.colorbar(cf
                   #,ticks=[-90,-45,0,45,90]
                   ,ticks=[-1e-4, 0, 1e-4]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   ,format='%.1e'
                   )

plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)

clb.set_label(r'$\partial I / \partial B_c$ (analytical solution)  [deg/nT]', fontsize = 20)
#clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/GI_Johnson1997.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()

fig = plt.figure(figsize=(10,6))
Zs = GDc
Zmax = np.nanmax(abs(Zs))
Zmin = -Zmax
#Zmax = 1
#Zmin = -1
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               ,cmap='PiYG' 
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )

plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)


clb = plt.colorbar(cf
                   #,ticks=[-90,-45,0,45,90]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   ,format='%.1e'
                   )
clb.set_label(r'$G_D$  [deg/nT]', fontsize = 20)
#clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/GD_Johnson1997.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()


fig = plt.figure(figsize=(10,6))
Zs = GFc
Zmax = np.nanmax(abs(Zs))
Zmin = -Zmax
#Zmax = 1
#Zmin = -1
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               ,cmap='PiYG' 
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )

plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)


clb = plt.colorbar(cf
                   #,ticks=[-90,-45,0,45,90]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   ,format='%.1e'
                   )
clb.set_label(r'$G_F$  [nT/nT]', fontsize = 20)
#clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/GF_Johnson1997.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()

fig = plt.figure(figsize=(10,6))
Zs = GFc
Zmax = np.nanmax(abs(Zs))
Zmin = -Zmax
#Zmax = 1
#Zmin = -1
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global() # important (apparently)
ax.coastlines()
ax.gridlines()    

cf = plt.contourf(lons*180/np.pi,lats*180/np.pi,Zs
               ,cmap='PiYG' 
               ,levels=np.linspace(Zmin,Zmax,41)
               ,transform=ccrs.PlateCarree()
               )

plt.scatter(lon_SUL,90-colat_SUL,s=300,edgecolor='k',color='yellow',alpha=1,marker='*'
             , zorder=100)


clb = plt.colorbar(cf
                   #,ticks=[-90,-45,0,45,90]
                   ,fraction=0.05
                   ,orientation='horizontal'
                   ,aspect = 40
                   ,pad = 0.05
                   ,format='%.1e'
                   )
clb.set_label(r'$G_F$  [nT/nT]', fontsize = 20)
#clb.ax.set_xticklabels(['-90','-45','0','45','90'])
clb.ax.tick_params(labelsize=20)
#plt.title('AACGM latitudes',fontsize=20)
plt.savefig(folder_quantity + '/GF_Johnson1997.pdf',bbox_inches='tight',pad_inches=0.1)
plt.show()
