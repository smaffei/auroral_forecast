#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:42:03 2021

@author: stefanomaffei
modify the magmodel file
"""



import numpy as np

folder = '/Users/stefanomaffei/Documents/geomag/models/magmodel/'

magfile = folder +'magmodel_1590-2020.txt'
 
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

'''
# modify to have only the dipole
coeffs_MF[:,4:] = 0
coeffs_SV[3:]=0
'''

# modify to have only the axial dipole
coeffs_MF[:,2:] = 0
coeffs_SV[1:]=0

coeffs_MF[:,2:4] = 0.1


dipole_file = folder + 'magmodel_1590-2020_axial_dipole_smallH.txt'

f = open(dipole_file,"w")
f.write(header1)
f.write(header2)
f.write(header3)
f.write(header4)

# write the coefficients
Nd = coeffs_MF.shape[1] # number of columns in MFcoeffs
l=1
m=0
cs='g'
for j in range(1,Nd): # cycle over the columns degree and order of SH
    #print(l,m,cs)
    f.write( cs+'  '+str(l)+'  '+str(m)+'  ' )
    for i in range(0,coeffs_MF.shape[0]): # print every 5 years, from 1590 to 2020
        if i < coeffs_MF.shape[0]:
            f.write(str(coeffs_MF[i,j]) + '  ')       
        else:
            f.write(str(0) + '  ')
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

f.close()
