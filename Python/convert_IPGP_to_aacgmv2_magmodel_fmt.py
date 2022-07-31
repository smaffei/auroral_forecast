#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:42:03 2021

@author: stefanomaffei
convert IPGP model in a fomat that aacgmv2 accepts

after many trial-and-errors, I think we need to copy-paste the magmodel file 
in the aacgmv2 repository and modify it with our coefficents. 
aacgmv2 somehow reads things in the header that I am not aware of
"""



import numpy as np

folder = './'

MF_file = folder +'gauss_MF_2015-2115'
SV_file = folder + 'gauss_SV_2015-2115'
 
out_file = folder + 'ipgpMF4aacgmv2_from1590_axial_dipole.txt'

# read coeffs
MFcoeffs = np.loadtxt(MF_file)
SVcoeffs = np.loadtxt(SV_file)

# to convince aacgmv2 to accept these data we need to time them so that the series starts in 1590 and does not go further than 2020
t0 = 1590.0 # start of the igrf13 model file
tf = 2020.0 # end of the igrf13 model file

# new year array:
years = MFcoeffs[:,0]-MFcoeffs[0,0]+t0



f = open(out_file,'w')
# header
f.write('# 2015 IPGP data assimilation forecast\n')
f.write('# units of nT \n')
# 1st row
# copy it from the magmodel file
f.write('c/s d o  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1  GUFM1   IGRF   IGRF   IGRF   IGRF   IGRF   IGRF   IGRF   IGRF   IGRF   DGRF   DGRF   DGRF   DGRF   DGRF   DGRF   DGRF   DGRF   DGRF   DGRF   DGRF   DGRF      DGRF      DGRF     DGRF     IGRF      SV\n')

# 2nd row, years
# copy it from the magmodel file
f.write('g/h n m 1590.0 1595.0 1600.0 1605.0 1610.0 1615.0 1620.0 1625.0 1630.0 1635.0 1640.0 1645.0 1650.0 1655.0 1660.0 1665.0 1670.0 1675.0 1680.0 1685.0 1690.0 1695.0 1700.0 1705.0 1710.0 1715.0 1720.0 1725.0 1730.0 1735.0 1740.0 1745.0 1750.0 1755.0 1760.0 1765.0 1770.0 1775.0 1780.0 1785.0 1790.0 1795.0 1800.0 1805.0 1810.0 1815.0 1820.0 1825.0 1830.0 1835.0 1840.0 1845.0 1850.0 1855.0 1860.0 1865.0 1870.0 1875.0 1880.0 1885.0 1890.0 1895.0 1900.0 1905.0 1910.0 1915.0 1920.0 1925.0 1930.0 1935.0 1940.0 1945.0 1950.0 1955.0 1960.0 1965.0 1970.0 1975.0 1980.0 1985.0 1990.0 1995.0 2000.0    2005.0    2010.0   2015.0   2020.0 2020-25\n')

# the magmodel file has this format:
Nr = int(1+(2020-1590)/5) # number of rows   (excluding the SV coeffs) 
Nc = 196         # number of coefficients + one column for the year (leave it empty)
coeffs = np.zeros((Nr,Nc))



# write the coefficients
Nd = MFcoeffs.shape[1] # number of columns in MFcoeffs
l=1
m=0
cs='g'
for j in range(1,Nd): # cycle over the columns degree and order of SH
    #print(l,m,cs)
    f.write( cs+'  '+str(l)+'  '+str(m)+'  ' )
    for i in range(0,int(1+(2020-1590)),5): # print every 5 years, from 1590 to 2020
        if i < MFcoeffs.shape[0]:
            f.write(str(MFcoeffs[i,j]) + '  ')       
        else:
            f.write(str(0) + '  ')
    f.write(str(SVcoeffs[-1,j]) + '\n')    # keep the last value of the forecast
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