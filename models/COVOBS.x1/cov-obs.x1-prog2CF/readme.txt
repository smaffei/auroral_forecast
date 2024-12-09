# NG, 2013/04/23

- the program cov-obs-coefs.f produces times series of the Gauss coefficients and predictions at a given location for the COV-OBS.x1 ensemble of models

- for details about the COV-OBS.x1 field models see 

Gillet, N., Jault, D., Finlay, C. C., & Olsen, N. (2013). Stochastic modelling of the Earth's magnetic field: inversion for covariances over the observatory era. Geochemistry, Geophysics, Geosystems., doi:10.1002/ggge.20041

Gillet, N., Barrois, O., & Finlay, C. C. (2015). Stochastic forecasting of the geomagnetic field from the COV-OBS. x1 geomagnetic field model, and candidate models for IGRF-12. Earth, Planets and Space, 67(1), 71.

- the program offers 4 possibilities (change IFUNC in param.dat)

c	[1] evaluate the MF and SV Gauss coefficients (int and ext) for the average model,
c		from ts to te with fsamp sampling rate  
c	- outputs : files mf_int_coefs.dat, sv_int_coefs.dat, mf_ext_coefs.dat, sv_ext_coefs.dat
c	- storage per line: epoch, g10, g11, h11, g20, g21, h21, g22, h22, g30...
c
c	[2] evaluate the (X,Y,Z) prediction (MF and SV) at a given location
c		at the Earth's surface
c	- careful ! no constraint on the observatory bias !
c	- altitude, colatitude and longitude of the location are entered from the keyboard		
c	- outputs : file svpred.dat
c	- storage per line: epoch, Xint, Yint, Zint, Xext, Yext, Zext, Xind, Yind, Zind 
c
c	[3] evaluate an ensemble of NREAL (max provided=100) MF and SV Gauss coefficients (int and ext) 
c		generated from the posterior covariance matrix,	from ts to te with fsamp sampling rate  
c	- outputs : NREAL files real***mf_int_coefs.dat, real***sv_int_coefs.dat, real***mf_ext_coefs.dat, real***sv_ext_coefs.dat
c	- same storage as in [1]
c
c	[4] evaluate an ensemble of NREAL (max provided=100) (X,Y,Z) prediction (MF and SV) 
c		at a given location at the Earth's surface
c	- careful ! no constraint on the observatory bias !
c	- outputs : NREAL files real***svpred.dat
c	- same storage as in [2]

- "int", "ext" and "ind" refer to respectively internal, external and induced ;
"sv" and "mf" refer to respectively secular variation and main field

- the ensemble of field models have been obtained from the posterior covariance matrix, using its choleski decomposition as detailed in the g-cubed paper

- in the file param.dat, you can choose the period and sampling rate you wish for the outputs, as well as the number of realisations you wish for [3] and [4]. 

- to compile and run the program : from a terminal 

gfortran cov-obs-coefs.f
./a.out

- !! notice that there is a typo in figure 8 of the 2013 paper : the legend and ylabel of the plot should refer to q10 (and not q10 in dipole coordinates) !!
