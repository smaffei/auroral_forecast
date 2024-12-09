
README
======
In this directory you will find a simple example of generating predictions
of the CHAOS field model (both internal and external field parts). This
example should be modified by users to fit the task required.

Inputs
------
CHAOS-6-x9.mat : Latest CHAOS model in .mat format
RC_1997-2019_April_v3.dat : RC index  (needed for hourly 
    time-dependence of q10 in SM co-ords; get latest version** if wish to
    use after April, 2019).

    ** The most recent version of the RC index is available at
    http://www.spacecenter.dk/files/magnetic-models/RC/current/

input.dat : Example input file,
    Each row defines the location and time of one datum to be predicted.
    The 4 columns are geoc. co-lat. (degs), geoc. long. (degs), alt. (km)
    and time (mjd2000)

Script
------
CHAOS_example.m : Reads CHAOS model, RC, and data locations/times
    from input.dat. Calculates predictions of internal and external parts
    of the CHAOS model and outputs predictions

Outputs
-------
example_output.dat : Output results as acii file.  Each row corresponds to
    an input point. Columns are, Time, Radius, Co-lat (deg), Long (deg),
    CHAOS model predictions (total, internal, external) for B_r,
    B_theta, B_phi (nT).

Needs in path/or same directory:  

design_SHA_sm.m, design_SHA_gsm.m, design_SHA.m, geo2sm.m,
GSM_SM_SHA.mat, mat_mul_mat.m, sunGEI.m, synth_values.m
synth_values_CHAOS_ext.m

Changelog
---------
March 26, 2019 :

(1) Corrected bug in synth_values_CHAOS_ext which raised an error when exactly
    6 input arguments were given. Dipole and GSM/SM transformation
    coefficients are now handled by the underlying functions: geo2sm and
    design_SHA_sm.

March 5, 2019 :

(1) New coefficient file GSM_SM_SHA_IGRF12_2015.mat. Changed the main script
    name to CHAOS_example.m. It return example_output.dat.

October 8, 2018 :

(1) The current version has been updated to include the latest CHAOS-6 model
    and RC-index file. Note that some functions** have been slightly changed to
    accept the IGRF dipole direction as input argument.
    This ensures that the coordinate transformations are based on the same
    IGRF (currently IGRF-12 epoch 2015). Similarly, the path to the
    coefficient file GSM_SM_SHA is now an optional input argument.

    **design_SHA_sm, design_SHA_gsm, geo2sm, synth_values_CHAOS_ext

(2) A bug was fixed in "geo2sm.m" which caused the coordinate
    transformation to return wrong results when exactly 3 points were given as
    input. This undesired behavior resulted from the fact that MATLAB computes
    the cross product over the first dimension where it finds 3 elements. It
    therefore mixed up the 3 points with the 3 vector components at each point.

Clemens Kloss

Contact
-------
If you have questions, please contact Chris Finlay at cfinlay@space.dtu.dk
or Clemens Kloss at ancklo@space.dtu.dk

