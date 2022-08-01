# auroral_forecast

## Overview

This repository contains estimates for the auroral zones boundaries for the years 2020 and 2070. This repository is a companion to the paper "Climatological Predictions of the Auroral Zone Locations Driven by Low- and High-Energy Space Weather Events" by Maffei et al.

These estimates are based on existing geomagnetic field models and forecasts. In particular the configuration of the geomagnetic field in 2020 is taken from the IGRF13 (Alken et al., 2021), while estimates for the geomagnetic field in 2070 are taken from three different sources: 

1) the data assimilation forecast from Sanchez et al., 2020 (MPG forecast);
2) a simple linear extrapolation of the IGRF13 model coefficients in time (IGRF13 forecast). Original dataset available at https://www.ncei.noaa.gov/products/international-geomagnetic-reference-field;
3) the single-epoch inversion forecast from Aubert, 2015 (IPGP forecast). Original dataset available at http://www.ipgp.fr/~aubert/tools.html#data.


## Content of the repository

### Dataset


The folder ```datasets``` contains all the raw data necessary to reproduce the study described in Maffei et al.

- ```datasets/Gorgon```: results from the magnetospheric MHD numerical simulations performed in Gorgon with the purpose of locating the poleward edge of auroral zones under different solar wind conditions (see Figure 4 of Maffei et al.);
-  ```datasets/models```: Gauss coefficients of the geomagnetic field forecasts used in this study. All datafiles are samples of the original datasets (see above) at 5 years intervals and in a format that is easily accepted by ```aacgmv2```. In all files each row indicates a specific Gauss coefficient's time series. In all models the time index begins in 1590 to be accepted by ```aacgmv2```. All files are generated from the original models via a Python script similar to ````convert_IPGP_to_aacgmv2_magmodel_fmt.py```, foun in ```Python/``` (see later):
    - ```datasets/models/ipgpMF4aacgmv2_from1590.txt```: IPGP forecast (for epochs 2015 onward);
    - ```datasets/models/magmodel_1590-2020.txt```: ```magmodel``` file, also provided in the ```aacgmv2``` package and reproduced here for convenience. It is a combination of the ```gufm1``` (Jackson et al., 2000) and ```IGRF13``` models and is valid for epochs 2015 onward. This is not a forecast;
    - ```datasets/models/magmodel_2020-2140.txt```: the IGRF13 forecast (valid from 2020 onward). This is a simple linear extrapolation of the IGRF13 coefficients evolution;
    - ```datasets/models/mpgMF4aacgmv2_from1590_L13.txt```: MPG averaged forecast (valid from 2020 onward);
    -  ```datasets/models/Ensemble```: MPG forecast ensemble members.
- ```datasets/simplemaps_worldcities_basicv1.73```: a list of world cities and their coordinates. Taken from https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwji1rqjnqX5AhUaHuwKHT6VCocQFnoECAkQAQ&url=https%3A%2F%2Fgitlab.huma-num.fr%2Fnlambert%2Fresources%2F-%2Ftree%2Fmaster%2Fdatasets%2Fsimplemaps_worldcities_basicv1.73&usg=AOvVaw34lTEr0D5PXeySE2pXxK4s


The folder ```Python``` contains the data on the auroral boundary (both for 2020 and 2070) and the Python scripts used to produce figure 6 of Maffei et al.:
- ```Python/coords2020_bisection``` and ```Python/coords2070_bisection``` contain files in which poleward and equatorward bounds are given as a function of longitude for both the auroral zones and danger zones. Bounds are computed for 2020 from the IGRF13 model and for 2070 from the forecasts mentioned above.
- ```Python/convert_IPGP_to_aacgmv2_magmodel_fmt.py``` is a Python scripts used to modify the original IPGP forecast coefficients file (see above) to prepare the file ```datasets/models/ipgpMF4aacgmv2_from1590.txt```. With small modifications, it can be used to convert any geomagnetic field model into a format easily accepted by ```aacgmv2```
- ```Python/auroral_forecast.py``` is the script used to reproduce figure 6 of Maffei et al. (also placed in ```Python/figures```). Given that the auroral bounds are already present in the repository (under ```Python/coords2020_bisection``` and ```Python/coords2070_bisection```) the script will not recalculate them. If the folders ```Python/coords2020_bisection``` and ```Python/coords2070_bisection``` are removed, the auroral bounds will be recalculated from the geomagnetic field models and forecast located in ```datasets/models/```.
- ```Python/aacgm_functions.py``` a collection of Python functions needed by ```Python/auroral_forecast.py```.

The above Python scripts have been written in Python 3.8



## References

Aubert, J.: Geomagnetic forecasts driven by thermal wind dynamics in Earth’s core, _Geophys. J. Int._ 203, 1738-1751, 2015, doi: 10.1093/gji/ggv394 

Alken, P., Thébault, E., Beggan, C.D. et al. International Geomagnetic Reference Field: the thirteenth generation. _Earth Planets Space_ 73, 49 (2021). https://doi.org/10.1186/s40623-020-01288-x

Jackson, A., Jonkers, A. R., & Walker, M. R. (2000). Four centuries of geomagnetic secular variation from historical records. Philosophical Transactions of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences, 358(1768), 957-990.

Maffei et al., Climatological Predictions of the Auroral Zone Locations Driven by Low- and High-Energy Space Weather Events. Submitted to _Scientific Reports_

Sanchez, S., Wicht, J. & Bärenzung, J. Predictions of the geomagnetic secular variation based on the ensemble sequential assimilation of geomagnetic field models by dynamo simulations. _Earth Planets Space_ 72, 157 (2020). https://doi.org/10.1186/s40623-020-01279-y



