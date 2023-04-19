# SMBHB-merger-delay
This is a combined Python and Fortran code for GW Population Analysis using hierarchical Bayesian inference to infer the delay time distribution of SMBH binary mergers. 

Author: Yun Fang & Huan Yang

This data analysis make use of the information of the marginalized joint distribution of binary black hole mass and redshift in GW events. The delay time distribution is estimated by combining the merger rate of SMBH binaries and galaxies, together with the fundamental relationship between SMBHs and galaxies. 

The uncertainties from the galaxy merger rate and fundamental relationship have been considered in the main code main_gaussian.py. The evaluation of likelihood is sped up with Fortran subroutine and imported by pyhton. 

Compile the Fortran code with the command:<br />
f2py -c -m bhmergerrate MergerRateSMBH.f90<br />
f2py -c -m bhmergernumber MergerNumberSMBH.f90

The interpolation data for the merger rate in the Power-law delay model can be downloaded from the following link:<br />
https://www.dropbox.com/s/7f1dyp1bhfhhn3p/powerlaw_delay.zip?dl=0

Citation of the reference will be appreciated when you use the code/data:<br />
Yun Fang & Huan Yang<br />
``Probing the Delay Time of Supermassive Black Hole Binary Mergers With Gravitational Waves"<br />
arXiv: 2209.14509<br />
