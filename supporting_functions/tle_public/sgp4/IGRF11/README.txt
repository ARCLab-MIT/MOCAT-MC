To execute [B]=igrf11syn(fyears,alt,nlat,elong)
first execute GetIGRF11_Coefficients(1)

GetIGRF11_Coefficients(1) loads the global variable gh with 
spherical harmonic coefficients
Initially the coefficients are defined in the spread sheet  igrf11coeffs.xls 
If the *.mat file GHcoefficients.mat is not found it is written and used for
subsequent calls.

NOTE:  fyears=date in fractional years
       alt   = altitude in km
       nlat  = latitude positive North (deg)
       elong = longitude positive East (deg)

       B     =[B_North; B_East; B_up]   
 