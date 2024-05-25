The Northern Adria Crust (NAC) model

Authors:
A. Magrin, G. Rossi
OGS, Centro Ricerche Sismologiche, Istituto Nazionale di Oceanografia e di Geofisica Sperimentale, Via Treviso 55, 33100 Udine, Italy
OGS, Centro Ricerche Sismologiche, Istituto Nazionale di Oceanografia e di Geofisica Sperimentale, Borgo Grotta Gigante, 42/c, 34010 Trieste, Italy

Contact:
Andrea Magrin @ OGS (Udine - Italy) amagrin@inogs.it

The dataset provided here is the crustal model presented in Magrin and Rossi, Front. Earth Sci., 2020 (doi: 10.3389/feart.2020.00089). 

If you use this dataset, please cite the following paper:

Magrin, A., Rossi, G. (2020). Deriving a new crustal model of Northern Adria: the Northern Adria Crust (NAC) model. Front. Earth Sci, 10.3389/feart.2020.00089


NAC/NAC1 : version of the model with mo1 (Moho defined by a continuous surface)
NAC/NAC1/GRD : files in NETCDF4 format
               nac_interfaces_mo1.grd	-> depth of the interfaces with their error
               nac_parameters_mo1.grd	-> geophysical parameters of the crust with their error
NAC/NAC1/XYZ
               nac_interfaces_mo1.xyz	-> depth of the interfaces with their error
               nac_parameters_mo1.xyz	-> geophysical parameters of the crust with their error

NAC/NAC2 : version of the model with mo2 (segmented Moho)
NAC/NAC2/moho_fragments.txt : boundaries between AD, EU and PA
NAC/NAC2/GRD : files in NETCDF4 format
               nac_interfaces_mo2.grd	-> depth of the interfaces with their error
               nac_parameters_mo2.grd	-> geophysical parameters of the crust with their error
NAC/NAC2/XYZ
               nac_interfaces_mo2.xyz	-> depth of the interfaces with their error
               nac_parameters_mo2.xyz	-> geophysical parameters of the crust with their error


The Cartesian coordinates of the model are referred to the origin point 
(at 10.2° E and 44.6° N) and are obtained from geographical coordinates
using the UTM (Universal Transverse Mercator) projection (zone 33N). 
z increases with depth below sea level.

Interfaces are:
- T: topography
- BS: the bottom of the sedimentary layer
- MO1: Moho defined by a continuous surface
- MO2: segmented Moho

The geophysical parameters of the crust are:
- Vp: velocity of P-wave (km/s)
- Vs: velocity of S-wave (km/s)
- rho: density (10^3 kg/m^3)
- mu: shear modulus (10^9 N/m^2)
- young: young modulus (10^9 N/m^2)

The errors of interfaces and parameters are the total errors (see paper for a complete description).
