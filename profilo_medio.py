#!/usr/bin/env python3

import argparse
import numpy as np
import netCDF4 as nc
import sys
import os


# Extracting NAC parameters from the netcdf files
def import_nac_param(nac_path):
    nac_files = [ os.path.join(nac_path, 'NAC1', 'GRD', 'nac_parameters_mo1.grd'),
              os.path.join(nac_path, 'NAC1', 'GRD', 'nac_interfaces_mo1.grd') ]

    nac_param_nc = nc.Dataset( nac_files[0] )
    nac_inter_nc = nc.Dataset( nac_files[1] )

    nac_param = {}
    nac_param['x'] = nac_param_nc.variables['x'][:]
    nac_param['y'] = nac_param_nc.variables['y'][:]
    nac_param['z'] = nac_param_nc.variables['z'][:] 
    nac_param['lon'] = nac_param_nc.variables['Longitude'][:]
    nac_param['lat'] = nac_param_nc.variables['Latitude'][:]
    nac_param['rho'] = nac_param_nc.variables['Rho'][:]
    nac_param['vp'] = nac_param_nc.variables['Vp'][:]
    nac_param['vs'] = nac_param_nc.variables['Vs'][:]

    return nac_param

################################################################################################################################
# main
################################################################################################################################

# default value of nac_path
nac_path = '/home/magrin/NAC/ModelloFinale/NAC/'

# define parser
parser = argparse.ArgumentParser(description='Compute mean vertical profile from NAC model')
parser.add_argument('--x', help='define min and max cartesian coordinate x',nargs=2)
parser.add_argument('--y', help='define min and max cartesian coordinate y',nargs=2)
parser.add_argument('--z', help='define min and max cartesian coordinate z',nargs=2)
parser.add_argument('--out', help='define name of output file',nargs=1)
parser.add_argument('--nac_path', help='define path of NAC model files',nargs=1)
parser.add_argument('--lon', help='define min and max longitude',nargs=2)
parser.add_argument('--lat', help='define min and max latitude',nargs=2)

args = parser.parse_args()

cartesian = False
geographic = False
if  args.x or args.y:
   cartesian = True
if  args.lon or args.lat:
   geographic = True

if  cartesian and geographic:
    print('ERROR - only cartesian or only geographic coordinates should be used')
    sys.exit(1)
if  not geographic:
    cartesian = True

if  args.nac_path:
    nac_path = args.nac_path[0]

# import NAC paramaters
nac_param = import_nac_param(nac_path)

if  cartesian:
    minx = np.min(nac_param['x'])
    maxx = np.max(nac_param['x'])
    miny = np.min(nac_param['y'])
    maxy = np.max(nac_param['y'])
    if  args.x:
        minx, maxx = args.x
        minx = float(minx)
        maxx = float(maxx)
    if  args.y:
        miny, maxy = args.y
        miny = float(miny)
        maxy = float(maxy)
else:
    minlo = np.min(nac_param['lon'])
    maxlo = np.max(nac_param['lon'])
    minla = np.min(nac_param['lat'])
    maxla = np.max(nac_param['lat'])
    if  args.lon:
        minlo, maxlo = args.lon
        minlo = float(minlo)
        maxlo = float(maxlo)
    if  args.lat:
        minla, maxla = args.lat
        minla = float(minla)
        maxla = float(maxla)

# limits on vertical coordinate
minz = 0 # mean values above topography can't computed correctly using a simple mean
maxz = np.max(nac_param['z'])

if  args.z:
    minz, maxz = args.z
    minz = float(minz)
    maxz = float(maxz)

# output file
if  args.out:
    fout = args.out[0]
else: 
    fout = 'profile_1D.txt'

# compute mean profile in the region

params = ('vp','vs','rho')
type_profile1D = [('z','f')]
for p in params:
    type_profile1D.append((p,'f'))
    type_profile1D.append(('d'+p,'f'))

# define index for region
if cartesian:
    yy, xx = np.meshgrid(nac_param['y'],nac_param['x'],indexing='ij')
    ix = np.logical_and(np.logical_and(xx >= minx, xx <= maxx), np.logical_and(yy >= miny, yy <= maxy))
else:
    ix = np.logical_and(np.logical_and(nac_param['lon'] >= minlo, nac_param['lon'] <= maxlo), 
                        np.logical_and(nac_param['lat'] >= minla, nac_param['lat'] <= maxla))

iz = np.logical_and(nac_param['z']<=maxz,nac_param['z']>=minz)
nac_param['z'] = nac_param['z'][iz]
for p in params:
    if  p not in ['x','y','z']:
        nac_param[p] = nac_param[p][iz]

profile_1D = np.zeros(nac_param['z'].shape,
             dtype=type_profile1D)

for i,z in enumerate(nac_param['z']):
    profile_1D['z'][i] = z
    for p in params:
        profile_1D[p][i] = np.nanmean(nac_param[p][i,:,:][ix] )
        profile_1D['d'+p][i] = np.nanstd(nac_param[p][i,:,:][ix] )

# write output file
if cartesian:
    header = '1D profile for derived from NAC (x='+str(minx)+'-'+str(maxx)+', y='+str(miny)+'-'+str(maxy)+') (dVp, ... are standard deviations, not erros)\n'
else:
    header = '1D profile for derived from NAC (lon='+str(minlo)+'-'+str(maxlo)+', lat='+str(minla)+'-'+str(maxla)+') (dVp, ... are standard deviations, not erros)\n'
header = header + ' z(km)   Vp(km/s)  dVp(km/s)  Vs(km/s)  dVs(km/s)  rho(10^3 kg/m^3) drho(10^3 kg/m^3)'

np.savetxt(fout, np.column_stack(tuple([ profile_1D[i] for i in list(profile_1D.dtype.fields) ])),fmt='%9.2f',header=header)


