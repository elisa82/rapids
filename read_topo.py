#!/usr/bin/env python
from __future__ import print_function

import os
import sys
import fileinput
import string
import math

import cubit

print(sys.path)

#############################################################
# USER PARAMETERS

# topography file, data points per lon-increment
#inputFile = '/Users/elisa/Documents/Progetti/SimulazioniFriuli/run_ext_topo/SPEED/MESH/ptopo.mean.utm'
inputFile = '/Users/elisa/Documents/Progetti/rapids/DATA/topo/Friuli.utm'

# X coordinate in topography file repeats after line
#nstep = 529
nstep = 497

#############################################################

# converts xyz to utm
#os.system('./convert_lonlat2utm.pl Friuli.xyz 33 > Friuli.utm ')

print('#reading from file: ', inputFile)

cubit.cmd('reset')
cubit.cmd('echo off')
cubit.cmd('Graphics Pause')


# creates point vertices
print('#creating points...')
count = 0
for line in fileinput.input(inputFile):
  count = count + 1
  lineitems = line.split()
  x = lineitems[0]
  y = lineitems[1]
  z = lineitems[2]
  xyz = str(x) + ' ' + str(y) + ' ' + str(z)
  cubit.cmd('create vertex ' + xyz)
fileinput.close()
print('#done points: ' + str(count))
print('')
cubit.cmd('Display')

# creates smooth spline curves for surface
print('#creating curves...')
countcurves = 0
for i in range(1, count+1):
  if i > 1:
    if i % nstep == 0:
      countcurves = countcurves + 1
      cubit.cmd('create curve spline vertex ' + str(i-nstep+1) + ' to ' + str(i) + ' delete')
print('#done curves: '+str(countcurves))
print('')
cubit.cmd('Display')
cubit.cmd('pause')
print('')
# creates surface
print('#creating skin surface...')
cubit.cmd('create surface skin curve all')

cubit.cmd('Display')
cubit.cmd('pause')

# cleans up
cubit.cmd('merge all ')
cubit.cmd('delete vertex all')
cubit.cmd('delete curve all')

print('#done cleaning up')
cubit.cmd('Display')
cubit.cmd('echo on')

# saves and exports surface
# cubit file (uses ACIS file format)
#cubit.cmd('save as "/Users/elisa/Documents/Progetti/SimulazioniFriuli/run_ext_topo/SPEED/MESH/topo.cub" overwrite')
cubit.cmd('save as "/Users/elisa/Documents/Progetti/rapids/DATA/topo/Friuli.cub" overwrite')

print('#exporting done')
print('#finished')
cubit.cmd('exit')
