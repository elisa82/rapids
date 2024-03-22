#!/bin/bash
module load python/3.8.12--gcc--8.4.1
module load netcdf-c/4.8.1--gcc--10.2.0

source rre_crs_env/bin/activate

python3 -m venv rre_crs_env
pip3 install basemap
pip3 install matplotlib
pip3 install numpy
pip3 install scipy
pip3 install shapely
pip3 install pyproj
pip3 install obspy
