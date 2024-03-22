#!/bin/bash
#SBATCH --job-name=PostProcessing
#SBATCH -t 1:00:00
#SBATCH --account=IscrC_URASS-2
#SBATCH --partition=g100_usr_dbg
  
module load python/3.8.12--gcc--8.4.1
module load netcdf-c/4.8.1--gcc--10.2.0

source rre_crs_env/bin/activate
python -m rapids prova.ini speed --post
