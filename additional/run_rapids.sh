#!/bin/bash

current_date=$(date +"%Y-%m-%d")
current_time=$(date +"%H:%M:%S")
DIR='/g100_scratch/userexternal/ezuccolo/'$current_date'_'$current_time
echo $DIR

module load python/3.8.12--gcc--8.4.1
module load netcdf-c/4.8.1--gcc--10.2.0
module load gmt/6.2.0--gcc--10.2.0

source rre_crs_env/bin/activate

python launch_scenarios.py $DIR

python -m rapids rts.ini speed --input

cd $DIR'/SPEED'

JOB_STRING=$( sbatch run_speed_cineca.slurm )
JOBID=$( echo $JOB_STRING | cut -d " " -f4)

cd $HOME

sbatch --dependency=afterok:$JOBID post_processing.slurm
