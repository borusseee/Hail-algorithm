#!/bin/bash -l
### Job Name
#PBS -N ERA5_hail
### Charging account
#PBS -A P66770001
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=128GB
### Allow job to run up to 40 minutes
#PBS -l walltime=24:00:00
### Route the job to the casper queue
#PBS -q casper
### Join output and error streams into single file
#PBS -j oe


# to check process:  qstat -u $USER

# module load python/2.7.14
ncar_pylib
module load conda
conda activate npl
./ERA5_hail_proxies_extraction_newSRH03.py MM YYYY


