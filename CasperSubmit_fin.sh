#!/bin/bash -l
### Job Name
#PBS -N hail_testing
### Charging account
#PBS -A P66770001
### Request one chunk of resources with 1 CPU and 10 GB of memory
#PBS -l select=1:ncpus=1:mem=100GB
### Allow job to run up to 24 hrs
#PBS -l walltime=4:00:00
### Route the job to the casper queue
#PBS -q casper
### Join output and error streams into single file
#PBS -j oe


# module load python/2.7.14
module unload python
module load conda
conda activate npl
./Synthetic_Hail_Model_modified-final2.py 4 4 26 2022


