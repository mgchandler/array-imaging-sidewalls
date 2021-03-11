#!/bin/bash
#
#
#PBS -l select=1:ncpus=1:mem=3500mb
#PBS -l walltime=20:00:00
#PBS -J 1-100

# Load modules
module load apps/matlab/r2019a

# Change to working directory. This should be the same location as the control Matlab script.
cd ${PBS_O_WORKDIR}

# Get useful information

echo JOB ID: ${PBS_JOBID}
echo ARRAY ID: ${PBS_ARRAY_INDEX}
echo Working Directory: $(pwd)
echo Start Time: $(date)
echo Run location:
hostname

echo Sidewall SDH Scan

# Create and define variables for use in the job.
declare -i Job Num Views
Job=${PBS_ARRAY_INDEX}
Num=100

# Execute code
echo matlab -nodisplay -nodesktop -singleCompThread -nosplash -r "WORKING_sdh_sidewall_scan(${Job}, ${Num});quit;"
matlab -nodisplay -nodesktop -singleCompThread -nosplash -r "WORKING_sdh_sidewall_scan(${Job}, ${Num});quit;"

echo End Time: $(date)