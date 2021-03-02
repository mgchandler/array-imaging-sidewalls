#!/bin/bash
#
#
#PBS -l select=1:ncpus=1:mem=2000mb
#PBS -l walltime=12:00:00
#PBS -J 1-90

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

# Create and define variables for use in the job.
declare -i Job Num Views
Job=${PBS_ARRAY_INDEX}
Num=90
Views=2

# Execute code
echo matlab -nodisplay -nodesktop -singleCompThread -nosplash -r "WORKING_crack_rotation_sens(${Job}, ${Num}, ${Views});quit;"
matlab -nodisplay -nodesktop -singleCompThread -nosplash -r "WORKING_crack_rotation_sens(${Job}, ${Num}, ${Views});quit;"

echo End Time: $(date)