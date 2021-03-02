#!/bin/bash
#
#
#PBS -l select=1:ncpus=1:mem=3500mb
#PBS -l walltime=16:00:00

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

echo Backwall SDH High Res

# Execute code
echo matlab -nodisplay -nodesktop -singleCompThread -nosplash -r "WORKING_high_res_sdh();quit;"
matlab -nodisplay -nodesktop -singleCompThread -nosplash -r "WORKING_high_res_sdh();quit;"

echo End Time: $(date)