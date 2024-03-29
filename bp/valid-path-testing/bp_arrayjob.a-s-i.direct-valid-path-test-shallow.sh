#!/bin/bash
#
#
#PBS -l select=1:ncpus=1:mem=4000mb
#PBS -l walltime=03:59:59

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

echo Valid Path Testing - Direct SDH Shallow

# Execute code
echo matlab -nodisplay -nodesktop -singleCompThread -nosplash -r "WORKING_valid_path_direct_sdh_shallow();quit;"
matlab -nodisplay -nodesktop -singleCompThread -nosplash -r "WORKING_valid_path_direct_sdh_shallow();quit;"

echo End Time: $(date)