#!/bin/sh
#SBATCH --job-name=basic_r_job
#SBATCH --cpus-per-task=1  
#SBATCH --nodes=1
#SBATCH --array=1-51
#SBATCH --ntasks-per-node=1
#SBATCH --output=job.%J.out # tell it to store the output console text to a file called job.<assigned job number>.out
#SBATCH --exclude=node018   # ignore the gpu node.. it dosnt have the R libraries installed. 

#below this line is where we can place our commands, in this case it will just simply output the task ID of the array
# the array=1-51 command tells slurm to launch 51 INDEEPENDENT JOBS.  
# 

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running $SLURM_NTASKS tasks."
echo "Current working directory is `pwd`"
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
/bin/hostname
Rscript USA_UR_Deaths_LancetID.R $SLURM_ARRAY_TASK_ID
#srun -l $SLURM_ARRAY_TASK_ID