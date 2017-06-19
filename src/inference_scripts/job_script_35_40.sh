#!/bin/bash
# Job name:
#SBATCH --job-name=savio_test
#
# Partition:
#SBATCH --partition=savio2
#
# QoS:
#SBATCH --qos=savio_normal
#
# Account:
#SBATCH --account=fc_mhmm
#
# Request one node:
#SBATCH --nodes=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Number of Processors per Node:
#SBATCH --ntasks-per-node=24
#
# Wall clock limit:
#SBATCH --time=36:00:00
#
## Command(s) to run:
module load matlab/R2015a
# Make a temporary scratch directory for storing job
# and task information, to coordinate parallelization.
# This directory can then be referenced by assigning it to
# a 'parcluster.JobStorageLocation' property in your script.
mkdir -p /global/scratch/$USER/$SLURM_JOB_ID
matlab -nodisplay -nosplash -nodesktop < inf_multi_ap_35_40.m