#!/bin/bash
#
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=updateR
#SBATCH --output=updateR.out
sudo amazon-linux-extras install R4
sudo Rscript pkg_installs.R