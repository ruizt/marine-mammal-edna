#!/bin/bash
#
#SBATCH --array=0-4
#SBATCH --job-name=packages
#SBATCH --output=packages_%a.out
#SBATCH --mail-user truiz01@calpoly.edu
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
sudo Rscript --vanilla pkg_installs.R