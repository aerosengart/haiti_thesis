#!/bin/bash

#SBATCH --job-name=mif_panel_global-search
#SBATCH --mail-user=aelr@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --account=stats_dept1
#SBATCH --partition=standard

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1

## wall time hours:minutes:seconds
#SBATCH --time=01-18:00:00

###   Load software modules

module load R/4.1.0
module list

####  Commands your job should run follow this line

echo "Running on $SLURM_JOB_NODELIST"
echo "Running in $(pwd)"

Rscript --vanilla model1/scripts/panel_mif_search.R;