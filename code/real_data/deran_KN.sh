#!/bin/bash
#SBATCH -n 1                                     # Number of cores
#SBATCH -t 0-24:00                               # Runtime in D-HH:MM
#SBATCH -p serial_requeue                        # Partition to submit to
#SBATCH --mem=10000M                             # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END                          # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=lbywhj@mail.ustc.edu.cn   # Email
#SBATCH --array=1-1
#SBATCH -o out/err_deran_KN
#SBATCH -J deran_KN

# LOAD SOFTWARE ENV ##
module load gcc/7.1.0-fasrc01 R/3.6.3-fasrc01
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
input=deran_KN.R
cd /n/home09/cdai/FDR/real_data/code

for att in 1 2 3 4 5; do
#
echo "${att}"
export att
#

R CMD BATCH $input out/$input.$SLURM_ARRAY_TASK_ID.out

sleep 1                                          # pause to be kind to the scheduler
done
