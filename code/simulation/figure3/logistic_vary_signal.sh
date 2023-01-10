#!/bin/bash
#SBATCH -n 1                                     # Number of cores
#SBATCH -t 0-24:00                               # Runtime in D-HH:MM
#SBATCH -p serial_requeue                       # Partition to submit to
#SBATCH --mem=10000M                             # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-type=END                          # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=   # Fill in with your own email
#SBATCH --array=1-50
#SBATCH -o out/err_vary_signal
#SBATCH -J signal

## LOAD SOFTWARE ENV ##
module load Anaconda3/2019.10
## I download the requried package in the FDR environment, change this.
source activate FDR
cd ~/code/simulation/figure3

for att in 8 9 10 11 12; do
#
echo "${att}"
export att
#
python logistic_vary_signal.py

sleep 1                                          # pause to be kind to the scheduler
done
