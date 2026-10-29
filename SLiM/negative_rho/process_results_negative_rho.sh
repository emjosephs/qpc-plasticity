#!/bin/bash --login
#SBATCH --job-name=process_qpc_results_negative_rho
#SBATCH -o /mnt/research/josephslab/Husain/cond_Qpc/negative_rho/logs/process_sims.o
#SBATCH -e /mnt/research/josephslab/Husain/cond_Qpc/negative_rho/logs/process_sims.e
#SBATCH --nodes=1-1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=500G
#SBATCH -t 07:59:59
#SBATCH --mail-user=middlestauthor@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --account=josephsnodes

cd /mnt/research/josephslab/Husain/cond_Qpc/negative_rho/
module purge
module load R/4.4.1-gfbf-2023b

Rscript src/calc_qpc_negative_rho.r '/mnt/research/josephslab/Husain/cond_Qpc/negative_rho/sim_output/' 2610210001 2610210600
