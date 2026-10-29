#!/bin/bash --login
#SBATCH --job-name=slim_symmetric
#SBATCH -o /mnt/research/josephslab/Husain/cond_Qpc/symmetric/logs/qpc_plasticity_%a.o
#SBATCH -e /mnt/research/josephslab/Husain/cond_Qpc/symmetric/logs/qpc_plasticity_%a.e
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH -t 23:59:59
#SBATCH --mail-user=middlestauthor@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=josephsnodes
#SBATCH --account=josephsnodes

module purge
module load VCFtools/0.1.16-GCC-12.3.0
module load  Conda/3
conda activate myenv

array=$SLURM_ARRAY_TASK_ID
mypath="/mnt/research/josephslab/Husain/cond_Qpc/symmetric/"
cd ${mypath}

# Where to save output from SLiM simulations and processed trees
outpath="/mnt/research/josephslab/Husain/cond_Qpc/symmetric/sim_output"

# File containing parameters for each array task. Each line corresponds to a different set of 
# parameters for a SLiM simulation, and the line number corresponds to the array task ID. The
# array number starts at 2 to skip the header line in the params file. The parameters are:
# 1. metapopSide: the number of subpopulations on each side of the metapopulation
# 2. N: the population size of each subpopulation
# 3. optimum1: the phenotypic optimum in environment 1
# 4. optimum2: the phenotypic optimum in environment 2
# 5. env2prob: the probability of being in environment 2
# 6. rho: the genetic correlation between the trait in environment 1 and the trait in environment 2
# 7. sdW: the standard deviation of the fitness function
# 8. n_tot: the total population size of the metapopulation
# 9. seed: the set seed for the SLiM simulation
params="/mnt/research/josephslab/Husain/cond_Qpc/symmetric/src/params260817.txt"

# Extract parameters for the current array task from the params file
metapopSide=`awk -v var=$array 'NR==var' $params | awk '{print $1}'`
N=`awk -v var=$array 'NR==var' $params | awk '{print $2}'`
optimum1=`awk -v var=$array 'NR==var' $params | awk '{print $3}'`
optimum2=`awk -v var=$array 'NR==var' $params | awk '{print $4}'`
env2prob=`awk -v var=$array 'NR==var' $params | awk '{print $5}'`
rho=`awk -v var=$array 'NR==var' $params | awk '{print $6}'`
sdW=`awk -v var=$array 'NR==var' $params | awk '{print $7}'`
n_tot=`awk -v var=$array 'NR==var' $params | awk '{print $8}'`
seed=`awk -v var=$array 'NR==var' $params | awk '{print $9}'`

#run slim in background
slim -d path=\"${mypath}\" \
  -d metapopSide=${metapopSide} \
  -d N=${N} \
  -d optimum1=${optimum1} \
  -d optimum2=${optimum2} \
  -d env2prob=${env2prob} \
  -d rho=${rho} \
  -d seed=${seed} \
  -d sdW=${sdW} \
  /mnt/research/josephslab/Husain/cond_Qpc/symmetric/src/qpc_plasticity.slim \
  > ${outpath}/slimLogs/${seed}.o \
  2> ${outpath}/slimLogs/${seed}.e

python3 /mnt/research/josephslab/Husain/cond_Qpc/symmetric/src/tree_seq_pyslim.py \
  -s ${seed} \
  -r 1e-06 \
  -mu 1e-6 \
  -N ${n_tot} \
  -o ${outpath} \
  > ${outpath}/pytreeLogs/${seed}_pytree.o \
  2> ${outpath}/pytreeLogs/${seed}_pytree.e

vcftools --vcf ${outpath}/vcf/${seed}_mutated.vcf \
  --012 \
  --out ${outpath}/genoMat/${seed}_genoMat
