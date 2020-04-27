#!/bin/sh

IN=$1

echo "#!/bin/bash" > $IN.sh
echo "#SBATCH --job-name=corncob_sim_$IN" >> $IN.sh 
echo "#SBATCH --priority=normal" >> $IN.sh 
echo "#SBATCH --ntasks-per-core=1" >> $IN.sh 
echo "#SBATCH --ntasks-per-node=1" >> $IN.sh 
echo "#SBATCH --cpus-per-task=1" >> $IN.sh 

echo "module purge" >>$IN.sh
echo "module load R/3.5.2" >>$IN.sh

echo "
Rscript \
./eval_functions_call_corncob.R \
./data/Stool_16S_wMG_split/sim$IN.RData \
./data/Stool_16S_wMG_split/corncob/evals_sim$IN.RDS" >> $IN.sh
