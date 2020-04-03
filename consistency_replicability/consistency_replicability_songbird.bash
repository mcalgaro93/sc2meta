#!/bin/bash
#SBATCH --job-name=songbird   			# Job name
#SBATCH -p normal   					# priority
#SBATCH --mail-type=ALL         		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mcalgaro93@gmail.com # Where to send mail	
#SBATCH --nodes=1                   	# Use one node
#SBATCH --ntasks=1                  	# Run a single task
#SBATCH --mem-per-cpu=5gb           	# Memory per processor
#SBATCH --time=24:00:00             	# Time limit hrs:min:sec
#SBATCH --output=array_%A-%a.out    	# Standard output and error log
#SBATCH --array=1-20                 	# Array range
# This is an example script that combines array tasks with
# bash loops to process many short runs. Array jobs are convenient
# for running lots of tasks, but if each task is short, they
# quickly become inefficient, taking more time to schedule than
# they spend doing any work and bogging down the scheduler for
# all users.
#Set the number of runs that each SLURM task should do
PER_TASK=5
# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))
# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

module purge
module load Anaconda/2019.07
source activate songbird_env

# Run the loop of runs for this task.
run=$START_NUM
while [ $run -le $END_NUM ]; do 
	echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run; 
	
	biom convert -i './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset1_otutable.tsv' -o './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset1_otutable.biom' --to-hdf5 
	biom convert -i './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset2_otutable.tsv' -o './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset2_otutable.biom' --to-hdf5 
	biom convert -i './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset1_otutable.tsv' -o './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset1_otutable.biom' --to-hdf5 
	biom convert -i './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset2_otutable.tsv' -o './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset2_otutable.biom' --to-hdf5 
	biom convert -i './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset1_otutable.tsv' -o './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset1_otutable.biom' --to-hdf5 
	biom convert -i './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset2_otutable.tsv' -o './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset2_otutable.biom' --to-hdf5 
	
	mkdir './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset1'
	songbird multinomial --input-biom './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset1_otutable.biom' \
                         --metadata-file './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset1_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset1'
	
	mkdir './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset2'
	songbird multinomial --input-biom './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset2_otutable.biom' \
                         --metadata-file './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset2_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/16S_subgingival_supragingival_Comparison'$run'_Subset2'
	
	mkdir './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset1'
	songbird multinomial --input-biom './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset1_otutable.biom' \
                         --metadata-file './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset1_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset1'
	
	mkdir './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset2'
	songbird multinomial --input-biom './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset2_otutable.biom' \
                         --metadata-file './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset2_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/16S_gingiva_mucosa_Comparison'$run'_Subset2'
	
	mkdir './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset1'
	songbird multinomial --input-biom './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset1_otutable.biom' \
                         --metadata-file './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset1_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset1'
	
	mkdir './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset2'
	songbird multinomial --input-biom './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset2_otutable.biom' \
                         --metadata-file './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset2_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/16S_tonguedorsum_stool_Comparison'$run'_Subset2'
	
	############################################################################################################################
	
	biom convert -i './songbird/WMS_CRC_control_Comparison'$run'_Subset1_otutable.tsv' -o './songbird/WMS_CRC_control_Comparison'$run'_Subset1_otutable.biom' --to-hdf5 
	biom convert -i './songbird/WMS_CRC_control_Comparison'$run'_Subset2_otutable.tsv' -o './songbird/WMS_CRC_control_Comparison'$run'_Subset2_otutable.biom' --to-hdf5 
	biom convert -i './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset1_otutable.tsv' -o './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset1_otutable.biom' --to-hdf5 
	biom convert -i './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset2_otutable.tsv' -o './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset2_otutable.biom' --to-hdf5 
	biom convert -i './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset1_otutable.tsv' -o './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset1_otutable.biom' --to-hdf5 
	biom convert -i './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset2_otutable.tsv' -o './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset2_otutable.biom' --to-hdf5 
	
	mkdir './songbird/WMS_CRC_control_Comparison'$run'_Subset1'
	songbird multinomial --input-biom './songbird/WMS_CRC_control_Comparison'$run'_Subset1_otutable.biom' \
                         --metadata-file './songbird/WMS_CRC_control_Comparison'$run'_Subset1_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/WMS_CRC_control_Comparison'$run'_Subset1'
	
	mkdir './songbird/WMS_CRC_control_Comparison'$run'_Subset2'
	songbird multinomial --input-biom './songbird/WMS_CRC_control_Comparison'$run'_Subset2_otutable.biom' \
                         --metadata-file './songbird/WMS_CRC_control_Comparison'$run'_Subset2_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/WMS_CRC_control_Comparison'$run'_Subset2'
	
	mkdir './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset1'
	songbird multinomial --input-biom './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset1_otutable.biom' \
                         --metadata-file './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset1_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset1'
	
	mkdir './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset2'
	songbird multinomial --input-biom './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset2_otutable.biom' \
                         --metadata-file './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset2_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/WMS_schizophrenia_control_Comparison'$run'_Subset2'
	
	mkdir './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset1'
	songbird multinomial --input-biom './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset1_otutable.biom' \
                         --metadata-file './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset1_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset1'
	
	mkdir './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset2'
	songbird multinomial --input-biom './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset2_otutable.biom' \
                         --metadata-file './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset2_samdata.tsv' \
                         --formula "grp" \
                         --summary-dir './songbird/WMS_tonguedorsum_stool_Comparison'$run'_Subset2'
	
	run=$((run+1)); 
done