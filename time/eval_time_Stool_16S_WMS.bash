#!/bin/bash
#SBATCH --job-name=eval_time   			# Job name
#SBATCH -p normal   					# priority
#SBATCH --mail-type=ALL         		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mcalgaro93@gmail.com # Where to send mail	
#SBATCH --nodes=1                   	# Use one node
#SBATCH --ntasks=1                  	# Run a single task
#SBATCH --mem-per-cpu=5gb           	# Memory per processor
#SBATCH --time=24:00:00             	# Time limit hrs:min:sec
#SBATCH --output=array_%A-%a.out    	# Standard output and error log
#SBATCH --array=1-30                 	# Array range

#Set the number of runs that each SLURM task should do
PER_TASK=1
# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))
# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

module purge
module load R/3.5.2

# Run the loop of runs for this task.
run=$START_NUM
while [ $run -le $END_NUM ]; do 
	echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run; 
	Rscript eval_time_Stool_16S_WMS_script.R ${run};
	run=$((run+1)); 
done
