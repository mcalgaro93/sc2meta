#!/bin/bash
#SBATCH --job-name=16S_concordance   			# Job name
#SBATCH -p normal   					# priority
#SBATCH --mail-type=ALL         		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=mcalgaro93@gmail.com # Where to send mail	
#SBATCH --nodes=1                   	# Use one node
#SBATCH --ntasks=1                  	# Run a single task
#SBATCH --mem-per-cpu=5gb           	# Memory per processor
#SBATCH --time=24:00:00             	# Time limit hrs:min:sec
#SBATCH --output=concordance_%A-%a.out    	# Standard output and error log

module purge
module load R/3.5.2

Rscript compute_concordance16S.R

