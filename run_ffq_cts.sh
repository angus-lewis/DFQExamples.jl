#!/bin/bash
#SBATCH -p batch            	                            # partition (this is the queue your job will be added to) 
#SBATCH -N 1               	                                # number of nodes (no MPI, so we only use a single node)
#SBATCH -n 11            	                                # number of cores
#SBATCH --time=2-20:00:00  	                                # walltime allocation, which has the format (D-HH:MM:SS), here set to 1 hour
#SBATCH --mem=4GB         	                                # memory required per node (here set to 4 GB)

# Notification configuration 
#SBATCH --mail-type=END					    	# Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL   					# Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=angus.lewis@adelaide.edu.au  	# Email to which notifications will be sent

# Execute the program
module load Julia/1.6.0

julia -p 10 empirical_cdfs.jl fluidfluid_continuous/model_def.jl ffq_cont_first_return_sim_defs.jl 50000000000 1000 16860 ffq_cont_first_return_50_000_000_000.json ffq_cont_first_return_50_000_000_000_boot.json
