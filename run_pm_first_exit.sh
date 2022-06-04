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

julia -p 10 empirical_cdfs.jl hitting_times_model/reflecting_model/model_def.jl pm_initial_condition_exit_times_sim_defs.jl 50000000000 1000 16860 pm_ic_exit_time_50_000_000_000.json pm_ic_exit_time_50_000_000_000_boot.json
