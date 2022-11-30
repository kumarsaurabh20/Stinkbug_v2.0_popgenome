#!/bin/bash -l
#SBATCH --partition=defq
#SBATCH --job-name=qm_SB1                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=10                              # Number of CPU cores per task
#SBATCH --mem=10gb                                       # Job memory request
#SBATCH --time=40:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=qm_dbg2olc_flye_%j.log                        # Standard output and error log
#SBATCH --account=c.bass
pwd; hostname; date
export OMP_NUM_THREADS=10
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
module load Workspace/v1 quickmerge/0.3
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
cd /nobackup/beegfs/workspace/ks575/Data/Stinkbug/qm_dbg2olc_flye/default
work_dir=/nobackup/beegfs/workspace/ks575/Data/Stinkbug/qm_dbg2olc_flye/default
merge_wrapper.py $work_dir/final_LR+SR_polishing.fasta $work_dir/final_SR_only_polishing_3R.fasta
date
