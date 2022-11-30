#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=psmc.a                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=16                         # Number of CPU cores per task
#SBATCH --mem=10gb                                       # Job memory request
#SBATCH --time=340:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=psmc.a_parallel_%j.log                        # Standard output and error log
#SBATCH --account=c.bass

pwd; hostname; date
export OMP_NUM_THREADS=16
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module load DefaultModules shared Workspace/v1 slurm/18.08.4  Conda/Python2/2.7.15 ks575/BioPython/v1.74 VCFtools/0.1.16 ActivePerl/5.24.3.2404 bcftools/1.9 samtools/1.10
#
export PATH=$PATH:/nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/VCF/final/filtering/PSMC/psmc:/nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/VCF/final/filtering/PSMC/psmc/utils
work_dir=$WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/psmc/data
cd $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/psmc/data
#

psmc -N30 -t 15 -r5 -p "4+30*2" -o $work_dir/SPA10.psmc $work_dir/SPA10mod.psmcfa
#psmc -N25 -t 15 -r5 -p "4+30*2+4+6+10" -o $work_dir/SPA10.psmc $work_dir/SPA10mod.psmcfa
#psmc -N25 -t 15 -r5 -p "6+2*4+3+13*2+3+2*4+6" -o $work_dir/SPA10.psmc $work_dir/SPA10mod.psmcfa
