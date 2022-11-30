#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=iqtree                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=64                         # Number of CPU cores per task
#SBATCH --mem=5gb                                       # Job memory request
#SBATCH --time=500:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=iqtree_parallel_%j.log                        # Standard output and error log
#SBATCH --account=c.bass
##
pwd; hostname; date
export OMP_NUM_THREADS=64
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
##
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module load DefaultModules shared Workspace/v1 slurm/18.08.4  Conda/Python2/2.7.15 ks575/BioPython/v1.74 VCFtools/0.1.16 
module load bcftools/1.9 samtools/1.9 ks575/IQ-tree/v2.0.3
##
cd /nobackup/beegfs/workspace/ks575/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/gstacks/demography/filtered
work_dir="/nobackup/beegfs/workspace/ks575/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/gstacks/demography/filtered"
iqtree -s $work_dir/SB.SNP.sorted.renamed.filtered.min79.phy -B 1000 -alrt 1000 -bnni -T AUTO
