#!/bin/bash -l
#SBATCH --partition=defq
#SBATCH --job-name=purge13                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=64                         # Number of CPU cores per task
#SBATCH --mem=10gb                                       # Job memory request
#SBATCH --time=240:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=purge_parallel_%j.log                        # Standard output and error log
#SBATCH --account=c.bass

pwd; hostname; date
export OMP_NUM_THREADS=64
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
##
##/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module load slurm/18.08.4 DefaultModules shared gcc/8.2.0 Workspace/v1 Oracle_Java/8u192
module load ks575/minimap2/v2.17
module load samtools/1.9
module load ks575/Purge-haplotigs/v1.0.4
base=/nobackup/beegfs/workspace/ks575/Data/Stinkbug/qm_dbg2olc_flye/default/Polish/LR+SR
work_dir=$base/purge
genome=$base/final_LR+SR_polished.fasta
pacbio=$WORKSPACE/Data/Stinkbug/SB_all_reads.fasta
minimap2 -t 64 -ax map-pb $genome $pacbio | samtools view -hF 256 - | samtools sort -@ 64 -m 1G -o $work_dir/aligned.bam -T tmp.ali
##
purge_haplotigs readhist -t 64 -b $work_dir/aligned.bam -g $genome
