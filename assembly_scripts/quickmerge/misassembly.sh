#!/bin/bash -l
#SBATCH --partition=defq
#SBATCH --job-name=emiscorr                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=64                         # Number of CPU cores per task
#SBATCH --mem=1gb                                       # Job memory request
#SBATCH --time=240:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=emiscorr_%j.log                        # Standard output and error log
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
work_dir=$WORKSPACE/Data/Stinkbug/qm_dbg2olc_flye/default/Polish/LR+SR/purge/misassembly
genome=$work_dir/curated.fasta
pacbio=$WORKSPACE/Data/Stinkbug/SB_all_reads.fasta
minimap2 -t 64 -ax map-pb $genome $pacbio | samtools view -Sb -h -F 4 - | samtools sort -@ 64 -m 1G -o $work_dir/aligned.bam -T tmp.ali
##
