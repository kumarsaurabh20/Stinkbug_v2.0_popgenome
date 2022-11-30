#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=SB_flye                   # Job name
#SBATCH --mail-type=FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk       # Where to send mail
#SBATCH --ntasks=1                           # Run a single task
#SBATCH --cpus-per-task=20                  # Number of CPU cores per task
#SBATCH --mem=20gb
#SBATCH --time=300:05:00                     # Time limit hrs:min:sec
#SBATCH --output=sbflye1_parallel_%j.log      # Standard output and error log
#SBATCH --account=c.bass

pwd; hostname; date
export OMP_NUM_THREADS=20
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module load shared slurm Workspace/v1 Conda/Python2/2.7.15 ks575/Flye/v2.6.0 ks575/minimap2/v2.17
export PATH=$PATH:$WORKSPACE/Data/software/ra/build/bin
cd $WORKSPACE/Data/Stinkbug/SB_pacbio2
flye --pacbio-raw $WORKSPACE/Data/Stinkbug/split_pacbio/SB_all_reads.part_002.fasta.gz --out-dir $WORKSPACE/Data/Stinkbug/SB_pacbio2 --genome-size 890m --threads 20 --asm-coverage 40
