#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=racon_SB                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=64                              # Number of CPU cores per task
#SBATCH --mem=5gb                                       # Job memory request
#SBATCH --time=240:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=polish_racon_%j.log                        # Standard output and error log
#SBATCH --account=c.bass

pwd; hostname; date
export OMP_NUM_THREADS=64
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module load slurm/18.08.4 DefaultModules shared Workspace/v1
module load ks575/Racon/v1.4.10
export PATH=$PATH:/cm/shared/user-apps/ks575/minimap2/v2.17/bin
cd $WORKSPACE/Data/Stinkbug/flye/Polish
base=$WORKSPACE/Data/Stinkbug/flye/Polish
genome=$base/assembly.fasta
pacbio=$WORKSPACE/Data/Stinkbug/SB_all_reads.fasta
run1=$base/racon/run1
run2=$base/racon/run2
run3=$base/racon/run3
out="racon"
##
minimap2 -t 64 -ax map-pb $genome $pacbio > $run1/aligned.sam
racon -t 64 $pacbio $run1/aligned.sam $genome > $run1/${out}_1_polished.fasta
##
minimap2 -t 64 -ax map-pb $run1/${out}_1_polished.fasta $pacbio > $run2/aligned.sam
racon -t 64 $pacbio $run2/aligned.sam $run1/${out}_1_polished.fasta > $run2/${out}_2_polished.fasta
##
minimap2 -t 64 -ax map-pb $run2/${out}_2_polished.fasta $pacbio > $run3/aligned.sam
racon -t 64 $pacbio $run3/aligned.sam $run2/${out}_2_polished.fasta > $run3/${out}_3_polished.fasta
##
cp $run3/${out}_3_polished.fasta $base/final_LR_polished_3R.fasta
