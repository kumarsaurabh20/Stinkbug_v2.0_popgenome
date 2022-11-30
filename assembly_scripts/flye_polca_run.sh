#!/bin/bash -l
#SBATCH --partition=defq
#SBATCH --job-name=polca_SB                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=64                              # Number of CPU cores per task
#SBATCH --mem=2gb                                       # Job memory request
#SBATCH --time=24:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=polish_polca_%j.log                        # Standard output and error log
#SBATCH --account=c.bass
##
pwd; hostname; date
export OMP_NUM_THREADS=64
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module load slurm/18.08.4 DefaultModules shared Workspace/v1
module load bwa/0.7.17 samtools/1.9
module load ks575/MaSurCA/v3.4.2
##export PATH=$PATH:$WORKSPACE/Data/Tuta/software/MaSuRCA-3.4.2/global-1/PacBio/src_reconcile
base=$WORKSPACE/Data/Stinkbug/flye/Polish
run1=$base/polca/run1
run2=$base/polca/run2
run3=$base/polca/run3
genome=$base/assembly.fasta
##read1=$WORKSPACE/Data/Tuta/LIB30793_R1.fastq
##read2=$WORKSPACE/Data/Tuta/LIB30793_R2.fastq
##Run1
cd $run1
polca.sh -a $genome -r '/nobackup/beegfs/workspace/ks575/Data/Stinkbug/2730_1B_trimmed_r1.fastq.gz /nobackup/beegfs/workspace/ks575/Data/Stinkbug/2730_1B_trimmed_r2.fastq.gz' -t 64 -m 1G
cp $run1/assembly.fasta.PolcaCorrected.fa $run1/final_PC_1.fa
cp $run1/final_PC_1.fa $run2/final_PC_1.fa
##Run2
cd $run2
genome2=$run2/final_PC_1.fa
polca.sh -a $genome2 -r '/nobackup/beegfs/workspace/ks575/Data/Stinkbug/2730_1B_trimmed_r1.fastq.gz /nobackup/beegfs/workspace/ks575/Data/Stinkbug/2730_1B_trimmed_r2.fastq.gz' -t 64 -m 1G
cp $run2/final_PC_1.fa.PolcaCorrected.fa $run2/final_PC_2.fa
cp $run2/final_PC_2.fa $run3/final_PC_2.fa
##Run 3
cd $run3
genome3=$run3/final_PC_2.fa
polca.sh -a $genome3 -r '/nobackup/beegfs/workspace/ks575/Data/Stinkbug/2730_1B_trimmed_r1.fastq.gz /nobackup/beegfs/workspace/ks575/Data/Stinkbug/2730_1B_trimmed_r2.fastq.gz' -t 64 -m 1G
cp $run3/final_PC_2.fa.PolcaCorrected.fa $run3/final_PC_3.fa
cp $run3/final_PC_3.fa $base/final_SR_only_polishing_3R.fasta
