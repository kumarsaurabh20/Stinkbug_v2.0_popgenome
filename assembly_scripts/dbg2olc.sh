#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=dbg2_SB                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=45                           # Number of CPU cores per task
#SBATCH --mem=10gb                                       # Job memory request
#SBATCH --time=240:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=dbg2olc_%j.log                        # Standard output and error log
#SBATCH --account=c.bass

pwd; hostname; date
export OMP_NUM_THREADS=45
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#
module purge
module load slurm/18.08.4 shared DefaultModules Workspace/v1 
module load Conda/Python2/2.7.15
module load ks575/DBG2OLC/v20180222
module load BLASR/5.3.2
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
cd $WORKSPACE/Data/Stinkbug/DBG2OLC
work_dir=$WORKSPACE/Data/Stinkbug/DBG2OLC
##
SparseAssembler LD 0 k 51 g 15 NodeCovTh 2 EdgeCovTh 1 GS 1130000000 i1 /nobackup/beegfs/workspace/ks575/Data/Stinkbug/2730_1B_trimmed_r1.fastq i2 /nobackup/beegfs/workspace/ks575/Data/Stinkbug/2730_1B_trimmed_r2.fastq
##
DBG2OLC k 17 AdaptiveTh 0.01 KmerCovTh 5 MinOverlap 50 RemoveChimera 1 Contigs $work_dir/Contigs.txt f $WORKSPACE/Data/Stinkbug/SB_all_reads.fasta
##
mkdir -p $work_dir/consensus_dir
cat $WORKSPACE/Data/Stinkbug/SB_all_reads.fasta $work_dir/Contigs.txt > $work_dir/Contigs_pb.fasta
##
split_and_run_sparc.2.sh $work_dir/backbone_raw.fasta $work_dir/DBG2OLC_Consensus_info.txt $work_dir/Contigs_pb.fasta $work_dir/consensus_dir 2
