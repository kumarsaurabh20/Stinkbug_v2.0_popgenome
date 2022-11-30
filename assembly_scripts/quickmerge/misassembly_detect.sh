#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=emiscorr                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=64                              # Number of CPU cores per task
#SBATCH --mem=2gb                                       # Job memory request
#SBATCH --time=400:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=missassembly_detect_%j.log                        # Standard output and error log
#SBATCH --account=c.bass
pwd; hostname; date
export OMP_NUM_THREADS=64
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module load slurm/18.08.4 DefaultModules shared Workspace/v1
module load samtools/1.9
module load ks575/BedTools/v2.9.1
module load ks575/HTS-box/v346
work_dir=/nobackup/beegfs/workspace/ks575/Data/Stinkbug/qm_dbg2olc_flye/default/Polish/LR+SR/purge/misassembly
assembly=curated.fasta
bamfile=aligned.bam
##xargs (bash) wont accept local variables so hardcode the bam file in the main script
cd $work_dir
bedtools makewindows -g $work_dir/${assembly}.fai -w 2000 | join - $work_dir/${assembly}.fai | tr ' ' '\t' | cut -f1-4 | awk '{{if ($2 > 2000 && $3 < $4 - 2000 && $4 > 50000) print $0}}' | xargs -P 64 -l bash -c ' htsbox samview aligned.bam $0:$1-$1 -p | cut -f8,9 | awk "{{if (\$1 < $1 - (2000/2) && \$2 > $1 + (2000/2)) print \$0}}" | wc -l | paste <(echo $0) <(echo $1) - ' | awk '{{if ($3 < 2) print $0}} ' > $work_dir/curated.missassemblies.tsv
