#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=stacks                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=8                           # Number of CPU cores per task
#SBATCH --mem=80g                                       # Job memory request
#SBATCH --time=240:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=bugs_%j.log                        # Standard output and error log
#SBATCH --account=c.bass

pwd; hostname; date
export OMP_NUM_THREADS=8
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#
module load Workspace/v1 bwa/0.7.17 samtools/1.9
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
work_dir=$WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus
index=$work_dir/index/SB_corrected
out_dir=$work_dir/mapping/bugs
#
#find $(pwd) -name "*_R1_val_1.fq" > file1.txt
#find $(pwd) -name "*_R2_val_2.fq" > file2.txt
#paste file1.txt file2.txt > file.txt
while IFS=$'\t' read -r f1
do
        echo "Started the sample $f1"
        outfile=$(echo $f1 | cut -d'/' -f11 | cut -d'.' -f1)
        bwa mem -t 8 -M $index "$f1" | samtools view -bS -F 4 - | samtools sort -@ 8 -m 10G -o $out_dir/${outfile}.bam -
done < $work_dir/bugs.list.txt 
date
