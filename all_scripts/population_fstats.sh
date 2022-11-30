#!/bin/bash -l
#SBATCH --partition=gpuq
#SBATCH --job-name=stacks                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=24                           # Number of CPU cores per task
#SBATCH --mem=8g                                       # Job memory request
#SBATCH --time=240:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=bugs_%j.log                        # Standard output and error log
#SBATCH --account=c.bass

pwd; hostname; date
export OMP_NUM_THREADS=24
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#
module load Workspace/v1 bwa/0.7.17 samtools/1.9
module load ks575/Stacks/v2.55
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
work_dir=$WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus
index=$work_dir/index/SB_corrected
out_dir=$work_dir/mapping_new/fstats
input_folder=$work_dir/mapping_new/fstats
map_file=$work_dir/mapping_new/filtered_sets.txt
log_file=$out_dir/populations_fstats.oe
#
#find $(pwd) -name "*_R1_val_1.fq" > file1.txt
#find $(pwd) -name "*_R2_val_2.fq" > file2.txt
#paste file1.txt file2.txt > file.txt
#gstacks -I $input_folder -S ".RG.bam" -M $map_file -O $out_dir -t $OMP_NUM_THREADS --max-clipped 0.1
#populations -P $input_folder -M $map_file -t $OMP_NUM_THREADS -p 2 -r 0.6 --ordered-export --fasta-loci --vcf
#populations -P $input_folder -M $map_file -t $OMP_NUM_THREADS -p 2 -r 0.8 --min-maf 0.05 --max-obs-het 0.70 --ordered-export --fasta-loci --vcf --genepop
populations -P $input_folder -M $map_file -t $OMP_NUM_THREADS -O $out_dir -p 2 -r 0.80 --min-maf 0.05 --max-obs-het 0.70 --fstats -k --sigma 100000
