#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=SBblobs1                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=64                         # Number of CPU cores per task
#SBATCH --mem=10gb                                       # Job memory request
#SBATCH --time=240:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=blobs_%j.log                        # Standard output and error log
#SBATCH --account=c.bass
##
pwd; hostname; date
export OMP_NUM_THREADS=64
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module purge
module load slurm/18.08.4 DefaultModules shared Workspace/v1
module load ks575/minimap2/v2.17 samtools/1.9
module load ks575/BlobTools/v1.1.1
module load BLAST+/2.8.1
##
work_dir=$WORKSPACE/Data/Stinkbug/qm_dbg2olc_flye/default/Polish/LR+SR/purge/misassembly/blobs
##
genome=$WORKSPACE/Data/Stinkbug/qm_dbg2olc_flye/default/Polish/LR+SR/purge/misassembly/curated.fasta
pacbio=$WORKSPACE/Data/Stinkbug/SB_all_reads.fasta
blastDB=$WORKSPACE/Data/Tools/blast_dbs/nt_db/nt
nodes=$WORKSPACE/Data/Tools/blast_dbs/tax_db/nodes.dmp
names=$WORKSPACE/Data/Tools/blast_dbs/tax_db/names.dmp
##
##minimap2 -t 64 -ax map-pb $genome $pacbio | samtools view -Sb -h -F 4 - | samtools sort -@ 64 -m 1G -o $work_dir/aligned.bam -T tmp.ali
##
blastn -task megablast -query $genome -db $blastDB -culling_limit 2 -out $work_dir/curated.corrected.blastn -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -negative_gilist $work_dir/sequence.gi -evalue 1e-25 -num_threads 64
##
blobtools create -i $genome -b $work_dir/aligned.bam -t $work_dir/curated.corrected.blastn --nodes $nodes --names $names
##
blobtools view -i $work_dir/blobDB.json --hits --rank all > $work_dir/blobDB.table
##
blobtools blobplot -i $work_dir/blobDB.json
