#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=emiscor2                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=64                              # Number of CPU cores per task
#SBATCH --mem=2gb                                       # Job memory request
#SBATCH --time=400:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=missassembly_correct_%j.log                        # Standard output and error log
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
work_dir=$WORKSPACE/Data/Stinkbug/qm_dbg2olc_flye/default/Polish/LR+SR/purge/misassembly
assembly=$work_dir/curated.fasta
index=$work_dir/curated.fasta.fai
misassembly_file=$work_dir/curated.missassemblies.tsv
outfile="curated.corrected.fasta"
##xargs (bash) wont accept local variables so hardcode the bam file in the main script
cd $work_dir
if [ -s $misassembly_file ] #if the input is nonempty...
                then
                        cat $misassembly_file | grep -v '^#' | sort -k1,1 -k2,2g | join - <(cat $index | sort -k1,1 -k2,2g) | \
                        awk '{{
                                if ($1 == prev_tig){{
                                 print($1,prev_coord,$2)
                                 }}
                                else{{
                                        if (prev_len > 0){{
                                                print(prev_tig,prev_coord,prev_len)
                                        }}
                                        print($1,"1",$2)
                                        }}
                                prev_tig = $1
                                prev_coord = $2
                                prev_len = $4
                                }}
                                END {{ print(prev_tig,prev_coord,prev_len)
                        }}' | sed "s/\(.*\) \(.*\)\ \(.*\)/\\1:\\2-\\3/g" | xargs samtools faidx /nobackup/beegfs/workspace/ks575/Data/Stinkbug/qm_dbg2olc_flye/default/Polish/LR+SR/purge/misassembly/curated.fasta | cut -f1 -d ':' | awk '(/^>/ && s[$0]++){{$0=$0"_"s[$0]}}1;' > $work_dir/curated.corrected.fasta
		#	cut -f1 $work_dir/curated.corrected.fasta > $work_dir/curated.tigs.toremove
		#	grep -vf $work_dir/curated.tigs.toremove $index | cut -f1 | xargs samtools faidx $assembly >> $work_dir/curated.corrected.fasta
		#	rm $work_dir/curated.tigs.toremove
                else
                        cp $assembly $work_dir/curated.corrected.fasta
                fi
