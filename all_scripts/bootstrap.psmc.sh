#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=psmc.b                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=32                         # Number of CPU cores per task
#SBATCH --mem=10gb                                       # Job memory request
#SBATCH --time=340:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=psmc.a_parallel_%j.log                        # Standard output and error log
#SBATCH --account=c.bass

pwd; hostname; date
export OMP_NUM_THREADS=32
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module load DefaultModules shared Workspace/v1 slurm/18.08.4  Conda/Python2/2.7.15 ks575/BioPython/v1.74 VCFtools/0.1.16 ActivePerl/5.24.3.2404 bcftools/1.9 samtools/1.10
#
export PATH=$PATH:/nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/VCF/final/filtering/PSMC/psmc:/nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/VCF/final/filtering/PSMC/psmc/utils
work_dir=$WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/psmc/data
cd $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/psmc/data
#
splitfa $work_dir/BASD31mod.psmcfa > $work_dir/split.psmcfa
psmc -N30 -t15 -r5 -p "64*1" -o $work_dir/BASD31.psmc $work_dir/BASD31mod.psmcfa
#seq 100 | xargs -n 1 -p 32 -i echo psmc -N30 -t15 -r5 -b -p "64*1" -o $work_dir/round-{}.psmc $work_dir/split.psmcfa | sh
seq 1 100 | xargs -n 1 -P 32 -i $work_dir/bash.sh $work_dir/split.psmcfa $work_dir/round-{}.psmc
cat $work_dir/BASD31.psmc $work_dir/round-*.psmc > $work_dir/combined.psmc
rm -vrf $work_dir/split.psmcfa $work_dir/BASD31.psmc $work_dir/round-*
mv $work_dir/combined.psmc $work_dir/BASD31.psmc
mv $work_dir/BASD31.psmc $work_dir/psmc/
mv $work_dir/BASD31mod.psmcfa $work_dir/psmcfa/
mv $work_dir/BASD31.psmcfa $work_dir/psmcfa/
##
splitfa $work_dir/GOPB22mod.psmcfa > $work_dir/split.psmcfa
psmc -N30 -t15 -r5 -p "64*1" -o $work_dir/GOPB22.psmc $work_dir/GOPB22mod.psmcfa
#seq 100 | xargs -n 1 -p 32 -i echo psmc -N30 -t15 -r5 -b -p "64*1" -o $work_dir/round-{}.psmc $work_dir/split.psmcfa | sh
seq 1 100 | xargs -n 1 -P 32 -i $work_dir/bash.sh $work_dir/split.psmcfa $work_dir/round-{}.psmc
cat $work_dir/GOPB22.psmc $work_dir/round-*.psmc > $work_dir/combined.psmc
rm -vrf $work_dir/split.psmcfa $work_dir/GOPB22.psmc $work_dir/round-*
mv $work_dir/combined.psmc $work_dir/GOPB22.psmc
mv $work_dir/GOPB22.psmc $work_dir/psmc/
mv $work_dir/GOPB22mod.psmcfa $work_dir/psmcfa/
mv $work_dir/GOPB22.psmcfa $work_dir/psmcfa/
##
splitfa $work_dir/MTR39mod.psmcfa > $work_dir/split.psmcfa
psmc -N30 -t15 -r5 -p "64*1" -o $work_dir/MTR39.psmc $work_dir/MTR39mod.psmcfa
#seq 100 | xargs -n 1 -p 32 -i echo psmc -N30 -t15 -r5 -b -p "64*1" -o $work_dir/round-{}.psmc $work_dir/split.psmcfa | sh
seq 1 100 | xargs -n 1 -P 32 -i $work_dir/bash.sh $work_dir/split.psmcfa $work_dir/round-{}.psmc
cat $work_dir/MTR39.psmc $work_dir/round-*.psmc > $work_dir/combined.psmc
rm -vrf $work_dir/split.psmcfa $work_dir/MTR39.psmc $work_dir/round-*
mv $work_dir/combined.psmc $work_dir/MTR39.psmc
mv $work_dir/MTR39.psmc $work_dir/psmc/
mv $work_dir/MTR39mod.psmcfa $work_dir/psmcfa/
mv $work_dir/MTR39.psmcfa $work_dir/psmcfa/
##
splitfa $work_dir/PIB7mod.psmcfa > $work_dir/split.psmcfa
psmc -N30 -t15 -r5 -p "64*1" -o $work_dir/PIB7.psmc $work_dir/PIB7mod.psmcfa
#seq 100 | xargs -n 1 -p 32 -i echo psmc -N30 -t15 -r5 -b -p "64*1" -o $work_dir/round-{}.psmc $work_dir/split.psmcfa | sh
seq 1 100 | xargs -n 1 -P 32 -i $work_dir/bash.sh $work_dir/split.psmcfa $work_dir/round-{}.psmc
cat $work_dir/PIB7.psmc $work_dir/round-*.psmc > $work_dir/combined.psmc
rm -vrf $work_dir/split.psmcfa $work_dir/PIB7.psmc $work_dir/round-*
mv $work_dir/combined.psmc $work_dir/PIB7.psmc
mv $work_dir/PIB7.psmc $work_dir/psmc/
mv $work_dir/PIB7mod.psmcfa $work_dir/psmcfa/
mv $work_dir/PIB7.psmcfa $work_dir/psmcfa/
##
splitfa $work_dir/ROC31mod.psmcfa > $work_dir/split.psmcfa
psmc -N30 -t15 -r5 -p "64*1" -o $work_dir/ROC31.psmc $work_dir/ROC31mod.psmcfa
#seq 100 | xargs -n 1 -p 32 -i echo psmc -N30 -t15 -r5 -b -p "64*1" -o $work_dir/round-{}.psmc $work_dir/split.psmcfa | sh
seq 1 100 | xargs -n 1 -P 32 -i $work_dir/bash.sh $work_dir/split.psmcfa $work_dir/round-{}.psmc
cat $work_dir/ROC31.psmc $work_dir/round-*.psmc > $work_dir/combined.psmc
rm -vrf $work_dir/split.psmcfa $work_dir/ROC31.psmc $work_dir/round-*
mv $work_dir/combined.psmc $work_dir/ROC31.psmc
mv $work_dir/ROC31.psmc $work_dir/psmc/
mv $work_dir/ROC31mod.psmcfa $work_dir/psmcfa/
mv $work_dir/ROC31.psmcfa $work_dir/psmcfa/
