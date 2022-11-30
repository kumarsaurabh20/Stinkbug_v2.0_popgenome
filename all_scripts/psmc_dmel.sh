#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=psmc.a                         # Job name
#SBATCH --mail-type=FAIL                            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ks575@exeter.ac.uk                  # Where to send mail
#SBATCH --ntasks=1                                      # Run a single task
#SBATCH --cpus-per-task=16                         # Number of CPU cores per task
#SBATCH --mem=10gb                                       # Job memory request
#SBATCH --time=340:05:00                                 # Time limit hrs:min:sec
#SBATCH --output=psmc.a_parallel_%j.log                        # Standard output and error log
#SBATCH --account=c.bass

pwd; hostname; date
export OMP_NUM_THREADS=16
pwd; hostname; date
echo "Running a program on $SLURM_JOB_NODELIST"
#
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module load DefaultModules shared Workspace/v1 slurm/18.08.4  Conda/Python2/2.7.15 ks575/BioPython/v1.74 VCFtools/0.1.16 ActivePerl/5.24.3.2404 bcftools/1.9 samtools/1.10
#
export PATH=$PATH:/nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/VCF/final/filtering/PSMC/psmc:/nobackup/beegfs/workspace/ks575/Data/Myzus/Myzus_mapping/Trimmed_Resequencing_Data_18_03_2019/VCF/final/filtering/PSMC/psmc/utils
cd $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC
#
fq2psmcfa -q20 $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/MGC1.diploid.fq.gz > $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/MGC1.dm.psmcfa
psmc -N25 -t 15 -r5 -p "4+25*2+4+6" -o $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/MGC1.dm.psmc $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/MGC1.dm.psmcfa
#
fq2psmcfa -q20 $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/MTR39.diploid.fq.gz > $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/MTR39.dm.psmcfa
psmc -N25 -t 15 -r5 -p "4+25*2+4+6" -o $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/MTR39.dm.psmc $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/MTR39.dm.psmcfa
#
fq2psmcfa -q20 $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/GOPB10.diploid.fq.gz > $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/GOPB10.dm.psmcfa
psmc -N25 -t 15 -r5 -p "4+25*2+4+6" -o $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/GOPB10.dm.psmc $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/GOPB10.dm.psmcfa
#
fq2psmcfa -q20 $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/BASD31.diploid.fq.gz > $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/BASD31.dm.psmcfa
psmc -N25 -t 15 -r5 -p "4+25*2+4+6" -o $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/BASD31.dm.psmc $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/BASD31.dm.psmcfa
#
fq2psmcfa -q20 $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/GOPB16.diploid.fq.gz > $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/GOPB16.dm.psmcfa
psmc -N25 -t 15 -r5 -p "4+25*2+4+6" -o $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/GOPB16.dm.psmc $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/GOPB16.dm.psmcfa
#
fq2psmcfa -q20 $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/RSSB33.diploid.fq.gz > $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/RSSB33.dm.psmcfa
psmc -N25 -t 15 -r1.12 -p "4+25*2+4+6" -o $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/RSSB33.dm.psmc $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/RSSB33.dm.psmcfa
#
fq2psmcfa -q20 $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/PRIV36.diploid.fq.gz > $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/PRIV36.dm.psmcfa
psmc -N25 -t 15 -r5 -p "4+25*2+4+6" -o $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/PRIV36.dm.psmc $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/PRIV36.dm.psmcfa
#
fq2psmcfa -q20 $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/RSSB15.diploid.fq.gz > $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC//RSSB15.dm.psmcfa
psmc -N25 -t 15 -r5 -p "4+25*2+4+6" -o $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC//RSSB15.dm.psmc $WORKSPACE/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/RSSB15.dm.psmcfa
