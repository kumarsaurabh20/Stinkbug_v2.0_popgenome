#!/bin/bash -l
#SBATCH --partition=hmq
#SBATCH --job-name=psmc8
#SBATCH --mail-type=FAIL 
#SBATCH --mail-user=ks575@exeter.ac.uk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --time=240:05:00
#SBATCH --output=psmc8_parallel_%j.log
#SBATCH --account=c.bass
##
pwd; hostname; date
export OMP_NUM_THREADS=1
pwd; hostname; date
echo "Running a program on "
##
#/path/to/my/application
echo "Hello $USER, this is running on the $SLURM_CLUSTER_NAME cluster at `hostname` using PI account = $SLURM_JOB_ACCOUNT"
module load DefaultModules shared Workspace/v1 slurm/18.08.4  Conda/Python2/2.7.15 ks575/BioPython/v1.74 VCFtools/0.1.16 ActivePerl/5.24.3.2404 bcftools/1.9 samtools/1.10
samtools mpileup -C50 -uf /nobackup/beegfs/workspace/ks575/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/Final_corrected.fasta /nobackup/beegfs/workspace/ks575/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/RSSB15.RG.bam | bcftools call -c - | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > /nobackup/beegfs/workspace/ks575/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments/PSMC/RSSB15.diploid.fq.gz
