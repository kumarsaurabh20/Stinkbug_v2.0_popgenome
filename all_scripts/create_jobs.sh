#!/bin/bash
work_dir="/nobackup/beegfs/workspace/ks575/Data/Stinkbug/SB_population_data/Euschistus/mapping_new/alignments"
i=1
for each in $(find $(pwd) -name "*.bam")
do
	touch jobs.txt
	sample=$(echo $each | cut -d'/' -f12 | cut -d'.' -f1) 	
	echo "gatk AddOrReplaceReadGroups --INPUT=$each --OUTPUT="$work_dir/${sample}.RG.bam" --RGID=$i --RGLB=$sample --RGPL=illumina --RGPU=unit1 --RGSM=$sample --SORT_ORDER=coordinate" >> jobs.txt
	((i++))
done
