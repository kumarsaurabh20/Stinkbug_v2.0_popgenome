while IFS=$'\t' read -r f1 f2
do
        echo "Started the sample $f1"
        outfile=$(echo $f1 | cut -d'/' -f12 | cut -d'.' -f1)
        #bwa mem -t 32 -M $index "$f1" | samtools view -bS -F 4 - | samtools sort -@ 32 -m 2G -o $out_dir/${outfile}.bam -
	echo "$outfile\t$f2"
done < bugs.pstacks.list
