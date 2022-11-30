if [ -s curated.missassemblies.tsv ] #if the input is nonempty...
                then
                        cat curated.missassemblies.tsv | grep -v '^#' | sort -k1,1 -k2,2g | join - <(cat curated.fasta.fai | sort -k1,1 -k2,2g) | \
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
                        }}' | sed "s/\(.*\) \(.*\)\ \(.*\)/\\1:\\2-\\3/g" | xargs samtools faidx curated.fasta | cut -f1 -d ':' | awk '(/^>/ && s[$0]++){{$0=$0"_"s[$0]}}1;' > curated.corrected.fasta

                       # cut -f1 curated.corrected.fasta > curated.tigs.toremove
                       # grep -vf curated.tigs.toremove curated.fasta.fai | cut -f1 | xargs samtools faidx curated.fasta >> curated.corrected.fasta
                       # rm curated.tigs.toremove
                else
			cp curated.fasta curated.corrected.fasta
                fi

