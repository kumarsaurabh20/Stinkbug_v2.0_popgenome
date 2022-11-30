#!/bin/bash
# Launch populations
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script as it was run
SCRIPT=$0
NAME=$(basename $0)
SAMPLE_FOLDER="04-all_samples"
STACKS_FOLDER="05-stacks"
LOG_FOLDER="10-log_files"
INFO_FILES_FOLDER="01-info_files"
POP_MAP="population_map.txt"

# Number of CPUs
# TODO use gnu parallel
NUM_CPU="32"

#id=1
#ls -1 "$SAMPLE_FOLDER"/*.gz |
#while read file
#do
#    name=$(basename "$file")
#    echo
#    echo "#######################################"
#    echo "  Treating indivudual $id: ${name%.fq.gz}"
#    echo "#######################################"
#    echo

#    ustacks -f "$SAMPLE_FOLDER"/"$name" -o "$STACKS_FOLDER" -i $id -p 1 \
#        -m 4 -M 3 -N 5 -H --deleverage
#    let "id+=1"

#done

# Gnu Parallel version
ls -1 "$SAMPLE_FOLDER"/*.gz |
    parallel -j "$NUM_CPU" ustacks -f {} -o "$STACKS_FOLDER" -i {#} -p 1 -m 4 -M 3 -N 5 -H --deleverage
