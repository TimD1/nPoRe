#!/bin/bash

BAM="reads.bam"
REF="ref.fasta"

# run nPoRe
time python3 ../src/realign.py \
    $BAM \
    $REF \
    out \
    --stats_dir ../guppy5_stats

# update MD tags
../scripts/align.sh out.sam $REF
