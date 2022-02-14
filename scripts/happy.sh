#!/bin/bash

main_dir="/x/gm24385"
guppy="guppy_5_0_6"

start_chr=20
end_chr=22

callvcfs=(
    # "$main_dir/test/$guppy/g5-ra-hap/1_calls.vcf.gz"
    # "$main_dir/test/$guppy/g5-hap/2_calls.vcf.gz"
    # "$main_dir/test/$guppy/g5-ra-hap/2_calls.vcf.gz"
    # "$main_dir/test/$guppy/g5-orig/1_calls.vcf.gz"
    "$main_dir/test/$guppy/g5-orig/1_variants.vcf.gz"
)
callvcfnames=(
    # 'g5'
    # 'g5-hap'
    # 'g5-ra-hap'
    # 'g5-orig1'
    'g5-orig2'
)

truthvcfs=(
    "$main_dir/ref/HG002_GRCh38_1_22_v4.1_draft_benchmark.vcf.gz"
    # "$main_dir/test/$guppy/g5-ra-hap/ref/1_std.vcf.gz"
)

truthvcfnames=(
    'truth41'
    # 'truth41std'
)

fullbeds=(
    "$main_dir/reference/np_0.bed"
    # "$main_dir/reference/np_1.bed"
    # "$main_dir/reference/np_2.bed"
    # "$main_dir/reference/np_3.bed"
    # "$main_dir/reference/np_4.bed"
    # "$main_dir/reference/np_5.bed"
    # "$main_dir/reference/np_6.bed"
    "$main_dir/reference/np_all.bed"
    "$main_dir/reference/all.bed"
)
bednames=(
    'np_0'
    # 'np_1'
    # 'np_2'
    # 'np_3'
    # 'np_4'
    # 'np_5'
    # 'np_6'
    'np_all'
    'all'
)

# create list of bed file names for just eval chrs
beds=()
for bed in ${bednames[@]}; do
    beds+=($main_dir/reference/${bed}_${start_chr}_${end_chr}.bed)
done

# create BED region files for just eval chrs
for i in ${!fullbeds[@]}; do
    seq $start_chr $end_chr | \
        sed 's/^[0-9]/chr&/' | \
        grep -f - ${fullbeds[i]} > \
        ${beds[i]}
done

evalbeds=(
    "$main_dir/ref/HG002_GRCh38_1_22_v4.1_draft_benchmark.bed"
    # "$main_dir/reference/all.bed"
)
evalbednames=(
    'eval41'
    # 'evalall'
)

ref_fasta="GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
export HGREF="$main_dir/reference/$ref_fasta"
source ~/software/happy/venv2/bin/activate
mkdir -p results/data/par-happy

parallel --joblog happy.log -j25 \
    "python ~/software/happy/install/bin/hap.py \
        {3} \
        {1} \
        -r $main_dir/reference/$ref_fasta \
        --engine-vcfeval-template $main_dir/reference/${ref_fasta}.sdf \
        -T {7} \
        -R {5} \
        --roc QUAL \
        --write-counts \
        --engine vcfeval \
        -o results/data/par-happy/{2}-{6}-{4}-{8}" ::: \
            ${callvcfs[@]} :::+ ${callvcfnames[@]} ::: \
            ${truthvcfs[@]} :::+ ${truthvcfnames[@]} ::: \
            ${beds[@]} :::+ ${bednames[@]} ::: \
            ${evalbeds[@]} :::+ ${evalbednames[@]}

parallel -j25 \
    "gunzip -f results/data/par-happy/{2}-{6}-{4}-{8}.*.gz" ::: \
        ${callvcfs[@]} :::+ ${callvcfnames[@]} ::: \
        ${truthvcfs[@]} :::+ ${truthvcfnames[@]} ::: \
        ${beds[@]} :::+ ${bednames[@]} ::: \
        ${evalbeds[@]} :::+ ${evalbednames[@]}

deactivate
