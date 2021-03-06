#!/bin/bash

echo "> zipping VCF"
bgzip -f test_std_vcf.vcf
tabix -f -p vcf test_std_vcf.vcf.gz

echo "> standardizing VCF"
python3 ../src/standardize_vcf.py \
    test_std_vcf.vcf.gz \
    test_std_ref.fasta \
    test_std \
    --stats_dir ../guppy5_stats
