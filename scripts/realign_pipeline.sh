#!/bin/bash

# flags for running pipeline sections
guppy5=true
out="g5-ra-hap"
download_reads=false
align_reads=false
copy_reads=false
cand_call_reads=false
train_clair3=false
rephase_reads=false
standard_ref=false
realign_reads=true
var_call_reads=false
start_itr=1
stop_itr=1

# region of interest
train_chrs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19"
test_chrs="chr20 chr21 chr22"
all_chrs="$train_chrs $test_chrs"

# set filepaths
if $guppy5; then
    bc="guppy_5_0_6"
else
    bc="guppy_4_0_11"
fi
ref_dir="/x/gm24385/ref"
main_dir="/x/gm24385/test/$bc"
out_dir="$main_dir/$out"
mkdir -p $out_dir

ref_bed="$ref_dir/HG002_GRCh38_1_22_v4.1_draft_benchmark.bed"
ref_vcf="$ref_dir/HG002_GRCh38_1_22_v4.1_draft_benchmark.vcf.gz"
ref="$ref_dir/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"

home_dir="/home/timdunn"
realign_dir="$home_dir/nPoRe"
sw_dir="$home_dir/software"
clair3_dir="$sw_dir/clair3-v0.1-r9-phased"

# init
R='\033[0;31m'
Y='\033[0;33m'
G='\033[0;32m'
W='\033[0m'
fail=false
source ~/miniconda3/etc/profile.d/conda.sh
conda activate &>/dev/null
echo -e "\n${G}### RUNNING PIPELINE FOR '$out' ###${W}"

# initialize region strings
train_chrs_str=""
for chr in $train_chrs; do
    train_chrs_str="$train_chrs_str,$chr"
done
train_chrs_str="${train_chrs_str:1}"
test_chrs_str=""
for chr in $test_chrs; do
    test_chrs_str="$test_chrs_str,$chr"
done
test_chrs_str="${test_chrs_str:1}"

# initialize beds
bed="$ref_dir/region.bed"
rm -f $bed
for chr in $all_chrs; do
    grep $chr $ref_bed >> $bed
done
test_bed="$ref_dir/test_region.bed"
rm -f $test_bed
for chr in $test_chrs; do
    grep $chr $ref_bed >> $test_bed
done



function download() {

    echo -e "\n${G}[downloading fastqs]${W}"
    ids="20201026_1644_2-E5-H5_PAG07162_d7f262d5"
    prefix="s3://ont-open-data/gm24385_2020.11/analysis/r9.4.1"
    g4="guppy_v4.0.11_r9.4.1_hac_prom"
    g5="guppy_v5.0.6_r9.4.1_sup_prom"
    sids=""
    for id in $ids; do
        sid=${id: -17:8}
        sids="$sids $sid"

        if $guppy5; then
            for chr in $all_chrs; do
                echo -e "    ${G}[downloading $bc $chr $sid fastqs]${@}"
                time aws s3 --no-sign-request cp \
                    $prefix/$id/$g4/align_unfiltered/${chr}/$g5/basecalls.fastq.gz \
                    $main_dir/fastq/${chr}-${sid}.fastq.gz
            done
        else
            echo -e "    ${G}[downloading $bc $sid fastqs]${@}"
            time aws s3 --no-sign-request cp \
                $prefix/$id/$g4/basecalls.fastq.gz \
                $main_dir/fastq/${sid}.fastq.gz
        fi

        echo -e "${G1}done!${G2}\n"
    done
    sids="${sids:1}"                                                                 

    echo -e "\n${G}[merging fastqs]${W}"
    rm -f $main_dir/fastq/all.fastq
    for sid in $sids; do
        if $guppy5; then
            for chr in $all_chrs; do
                gunzip -f $main_dir/fastq/${chr}-${sid}.fastq.gz
                cat $main_dir/fastq/${chr}-${sid}.fastq >> \
                    $main_dir/fastq/all.fastq
            done
        else
            gunzip -f $main_dir/fastq/${sid}.fastq.gz
            cat $main_dir/fastq/${sid}.fastq >> \
                $main_dir/fastq/all.fastq
        fi
    done

    for sid in $sids; do
        if $guppy5; then
            for chr in $all_chrs; do
                rm -f $main_dir/fastq/${chr}-${sid}.fastq
            done
        else
            rm -f $main_dir/fastq/${sid}.fastq
        fi
    done

}


function align() {
    echo -e "\n${G}[aligning reads]${W}"
    mkdir -p $main_dir/bam

    time ./align \
        $main_dir/fastq/all.fastq \
        $ref \
        $main_dir/bam/all || fail=true
        if $fail; then echo -e "${R}failed$W"; exit 1; fi

    rm -f $main_dir/bam/all.sam
    rm -rf $main_dir/fastq
    ln -sf $main_dir/bam/all.bam $out_dir/0_reads.bam
    ln -sf $main_dir/bam/all.bam.bai $out_dir/0_reads.bam.bai
}



function std_ref() {
    itr=$1

    echo -e "\n${G}[creating standard reference]${W}"

    mkdir -p $out_dir/ref
    if [[ $itr == 0 ]]; then
        if [ -f $out_dir/ref/0_std.vcf.gz ]; then
            echo -e "\n\t${G}[std ref: exists, skipping]${W}"
        else
            echo -e "\n\t${G}[std ref: unphasing]${W}"
            whatshap unphase \
                $ref_vcf \
                > $out_dir/ref/0_std.vcf
            echo -e "\n\t${G}[std ref: zipping vcf]${W}"
            time bgzip -f $out_dir/ref/0_std.vcf || fail=true
            if $fail; then echo -e "${R}failed$W"; exit 1; fi
            echo -e "\n\t${G}[std ref: indexing]${W}"
            time tabix -f -p vcf $out_dir/ref/0_std.vcf.gz || fail=true
            if $fail; then echo -e "${R}failed$W"; exit 1; fi
        fi
        return
    fi

    echo -e "\n\t${G}[std ref: rephasing]${W}"
    time whatshap phase \
        $out_dir/ref/$((itr-1))_std.vcf.gz \
        $out_dir/${itr}_phased.bam \
        --output $out_dir/ref/${itr}_phased.vcf.gz \
        --reference $ref \
        --ignore-read-groups \
        --indels || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    echo -e "\n\t${G}[std ref: indexing vcf]${W}"
    time tabix -f -p vcf $out_dir/ref/${itr}_phased.vcf.gz || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    if $guppy5; then
        stats_dir=$realign_dir/guppy5_stats
    else
        stats_dir=$realign_dir/stats
    fi

    echo -e "\n\t${G}[std ref: standardizing]${W}"
    mkdir -p $out_dir/ref
    time python3 $realign_dir/src/standardize_vcf.py \
        $out_dir/ref/${itr}_phased.vcf.gz \
        $ref \
        $out_dir/ref/${itr}_std \
        --contigs "$train_chrs_str,$test_chrs_str" \
        --stats_dir $stats_dir || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

}



function cand_call() {
    itr=$1

    source ./clair3/venv/bin/activate

    # optionally train new model
    echo -e "\n${G}[retraining candidate caller]${W}"
    if $train_clair3; then

        echo -e "\n\t${G}[retraining: generating tensors]${W}"
        time ./clair3/generate_pileup_tensors.sh \
            $out_dir/ref/${itr}_std.vcf.gz \
            $out_dir/${itr}_reads.bam \
            $ref \
            $ref_bed \
            $train_chrs_str \
            $out_dir/${itr}_clair3 \
            $itr || fail=true
        if $fail; then echo -e "${R}failed$W"; exit 1; fi

        echo -e "\n\t${G}[retraining: training model]${W}"
        time ./clair3/train_pileup_model.sh \
            $out_dir/${itr}_clair3 \
            $itr || fail=true
        if $fail; then echo -e "${R}failed$W"; exit 1; fi

    fi

    echo -e "\n${G}[calling candidates]${W}"
    if [[ $itr == 0 ]]; then
        time $clair3_dir/run_clair3.sh \
            --bam_fn="$out_dir/${itr}_reads.bam" \
            --ref_fn="$ref" \
            --bed_fn="$bed" \
            --threads=`nproc` \
            --platform="ont" \
            --pileup_only \
            --model_path="$out_dir/${itr}_clair3/train_pileup" \
            --output="$out_dir/tmp" || fail=true
    else
        time $clair3_dir/run_clair3.sh \
            --bam_fn="$out_dir/${itr}_reads.bam" \
            --ref_fn="$ref" \
            --bed_fn="$bed" \
            --threads=`nproc` \
            --platform="ont" \
            --pileup_only \
            --haplotypes \
            --model_path="$out_dir/${itr}_clair3/train_pileup" \
            --output="$out_dir/tmp" || fail=true
    fi
    mv $out_dir/tmp/pileup.vcf.gz $out_dir/$((itr+1))_candidates.vcf.gz
    mv $out_dir/tmp/pileup.vcf.gz.tbi $out_dir/$((itr+1))_candidates.vcf.gz.tbi
    rm -rf $out_dir/tmp
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    echo -e "\n${G}[filtering candidates]${W}"
    time bcftools filter $out_dir/$((itr+1))_candidates.vcf.gz \
        --exclude 'GT="0/0"' \
        --output $out_dir/$((itr+1))_allcalls.vcf || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    echo -e "\n\t${G}[filtering candidates: zipping vcf]${W}"
    time bgzip -f $out_dir/$((itr+1))_allcalls.vcf || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    echo -e "\n\t${G}[filtering candidates: indexing vcf]${W}"
    time tabix -f -p vcf $out_dir/$((itr+1))_allcalls.vcf.gz || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    deactivate

    echo -e "\n\t${G}[filtering candidates: removing overlaps]${W}"
    time python3 $realign_dir/src/filter.py \
        $out_dir/$((itr+1))_allcalls.vcf.gz \
        $out_dir/$((itr+1))_calls.vcf.gz || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    echo -e "\n\t${G}[filtering candidates: indexing vcf]${W}"
    time tabix -f -p vcf $out_dir/$((itr+1))_calls.vcf.gz || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

}



function rephase() {
    itr=$1

    echo -e "\n${G}[rephasing calls]${W}"
    time whatshap phase \
        $out_dir/${itr}_calls.vcf.gz \
        $out_dir/$((itr-1))_reads.bam \
        --output $out_dir/${itr}_phased.vcf.gz \
        --reference $ref \
        --ignore-read-groups \
        --indels || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    echo -e "\n\t${G}[rephasing calls: indexing vcf]${W}"
    time tabix -f -p vcf $out_dir/${itr}_phased.vcf.gz || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    echo -e "\n${G}[rephasing reads]${W}"
    time whatshap haplotag \
        --output $out_dir/${itr}_phased.bam \
        --reference $ref \
        --ignore-read-groups \
        $out_dir/${itr}_phased.vcf.gz \
        $out_dir/$((itr-1))_reads.bam || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    echo -e "\n\t${G}[rephasing reads: indexing bam]${W}"
    time ./align \
        $out_dir/${itr}_phased.bam \
        $ref || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi
}



function realign() {
    itr=$1

    echo -e "\n${Y}[REMINDER: rebuild nPoRe]${W}"

    if $guppy5; then
        stats_dir=$realign_dir/guppy5_stats
    else
        stats_dir=$realign_dir/stats
    fi

    echo -e "\n${G}[realigning reads]${W}"
    time python3 $realign_dir/src/realign.py \
        $out_dir/${itr}_phased.bam \
        $ref \
        $out_dir/${itr}_reads \
        --plot \
        --contigs "${train_chrs_str},${test_chrs_str}" \
        --stats_dir $stats_dir || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    # time ./align \
    #     $out_dir/${itr}_reads.sam \
    #     $ref || fail=true
    # # rm $out_dir/${itr}_reads.sam
    # if $fail; then echo -e "${R}failed$W"; exit 1; fi

}



function var_call() {
    itr=$1

    source ./clair3/venv/bin/activate

    # optionally train new model
    if $train_clair3; then
        echo -e "\n${G}[retraining variant caller]${W}"

        echo -e "\n\t${G}[retraining: generating tensors]${W}"
        time ./clair3/generate_full_tensors.sh \
            $out_dir/ref/${itr}_std.vcf.gz \
            $out_dir/$((itr+1))_reads.bam \
            $ref \
            $ref_bed \
            $chr \
            $chr \
            $train_beg \
            $train_end \
            $out_dir/${itr}_clair3_full \
            ${itr} || fail=true
        if $fail; then echo -e "${R}failed$W"; exit 1; fi

        echo -e "\n\t${G}[retraining: training model]${W}"
        time ./clair3/train_full_model.sh \
            $out_dir/${itr}_clair3_full \
            $itr || fail=true
        if $fail; then echo -e "${R}failed$W"; exit 1; fi

    fi

    echo -e "\n${G}[calling variants]${W}"
    time ./clair3/run_full_model.sh \
       $out_dir \
       $itr \
       $ref \
       $chr
    mv $out_dir/merge_output.vcf.gz $out_dir/$((itr+1))_variant_cands.vcf.gz
    mv $out_dir/merge_output.vcf.gz.tbi $out_dir/$((itr+1))_variant_cands.vcf.gz.tbi
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    echo -e "\n${G}[filtering variants]${W}"
    time bcftools filter $out_dir/$((itr+1))_variant_cands.vcf.gz \
        --exclude 'GT="0/0"' \
        --output $out_dir/$((itr+1))_variants.vcf || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    echo -e "\n\t${G}[filtering variants: zipping vcf]${W}"
    time bgzip -f $out_dir/$((itr+1))_variants.vcf || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    echo -e "\n\t${G}[filtering variants: indexing vcf]${W}"
    time tabix -f -p vcf $out_dir/$((itr+1))_variants.vcf.gz || fail=true
    if $fail; then echo -e "${R}failed$W"; exit 1; fi

    deactivate

}



function main() {
    mkdir -p $out_dir
    if $download_reads; then download; fi
    if $align_reads; then align; fi
    if $copy_reads; then
        ln -sf $main_dir/bam/all.bam $out_dir/0_reads.bam
        ln -sf $main_dir/bam/all.bam.bai $out_dir/0_reads.bam.bai
    fi

    # run full candidate calling and realignment pipeline
    for iter in `seq $start_itr $stop_itr`; do

        echo -e "\n${G}### ITERATION $iter ###${W}"
        if [ $iter == 0 ]; then
            if $standard_ref; then std_ref 0; fi
        else
            if $rephase_reads; then rephase ${iter}; fi
            if $standard_ref; then std_ref ${iter}; fi
            if $realign_reads; then realign ${iter}; fi
        fi
        if $cand_call_reads; then cand_call $iter; fi
    done

    if $var_call_reads; then 
        if $rephase_reads; then rephase $((iter+1)); fi
        if $standard_ref && $train_clair3; then std_ref $((iter+1)); fi
        if $realign_reads; then realign $((iter+1)); fi
        var_call $stop_itr; 
    fi
}
main
