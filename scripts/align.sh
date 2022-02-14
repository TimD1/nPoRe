#!/bin/bash

reads="$1"
reference="$2"
prefix="${reads%.*}"
if ! [ -z "$3" ]; then prefix="$3"; fi

minimap2="/home/timdunn/software/minimap2/minimap2"
samtools="/home/timdunn/software/samtools-1.10/samtools"

$samtools faidx $reference

if [[ $reads == *.fasta* ]] || [[ $reads == *.fastq* ]]; then
	$minimap2 \
		-ax map-ont \
		--eqx \
		--cap-sw-mem=4g \
		-t $(nproc) \
		$reference \
		$reads |
    $samtools view \
        -@ $(nproc) \
        -b \
        -h \
        -F 2304 \
    $samtools sort \
        -@ $(nproc) |
    $samtools calmd \
        -@ $(nproc) \
        -b \
        -Q \
        - \
        $reference \
        > ${prefix}.bam

elif [[ $reads == *.bam || $reads == *.sam || $reads == "-" ]]; then
    $samtools view \
        -@ $(nproc) \
        -b \
        -h \
        -F 2304 \
        $reads |
    $samtools sort \
        -@ $(nproc) |
    $samtools calmd \
        -@ $(nproc) \
        -b \
        -Q \
        - \
        $reference \
        > ${prefix}.bam

else # force SAM
    echo "ERROR: unexpected file extension"
    exit 1
fi

$samtools index \
    -@ $(nproc) \
    ${prefix}.bam
