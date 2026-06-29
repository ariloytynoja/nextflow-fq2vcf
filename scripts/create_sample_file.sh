#!/bin/bash

fastq_root='s3://fastq_data'
output_root='s3://processed_ribbon'

sample_file=$1
sub_folder=$2

fastqs=$(s3cmd ls --recursive $fastq_root | grep fq.gz | grep -w -f $sample_file | awk '{print $4}')

echo samples:
for smp in $(cat $sample_file); do
	echo "  - sample_id: $smp"
        echo "    fastq:"
	echo $fastqs | xargs -n1 | grep -w $smp | sed 's/^/      - /'
	echo "    bam:  $output_root/$sub_folder/$smp/${smp}_markdup.cram"
	echo "    gvcf: $output_root/$sub_folder/$smp/${smp}.gvcf.gz"
	echo "    type: SR"
done
