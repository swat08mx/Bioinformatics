#!/bin/bash

s=("013" "022" "030" "036" "007" "044" "034" "023")  ## add the prefix of all the fastq files in this array and make sure rest of the file name is same for all the files.
for n in ${s[@]};
do
        samtools mpileup -B -q 1 -f /home/user1/test/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta /home/user1/test/variants/$n-p2.bqsr.bam /home/user1/test/variants/$n-p1.bqsr.bam /home/user1/test/variants/$n-c.bqsr.bam > $n.trio.mpileup
        varscan trio $n.trio.mpileup $n-trio.mpileup.output --min-coverage 10 --min-var-freq 0.20 --p-value 0.05 -adj-var-freq 0.05 -adj-p-value 0.15
        rm $n.trio.mpileup
done
