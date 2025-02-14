#!/bin/bash

cnvkit.py autobin 61-c.bqsr.bam -t 2ndbatchchild.bed -g access.hg38_new.bed --annotate refFlat.txt
#s=("61-c" "61-p2" "78-c" "78-p2" "80-c" "80-p2" "81-c" "81-p2" "87-c" "87-p2" "90-c" "90-p2" "92-c" "92-p2" "101-c" "101-p2" "105-c" "105-p2")
#for n in ${s[@]};
cnvkit.py coverage 61-c.bqsr.bam 2ndbatchchild.target.bed -o 61-c.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 61-c.bqsr.bam 2ndbatchchild.antitarget.bed -o 61-c.antitarget.cnn
echo "coverage2"
cnvkit.py autobin 78-c.bqsr.bam -t 2ndbatchchild.bed -g access.hg38_new.bed --annotate refFlat.txt
cnvkit.py coverage 78-c.bqsr.bam 2ndbatchchild.target.bed -o 78-c.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 78-c.bqsr.bam 2ndbatchchild.antitarget.bed -o 78-c.antitarget.cnn
echo "coverage2"
cnvkit.py autobin 80-c.bqsr.bam -t 2ndbatchchild.bed -g access.hg38_new.bed --annotate refFlat.txt
cnvkit.py coverage 80-c.bqsr.bam 2ndbatchchild.target.bed -o 80-c.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 80-c.bqsr.bam 2ndbatchchild.antitarget.bed -o 80-c.antitarget.cnn
echo "coverage2"
cnvkit.py autobin 81-c.bqsr.bam -t 2ndbatchchild.bed -g access.hg38_new.bed --annotate refFlat.txt
cnvkit.py coverage 81-c.bqsr.bam 2ndbatchchild.target.bed -o 81-c.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 81-c.bqsr.bam 2ndbatchchild.antitarget.bed -o 81-c.antitarget.cnn
echo "coverage2"
cnvkit.py autobin 87-c.bqsr.bam -t 2ndbatchchild.bed -g access.hg38_new.bed --annotate refFlat.txt
cnvkit.py coverage 87-c.bqsr.bam 2ndbatchchild.target.bed -o 87-c.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 87-c.bqsr.bam 2ndbatchchild.antitarget.bed -o 87-c.antitarget.cnn
echo "coverage2"
cnvkit.py autobin 90-c.bqsr.bam -t 2ndbatchchild.bed -g access.hg38_new.bed --annotate refFlat.txt
cnvkit.py coverage 90-c.bqsr.bam 2ndbatchchild.target.bed -o 90-c.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 90-c.bqsr.bam 2ndbatchchild.antitarget.bed -o 90-c.antitarget.cnn
echo "coverage2"
cnvkit.py autobin 92-c.bqsr.bam -t 2ndbatchchild.bed -g access.hg38_new.bed --annotate refFlat.txt
cnvkit.py coverage 92-c.bqsr.bam 2ndbatchchild.target.bed -o 92-c.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 92-c.bqsr.bam 2ndbatchchild.antitarget.bed -o 92-c.antitarget.cnn
echo "coverage2"
cnvkit.py autobin 101-c.bqsr.bam -t 2ndbatchchild.bed -g access.hg38_new.bed --annotate refFlat.txt
cnvkit.py coverage 101-c.bqsr.bam 2ndbatchchild.target.bed -o 101-c.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 101-c.bqsr.bam 2ndbatchchild.antitarget.bed -o 101-c.antitarget.cnn
echo "coverage2"
cnvkit.py autobin 105-c.bqsr.bam -t 2ndbatchchild.bed -g access.hg38_new.bed --annotate refFlat.txt
cnvkit.py coverage 105-c.bqsr.bam 2ndbatchchild.target.bed -o 105-c.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 105-c.bqsr.bam 2ndbatchchild.antitarget.bed -o 105-c.antitarget.cnn
echo "coverage2"
cnvkit.py autobin 61-p2.bqsr.bam -t 2ndbatchchild.bed -g access.hg38_new.bed --annotate refFlat.txt
cnvkit.py coverage 61-p2.bqsr.bam 2ndbatchchild.target.bed -o 61-p2.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 61-p2.bqsr.bam 2ndbatchchild.antitarget.bed -o 61-p2.antitarget.cnn
echo "coverage2"
cnvkit.py autobin *.bam -t 2ndbatchchild.bed -g access.hg38_new.bed --annotate refFlat.txt
cnvkit.py coverage 78-p2.bqsr.bam 2ndbatchchild.target.bed -o 78-p2.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 78-p2.bqsr.bam 2ndbatchchild.antitarget.bed -o 78-p2.antitarget.cnn
echo "coverage2"
cnvkit.py coverage 80-p2.bqsr.bam 2ndbatchchild.target.bed -o 80-p2.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 80-p2.bqsr.bam 2ndbatchchild.antitarget.bed -o 80-p2.antitarget.cnn
echo "coverage2"
cnvkit.py coverage 81-p2.bqsr.bam 2ndbatchchild.target.bed -o 81-p2.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 81-p2.bqsr.bam 2ndbatchchild.antitarget.bed -o 81-p2.antitarget.cnn
echo "coverage2"
cnvkit.py coverage 87-p2.bqsr.bam 2ndbatchchild.target.bed -o 87-p2.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 87-p2.bqsr.bam 2ndbatchchild.antitarget.bed -o 87-p2.antitarget.cnn
echo "coverage2"
cnvkit.py coverage 90-p2.bqsr.bam 2ndbatchchild.target.bed -o 90-p2.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 90-p2.bqsr.bam 2ndbatchchild.antitarget.bed -o 90-p2.antitarget.cnn
echo "coverage2"
cnvkit.py coverage 92-p2.bqsr.bam 2ndbatchchild.target.bed -o 92-p2.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 92-p2.bqsr.bam 2ndbatchchild.antitarget.bed -o 92-p2.antitarget.cnn
echo "coverage2"
cnvkit.py coverage 101-p2.bqsr.bam 2ndbatchchild.target.bed -o 101-p2.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 101-p2.bqsr.bam 2ndbatchchild.antitarget.bed -o 101-p2.antitarget.cnn
echo "coverage2"
cnvkit.py coverage 105-p2.bqsr.bam 2ndbatchchild.target.bed -o 105-p2.targetcoverage.cnn
echo "coverage1"
cnvkit.py coverage 105-p2.bqsr.bam 2ndbatchchild.antitarget.bed -o 105-p2.antitarget.cnn
echo "coverage2"

cnvkit.py reference *-p2.{,anti}targetcoverage.cnn --fasta resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -o myreference.cnn
echo "reference"
#for i in ${s[@]};
#cnvkit.py fix $i.targetcoverage.cnn $i.antitargetcoverage.cnn $i.reference.cnn -o $i.cnr
#cnvkit.py segment $i.cnr -o $i.cns
#cnvkit.py scatter $i.cnr -s $i.cns -o $i-scatter.pdf
#cnvkit.py diagram $i.cnr -s $i.cns -o $i-diagram.pdf
#echo "done"
