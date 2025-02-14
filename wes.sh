#!/bin/bash
# PRE-REQUISITES FOR RUNNING THIS PIPELINE:
# Ensure a conda environement with bwa, samtools, cnvkit and other tools used are picard, gatk,
#	 annovar, snpeff.
# Issues that you might face with the pipeline: JAVA version error for gatk and picard, this
# 	pipeline is written for a specific machine optimise it according to your needs.
#	Using a conda env with java tools will throw a java path error and version error.
# Files you will need as an input for the whole SNP. indel and CNV analysis:
#	fastq files
#	hg38 reference fasta file
#	for clinvar analysis the clinvar.vcf file
#	BED files (exome capture kit used by the company for sequencing)
#	access file generated from the hg38.fa file using this command: cnvkit.py access hg38.fa -o access.hg38.bed
#	refFlat.txt file downloaded from this link specifically for the updated annotations: https://groups.google.com/a/soe.ucsc.edu/g/genome/c/SIS3Weqx9bA
#	target and antitarget files automatically generated from the autobin command, no need to worry about that.
#	------------------------------------ This is written dynamically so just list the prefixes of the fastq files in the array named c the code will do the rest for all the files -------------------------

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate exp ## change "exp" to your own environment name and this was done to fix the java path error.



############################################# gatk good practises based WES pipeline #################################################

s=("023-c" "023-p2")  ## add the prefix of all the fastq files in this array and make sure rest of the file name is same for all the files.
for n in ${s[@]};
do
#	bwa index -p ref_genome reference/hg38.fa
	echo "bwa mem running... "
	conda activate exp
	bwa mem -M reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta /home/user1/test/$n-R1.fastq.gz /home/user1/test/$n-R2.fastq.gz > variants/$n.sam -t 7
	echo "sam done"
	echo "sam to bam"
	samtools view variants/$n.sam -o variants/$n.bam
	echo "bam done"
	rm variants/$n.sam
	echo "sorting"
	samtools sort variants/$n.bam -o variants/$n.sorted.bam
	echo "sorted"
	rm variants/$n.bam
	echo "Mark and remove duplicates"
#	conda deactivate
	java -jar /home/user1/test/tools/picard.jar MarkDuplicates I=variants/$n.sorted.bam O=variants/$n.dup.bam \
		M=variants/$n.marked_metrics.txt REMOVE_DUPLICATES=true
	echo "removed duplicates"
	rm variants/$n.sorted.bam
	echo "Add read groups"
	java -jar /home/user1/test/tools/picard.jar AddOrReplaceReadGroups I=variants/$n.dup.bam O=variants/$n.new.bam \
		RGID=4 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
	echo "Added read groups"
	rm variants/$n.dup.bam
	echo "indexing"
#	conda activate exp
	samtools index variants/$n.new.bam
	echo "Haplotype Caller"
	echo "Base recalibrator"
#	Before feeding the vcf files for the snps and indels you will need to index the files using this command : gatk IndexFeatureFile -I resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf
#	conda deactivate
	gatk-4.6.1.0/gatk BaseRecalibrator -I variants/$n.new.bam -R reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
		--known-sites resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
		--known-sites resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf -O variants/$n-recal.table
	echo "Applying BQSR"
	gatk-4.6.1.0/gatk ApplyBQSR -R reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I variants/$n.new.bam \
		--bqsr-recal-file variants/$n-recal.table -O variants/$n.bqsr.bam
	echo "BQSR done"
	rm variants/$n.new.bam
	gatk-4.6.1.0/gatk --java-options "-Xmx8g" HaplotypeCaller \
		-R /home/user1/test/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -D /home/user1/test/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
		-I variants/$n.bqsr.bam -O variants/$n.vcf
	echo "Variant calling done"
	conda activate exp
#	bcftools filter -i'GQ=99 && AD>=10 && SOR<=1.5 && QD>=10 && MQ>=59 && ReadPosRankSum<=3 && FMT/DP>=30 && FMT/DP<=70' variants/$n.vcf -O v -o variants/$n.filtered.vcf

#	bcftools merge --force-samples $n-p1.filtered.vcf $n-p2.filtered.vcf > $n.merged.vcf
#	bgzip $n.merged.vcf
#	bcftools index $n.merged.vcf
#	bcftools isec -p $n-isec_files $n-c.filtred $n.merged.vcf
done

################################################################################# GVCF version ####################################################################################

#s=("61")
#for f in ${s[@]};
#do

#	echo "Combine GVCFs"
#	/home/chandrakanta/tests/tools/gatk-4.6.0.0/gatk CombineGVCFs -R reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta --variant variants/$f-c.g.vcf.gz \
#		--variant variants/$f-p1.g.vcf.gz --variant variants/$f-p2.g.vcf.gz -O variants/$f-trio.g.vcf.gz
#
#	echo "Convert GVCFs to VCFs"
#	/home/chandrakanta/tests/tools/gatk-4.6.0.0/gatk --java-options "-Xmx8g" GenotypeGVCFs -R reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
#		-V variants/$f-trio.g.vcf.gz -O variants/$f-trio.vcf.gz
#	echo "VCF generated"
#	/home/chandrakanta/tests/tools/gatk-4.6.0.0/gatk MakeSitesOnlyVcf -I variants/$f-trio.vcf.gz -O variants/$f-trio.sites.vcf.gz
#	echo "Only variants VCF generated"
#	echo "VQSR for SNP"
#	/home/chandrakanta/tests/tools/gatk-4.6.0.0/gatk VariantRecalibrator -R reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -V variants/$f-trio.sites.vcf.gz \
#               --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /home/chandrakanta/tests/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
#               --resource:omni,known=false,training=true,truth=false,prior=12.0 /home/chandrakanta/tests/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz \
#               --resource:1000G,known=false,training=true,truth=false,prior=10.0 /home/chandrakanta/tests/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf \
#               --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /home/chandrakanta/tests/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
#               -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP \
#               -O variants/$f.SNP.recal --tranches-file variants/$f.SNP.tranches
#
#	/home/chandrakanta/tests/tools/gatk-4.6.0.0/gatk ApplyVQSR -R reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -V variants/$f-trio.sites.vcf.gz --mode SNP \
#		--truth-sensitivity-filter-level 99.0 --recal-file variants/$f.SNP.recal --tranches-file variants/$f.SNP.tranches -O variants/$f.recalibrated_snps.vcf
#
#
#	echo "VQSR for indels"
#	/home/chandrakanta/tests/tools/gatk-4.6.0.0/gatk VariantRecalibrator -R reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -V variants/$f.recalibrated_snps.vcf \
#               --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /home/chandrakanta/tests/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz \
#               --resource:omni,known=false,training=true,truth=false,prior=12.0 /home/chandrakanta/tests/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz \
#               --resource:1000G,known=false,training=true,truth=false,prior=10.0 /home/chandrakanta/tests/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
#               --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /home/chandrakanta/tests/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
#               -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode INDEL \
#               -O variants/$f.INDEL.recal --tranches-file variants/$f.INDEL.tranches
#
#	/home/chandrakanta/tests/tools/gatk-4.6.0.0/gatk ApplyVQSR -R reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -V variants/$f.recalibrated_snps.vcf --mode INDEL \
#                --truth-sensitivity-filter-level 99.0 --recal-file variants/$f.INDEL.recal --tranches-file variants/$f.INDEL.tranches -O variants/$f.recalibrated_snps_indels.vcf
#
#	echo "Joint genotyping on merged GVCFs"
#	/home/chandrakanta/tests/tools/gatk-4.6.0.0/gatk CalculateGenotypePosteriors -V variants/$f.recalibrated_snps_indels.vcf \
#		-ped variants/$f.trio.ped --supporting-callsets /home/chandrakanta/tests/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz -O variants/$f.trio.vcf
#	echo "Variant filtering"
#	/home/chandrakanta/tests/tools/gatk-4.6.0.0/gatk VariantFiltration -R reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -V variants/$f.trio.vcf \
#		--filterExpression "GQ<20" --filterName "lowGQ" -O variants/$f.trio_filtered.vcf
#	echo "DeNovo variant calling"
#	/home/chandrakanta/tests/tools/gatk-4.6.0.0/gatk VariantAnnotator -R reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -V variants/$f.trio_filtered.vcf \
#		-A PossibleDeNovo -ped variants/$f.trio.ped -O variants/$f.trio_filtered.denovos.vcf
#
#	conda activate exp

################################################## commands for comparing variants with known databases and finding novelunreported mutations ################################################################

#	awk 'FNR==NR {a[$1,$2,$4,$5]; next} {if (($1,$2,$4,$5) in a) print $0}' 0002.vcf ../../resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf > matched_variants.vcf

#####################################################################################################

#	echo "annotation using annovar"
#	perl /home/chandrakanta/tests/tools/annovar/convert2annovar.pl -format vcf4 $n.vcf > $n.avinput
#	perl /home/chandrakanta/tests/tools/annovar/annotate_variation.pl -geneanno -buildver hg38 $n.avinput /home/chandrakanta/tests/tools/annovar/humandb/
#
#	perl /home/chandrakanta/tests/tools/annovar/table_annovar.pl $n.avinput /home/chandrakanta/tests/tools/annovar/humandb/ -buildver hg38 -out $n -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -csvout -polish
#	conda deactivate
#	java -Xmx8g -jar /home/chandrakanta/tests/tools/snpEff/snpEff.jar -v -stats $f.html hg38 $f.vcf > $n.ann.vcf
#	conda activate exp
#	echo "annotation done with snpeff and annovar"



################################### Pathogenicity analysis using Clinvar vcf file ##########################################3

#	echo"Pathogenicity analysis started"
#	awk 'NR == FNR {a[$1,$2,$4,$5];next} (("chr"$1,$2,$4,$5) in a)' $n.vcf /home/chandrakanta/tests/clinvar.vcf > $n.clinvar_validated_variants.vcf
#	echo "total number of clinvar variants:"
#	awk '!/^#/' $n.clinvar_validated_variants.vcf | wc -l
#	grep -w 'CLNSIG' $n.clinvar_validated_variants.vcf > $n.clnid_available.vcf
#	grep -w 'CLNSIG=Pathogenic' $n.clinid_available.vcf > $n.pathogenic.vcf
#	echo "Number of pathogenic variants:"
#	awk '!/^#/' $n.pathogenic.vcf | wc -l
#       grep -w 'CLNSIG=Likely_pathogenic' $n.clinid_available.vcf > $n.likely_pathogenic.vcf
#       echo "Number of Likely pathogenic variants:"
#       awk '!/^#/' $n.likely_pathogenic.vcf | wc -l
#       grep -w 'CLNSIG=Pathogenic/Likely_pathogenic' $n.clinid_available.vcf > $n.pathogenic_likelypathogenic.vcf
#       echo "Number of pathogenic/likely pathogenic variants:"
#       awk '!/^#/' $n.pathogenic_likelypathogenic.vcf | wc -l
#       grep -w 'CLNSIG=Uncertain_significance' $n.clinid_available.vcf > $n.uncertain_significance.vcf
#       echo "Number of uncertain significance variants:"
#       awk '!/^#/' $n.uncertain_significance.vcf | wc -l
#       grep -w 'CLNSIG=Conflicting_interpretations_of_pathogenicity' $n.clinid_available.vcf > $n.conflicting.vcf
#       echo "Number of Conflicting interpretations of pathogenicity variants:"
#       awk '!/^#/' $n.conflicting.vcf | wc -l
#	echo "Pathogenicity analysis with clinvar done"


##########################################  CNV ANALYSIS  ############################################
#	echo "cnv analysis started"
#	cnvkit.py batch $n-c.sorted.bam --normal $n-p2.sorted.bam --targets 2ndbatchchild.bed --annotate refFlat.txt --fasta ../reference/hg38.fa --access access.hg38.bed --output-reference $n-my_reference.cnn --output-dir $n-p2-results/ --diagram --scatter

# for cnvkit individual run-(not configured properly, go for the batch command, above)

#cnvkit.py autobin 101-*.bqsr.bam -t my_targets.bed -g access.hg38_new.bed --annotate refFlat.txt

#cnvkit.py coverage 101-c.bqsr.bam my_targets.target.bed -o 101-c.targetcoverage.cnn

#cnvkit.py coverage 101-c.bqsr.bam my_targets.antitarget.bed -o 101-c.antitargetcoverage.cnn

#cnvkit.py coverage 101-p2.bqsr.bam my_targets.target.bed -o 101-p2.targetcoverage.cnn

#cnvkit.py coverage 101-p2.bqsr.bam my_targets.antitarget.bed -o 101-p2.antitargetcoverage.cnn

#cnvkit.py reference 101-p2.{,anti}targetcoverage.cnn --fasta reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -o 101-my_reference.cnn

#cnvkit.py fix 101-c.targetcoverage.cnn 101-c.antitargetcoverage.cnn 101-my_reference.cnn -o 101.cnr

#cnvkit.py segment 101.cnr -o 101.cns

#cnvkit.py scatter 101.cnr -s 101.cns -o 101-scatter.pdf

#cnvkit.py diagram 101.cnr -s 101.cns -o 101-diagram.pdf


#done

echo "Done pipeline"
