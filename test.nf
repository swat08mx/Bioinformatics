nextflow.enable.dsl=2

params.reads = "/home/user1/test/023-c-R{1,2}.fastq.gz"
params.ref = "/home/user1/test/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta" 
params.outdir = "/home/user1/test/testing-result/" 
params.mills = "/home/user1/test/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf"
params.dbsnp = "/home/user1/test/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf"

workflow {
	read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
	aligned = align(read_pairs_ch)
        conversion = samtobam(aligned)
        sorted = sorting(conversion)
        duplicated = markdups(sorted)
	index = indexing(duplicated.bam)
	bqsrfile = BQSR(duplicated.bam)
	appbqsr = applybqsr(duplicated.bam, bqsrfile)
	variantcalling = varcalling(appbqsr)
}

process align {
	publishDir params.outdir, mode: 'copy'
	tag "$sample_id"

	input:
	tuple val(sample_id), path(reads)

	output:
	path "${sample_id}.sam"

	script:
	""" 
	bwa mem -R '@RG\\tID:4\\tLB:lib1\\tPL:ILLUMINA\\tPU:unit1\\tSM:20' -M ${params.ref} ${reads[0]} ${reads[1]} -o ${sample_id}.sam -t 4	
	"""
}

process samtobam {
	publishDir params.outdir, mode: 'copy' 
	tag "$processed_file" 
	
	input:
	path processed_file
	
	output:
	path "${processed_file.baseName}.bam"

	script:
	"""
	samtools view ${processed_file} -o ${processed_file.baseName}.bam
	"""
}

process sorting {
	publishDir params.outdir, mode: 'copy'
	tag "$unsorted_file"

	input:
	path unsorted_file

	output:
	path "${unsorted_file.simpleName}.sorted.bam"

	script:
	"""
	samtools sort ${unsorted_file} -o ${unsorted_file.simpleName}.sorted.bam
	"""
}

process markdups {
	publishDir params.outdir, mode: 'copy'
	tag "$sorted_bam"
	
	input:
	path sorted_bam

	output:
	path "${sorted_bam.simpleName}.marked_metrics.txt", emit: txt
	path "${sorted_bam.simpleName}.dup.bam", emit: bam
	
	script:
	"""
	java -jar /home/user1/test/tools/picard.jar MarkDuplicates I=${sorted_bam} O=${sorted_bam.simpleName}.dup.bam M=${sorted_bam.simpleName}.marked_metrics.txt REMOVE_DUPLICATES=true

	"""
}

process indexing {
        publishDir params.outdir, mode: 'copy'
        tag "$add_bam"

	input:
	path add_bam

	output:
	path "${add_bam.simpleName}.bam.bai"

	script:
	"""
	samtools index ${add_bam} -o ${add_bam.simpleName}.bam.bai
	"""
}

process BQSR {
	publishDir params.outdir, mode: 'copy' 
	tag "$sorted_bam"
	
	input:
	path sorted_bam
	
	output:
	path "${sorted_bam.simpleName}.recal.table" 
	
	script:
	"""
	/home/user1/test/gatk-4.6.1.0/gatk BaseRecalibrator -I ${sorted_bam} -R ${params.ref} --known-sites ${params.mills} --known-sites ${params.dbsnp} -O ${sorted_bam.simpleName}.recal.table

	"""
}

process applybqsr {
	publishDir params.outdir, mode: 'copy'
	tag "$sorted_bam"
	
	input:
	path sorted_bam 
	path recal_data
	
	output:
	path "${sorted_bam.simpleName}.bqsr.bam"
	
	script:
	"""
	/home/user1/test/gatk-4.6.1.0/gatk ApplyBQSR -R ${params.ref} -I ${sorted_bam} --bqsr-recal-file ${recal_data} -O ${sorted_bam.simpleName}.bqsr.bam
	"""
}

process varcalling {
	publishDir params.outdir, mode: 'copy'
	tag "$bqsr_bam"
	
	input:
	path bqsr_bam

	output:
	path "${bqsr_bam.simpleName}.vcf"

	script:
	"""
	/home/user1/test/gatk-4.6.1.0/gatk --java-options "-Xmx8g" HaplotypeCaller -R ${params.ref} -D ${params.dbsnp} -I ${bqsr_bam} -O ${bqsr_bam.simpleName}.vcf

	"""
}


