#!/bin/bash
## Pre-requisites: create your own metadata file and set your own trunc and trim parameters

echo "Pipeline started"

#s=("SRR24912546" "SRR24912547" "SRR24912548" "SRR24912549" "SRR24912550" "SRR24912551" "SRR24912552" "SRR24912553" "SRR24912554" "SRR24912555" "SRR24912556" "SRR24912557" "SRR24912558" "SRR24912559" "SRR24912560")

#for n in ${s[@]};
#do
	#fastqc $n-{1,2}.fastq -o qc/
	#fastp -i $n-1.fastq.gz -I $n-2.fastq.gz -o $n-trimmed-1.fastq.gz -O $n-trimmed-2.fastq.gz --detect_adapter_for_pe -h report.html
	#echo $n
#done
### create the manifest file using this code, edit according to file
echo "sample-id,absolute-filepath,direction" > manifest.csv
for i in *-1* ; do echo "${i/-1.fastq.gz},$PWD/$i,forward"; done >> manifest.csv
for i in *-2* ; do echo "${i/-2.fastq.gz},$PWD/$i,reverse"; done >> manifest.csv

## qiime tools import for importing the data into qiime

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.csv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33
echo "Import done"

## summarize the demultiplexed file
## check if your file is multiplxed or not

qiime demux summarize \
    --i-data paired-end-demux.qza \
    --o-visualization demux.qzv

## filter upon quality scores

#qiime quality-filter q-score \
#    --i-demux paired-end-demux.qza \
#    --o-filtered-sequences demux-filtered.qza \
#    --o-filter-stats demux-filter-stats.qza
#echo "Quality filtering done"

## deblur method to denoise the samples, set p-trim-length upon looking
## at the quality metrics

#qiime deblur denoise-16S \
#    --i-demultiplexed-seqs demux-filtered.qza \
#    --p-trim-length 150 \
#    --p-sample-stats \
#    --p-jobs-to-start 5 \
#    --o-stats deblur-stats.qza \
#    --o-representative-sequences rep-seqs-deblur.qza \
#    --o-table table-deblur.qza
#echo "Deblur denoising done"

## stats summary after denoising the samples

#qiime deblur visualize-stats \
#    --i-deblur-stats deblur-stats.qza \
#    --o-visualization deblur-stats.qzv

## for shorter runtimes use deblur which uses pre computed error profiles and for 
##     creating your own error profiles use DADA2, which comes with significantly long runtimes.

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 120 \
  --p-trunc-len-r 120 \
  --o-table table.qza \
  --o-representative-sequences rep_seqs.qza \
  --o-denoising-stats stats-dada2.qza

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv


## vizualize the representative sequences

qiime feature-table tabulate-seqs \
    --i-data rep-seqs.qza \
    --o-visualization rep-seqs-deblur.qzv

## vizualize the feature table along with meta data

qiime feature-table summarize \
    --i-table table-deblur.qza \
    --m-sample-metadata-file metadata_new.tsv \
    --o-visualization table-deblur.qzv

## phylogenetic tree creation for diversity analysis

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-deblur.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
echo "Phylogenetic tree created"

## reference based phylogenetics tree creation which takes more time than denovo

#wget -O "sepp-refs-gg-13-8.qza" \
#    "https://data.qiime2.org/2019.10/common/sepp-refs-gg-13-8.qza"

#qiime fragment-insertion sepp \
#    --i-representative-sequences rep-seqs-deblur.qza \
#    --i-reference-database sepp-refs-gg-13-8.qza \
#    --p-threads 4 \
#    --o-tree insertion-tree.qza \
#    --o-placements insertion-placements.qza
## Once tree creation is done, filtration should be done to have the sequences which were only
## theref or our dataset.

#qiime fragment-insertion filter-features \
#    --i-table table-deblur.qza \
#    --i-tree insertion-tree.qza \
#    --o-filtered-table filtered-table-deblur.qza \
#    --o-removed-table removed-table.qza



## core metrics file generation for diversity analysis

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-deblur.qza \
  --p-sampling-depth 1100 \
  --m-metadata-file metadata_new.tsv \
  --output-dir core-metrics-results
echo "Diversity core metrics done"

## alpha diversity

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata_new.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata_new.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv
echo "Alpha diversity done"

## beta diversity

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata_new.tsv \
  --m-metadata-column treatment \
  --o-visualization core-metrics-results/unweighted-unifrac-diseased-vs-healthy-significance.qzv \
  --p-pairwise
echo "Beta diversity done"

## emperor plot

#qiime emperor plot \
#  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
#  --m-metadata-file sample-metadata.tsv \
#  --p-custom-axes days-since-experiment-start \
#  --o-visualization core-metrics-results/unweighted-unifrac-emperor-days-since-experiment-start.qzv

#qiime emperor plot \
#  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
#  --m-metadata-file sample-metadata.tsv \
#  --p-custom-axes days-since-experiment-start \
#  --o-visualization core-metrics-results/bray-curtis-emperor-days-since-experiment-start.qzv

## alpha refraction plotting at multiple sample depths

qiime diversity alpha-rarefaction \
  --i-table table-deblur.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 44436 \
  --m-metadata-file metadata_new.tsv \
  --o-visualization alpha-rarefaction.qzv
echo "Alpha rarefraction done"

## Taxonomic analysis
## download the classifier

echo "Starting taxonomic analysis"

wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza"

## the classifier mostly work you will have to retrain it, following the code
echo "Retraining the classifier"
wget \
  -O "85_otus.fasta" \
  "https://data.qiime2.org/2024.10/tutorials/training-feature-classifiers/85_otus.fasta"

wget \
  -O "85_otu_taxonomy.txt" \
  "https://data.qiime2.org/2024.10/tutorials/training-feature-classifiers/85_otu_taxonomy.txt"

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 85_otus.fasta \
  --output-path 85_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path 85_otu_taxonomy.txt \
  --output-path ref-taxonomy.qza

## set the truc lens and min lens according to your needs

qiime feature-classifier extract-reads \
  --i-sequences 85_otus.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 120 \
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza

## classifier is trained on reference data locally, use the newly trained classifier

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs-deblur.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
echo "Taxonomic classification done"

## view the taxa composition in barplot

qiime taxa barplot \
  --i-table table-deblur.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata_new.tsv \
  --o-visualization taxa-bar-plots.qzv

## Volatility testing

#qiime longitudinal volatility \
#    --m-metadata-file metadata.tsv \
#    --m-metadata-file child-core-metrics-results/shannon_vector.qza \
#    --p-default-metric shannon \
#    --p-default-group-column delivery \
#    --p-state-column month \
#    --p-individual-id-column host_subject_id \
#    --o-visualization shannon-volatility.qzv


## Differential abundance testing with ANCOM-BC
echo "ANCOM-BC starting"

qiime composition add-pseudocount \
    --i-table table-deblur.qza \
    --o-composition-table comp-table-C6.qza

qiime composition ancom \
    --i-table comp-table-C6.qza \
    --m-metadata-file metadata_new.tsv \
    --m-metadata-column treatment \
    --o-visualization ancom-C6-delivery.qzv

echo "Done"

