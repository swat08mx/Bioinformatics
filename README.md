# Deathstar.py
Deathstar.py is a tool currently in development that can be used specifically to analyze the outputs from the tool CNVkit, this tool takes the cnr files, preferentially from a cohort and processess them together to generate plots for easier visualization of alterations in the genome, one has to put their gene list as in an array and feed it to the loop, along with the cnr files, most of the other tasks are automated, except in some places where manual ovveride is needed to filter amplifications and deletions based on the z score methods and to change the thresholds as needed by the user. This tool also uses a z score based filteration method to find out bins where the deletion or amplification is more than normal in that sample. Which helps us to identify potential regions which could harbour CNVs. 
This code might need additional packages such as num2words.
```bash
pip install num2words
```
# Metagenomics.sh
Metagenomics.sh is a metagenomics pipeline made out of qiime2 which has a lot of cutomizations that can be enabled and a lot of tools to choose from, like the denoisers or phylogenetic tree creation. This pipeline is still under development, with further improvements being done to provide customized plots based out of R and more downstream analysis for better understanding of the data.
# test.nf and wes.sh
These are the same pipelines in two different implementation for the WES data processing, made according to the GATK best practices and some more additonal analysis preseent in the bash file.
