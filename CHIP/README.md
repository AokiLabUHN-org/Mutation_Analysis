# Mutation_Analysis/CHIP
Description of scripts:

1. ChipByGeneAndVaf.Rmd: Analysis of CHIP in HCC cohort stratifying patients by CHIP gene and VAF of domninant clone

2. automateFiltering.R: R script to take in a filtered excel workbook (smMIP panel output) and automatically classify each mutation as a putative or non-putative CHIP driver event

3. identifyDrivers.R: R script which applies the same logic as automateFiltering.R but extends it to L-CHIP events as well. Can be applied to the output of parse_vcf_2.py following Mutect2/Funcotator analysis of any type of sequencing data (e.g., targeted panel, WES, WGS)

4. parse_vcf_2.py : Python script to parse Funcotated vcf files into small .tsv files containing all useful fields

