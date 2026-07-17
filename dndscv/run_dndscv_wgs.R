#code to determine genes/loci under selection in CPTAC and TCGA cohorts
#split cohorts by CHIP and/or CH-TUM status
first = F  #set to F after first pass to simply load processed object rather than re-processing large vcf file
if (first) {

	library(data.table)
	setDTthreads(4)
	bigvcf = fread("/cluster/projects/vannergroup/Stueckmann/out/dndscv/parsed_vcfs/GDC_VarScan2.merged.parsed.tsv")

	`%notin%` <- Negate(`%in%`)
	#intronic and intergenic mutations should not be considered by dndscv
	#removal will reduce runtime and storage of dataset
	cleaned_vcf = bigvcf[which(bigvcf$consequence %notin% c("intronic", "intergenic")),]
	saveRDS(cleaned_vcf, "/cluster/projects/vannergroup/Stueckmann/out/dndscv/all_coding_variants_merged.rds")

	#filter down large vcf file to only CPTAC cohort
	ix = str_which(cleaned_vcf$sample_id, "CPTAC")
	cptac_vcf = cleaned_vcf[ix,]
	mutations_cptac = data.frame(sampleID = cptac_vcf$sample_id, chr = gsub("chr", "", cptac_vcf$chr), pos = cptac_vcf$pos, ref = cptac_vcf$ref, alt = cptac_vcf$alt)
	saveRDS(mutations_cptac, "/cluster/projects/vannergroup/Stueckmann/out/dndscv/CPTAC_coding_variants_formatted.rds")

} else {
mutations_cptac = readRDS("/cluster/projects/vannergroup/Stueckmann/out/dndscv/CPTAC_coding_variants_formatted.rds")
}
#load required packages and metadata
library("dndscv")
library("stringr")

#create function to split case names from sample sheet
ssfunc = function(x) {
	return(strsplit(x, ",", fixed = T)[[1]][1])
}

cptac = readRDS("/cluster/projects/vannergroup/Stueckmann/input/chTum/cptac_summarized.rds")
manifest = data.table::fread("/cluster/projects/vannergroup/Stueckmann/input/chTum/gdc_manifest.2026-06-26.154622.txt")
sample_sheet = data.table::fread("/cluster/projects/vannergroup/Stueckmann/input/chTum/gdc_sample_sheet.2026-07-15.tsv")

#sample ID used to map between sample sheet and mutation files
sampleID_trimmed = gsub(".wgs.VarScan2.somatic_annotation.vcf.gz", "", sample_sheet$`File Name`, fixed = T)
caseNames = unlist(lapply(sample_sheet$`Case ID`, ssfunc))

#split sample ID list into CH+, CH-TUM, and CH- lists
chTum_cases = cptac$sample_id[which(cptac$chip_tumour == "chip_tum")]
chPb_cases = cptac$sample_id[which(cptac$chip_tumour == "no_chip_tum")]
control_cases = cptac$sample_id[which(cptac$chip_tumour == "no_chip")]

chTum_samples = sampleID_trimmed[which(caseNames %in% chTum_cases)]
#a few IDs are not present in Marco's files (no blood?) - they have been filtered out
#e.g., only 62/85 chTum cases were found in the sample sheet
chPb_samples = sampleID_trimmed[which(caseNames %in% chPb_cases)]
control_samples = sampleID_trimmed[which(caseNames %in% control_cases)]


#generate mutation files for CPTAC samples with and with CH

#dnds stats on entire CPTAC cohort (start with random sample first 100k mutations)
dndsout = dndscv(mutations_cptac, refdb = "hg38", max_muts_per_gene_per_sample = Inf, outmats = T)
saveRDS(dndsout, "/cluster/projects/vannergroup/Stueckmann/out/dndscv/CPTAC_dndsout.rds")

ix_chTum = which(mutations_cptac$sampleID %in% chTum_samples)
ix_chPb = which(mutations_cptac$sampleID %in% chPb_samples)
ix_control = which(mutations_cptac$sampleID %in% control_samples)


#all CPTAC samples with CH (both PB-only and CH-TUM)
dndsout = dndscv(mutations_cptac[c(ix_chPb, ix_chTum),], refdb = "hg38", max_muts_per_gene_per_sample = Inf, outmats = T)
saveRDS(dndsout, "/cluster/projects/vannergroup/Stueckmann/out/dndscv/CPTAC_allCH_dndsout.rds")

#all CPTAC samples with CH-TUM
dndsout = dndscv(mutations_cptac[ix_chTum,], refdb = "hg38", max_muts_per_gene_per_sample = Inf, outmats = T)
saveRDS(dndsout, "/cluster/projects/vannergroup/Stueckmann/out/dndscv/CPTAC_CH_TUM_dndsout.rds")

#all CPTAC samples without CH
dndsout = dndscv(mutations_cptac[ix_control,], refdb = "hg38", max_muts_per_gene_per_sample = Inf, outmats = T)
saveRDS(dndsout, "/cluster/projects/vannergroup/Stueckmann/out/dndscv/CPTAC_noCH_controls.rds")
