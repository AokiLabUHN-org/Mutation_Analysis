#code to determine genes/loci under selection in CPTAC and TCGA cohorts
#split cohorts by CHIP and/or CH-TUM status
library(data.table)
bigvcf = fread("/cluster/projects/vannergroup/Stueckmann/out/dndscv/parsed_mafs/CPTAC_MAF.merged.parsed.tsv")	

mutations_cptac = data.frame(sampleID = bigvcf$sample_id, chr = gsub("chr", "", bigvcf$chr), pos = bigvcf$pos, ref = bigvcf$ref, alt = bigvcf$alt)

saveRDS(mutations_cptac, "/cluster/projects/vannergroup/Stueckmann/out/dndscv/CPTAC_WES_variants_formatted.rds")

#load required packages and metadata
library("dndscv")
library("stringr")

#create function to split case names from sample sheet
ssfunc = function(x) {
	return(strsplit(x, ",", fixed = T)[[1]][1])
}

cptac = readRDS("/cluster/projects/vannergroup/Stueckmann/input/chTum/cptac_summarized.rds")
#load sample sheet and match on sample name
sample_sheet = fread("/cluster/projects/vannergroup/Stueckmann/out/CHIP/raw_mafs/inputs/gdc_sample_sheet.2026-05-29.tsv")
length(which(mutations_cptac$sampleID %in% sample_sheet$`File ID`))

#split sample ID list into CH+, CH-TUM, and CH- lists
chTum_cases = cptac$sample_id[which(cptac$chip_tumour == "chip_tum")]
chPb_cases = cptac$sample_id[which(cptac$chip_tumour == "no_chip_tum")]
control_cases = cptac$sample_id[which(cptac$chip_tumour == "no_chip")]

chTum_samples = sampleID_trimmed[which(caseNames %in% chTum_cases)]
chPb_samples = sampleID_trimmed[which(caseNames %in% chPb_cases)]
control_samples = sampleID_trimmed[which(caseNames %in% control_cases)]

caseNames = unlist(lapply(sample_sheet$`Case ID`, ssfunc))


chTum_samples = sample_sheet$`File ID`[which(caseNames %in% chTum_cases)]
#a few IDs are not present in Marco's files (no blood?) - they have been filtered out
#e.g., only 62/85 chTum cases were found in the sample sheet
chPb_samples = sample_sheet$`File ID`[which(caseNames %in% chPb_cases)]
control_samples = sample_sheet$`File ID`[which(caseNames %in% control_cases)]


#dnds stats on entire CPTAC cohort (start with random sample first 100k mutations)
dndsout = dndscv(mutations_cptac, refdb = "hg38", max_muts_per_gene_per_sample = Inf, outmats = T)
saveRDS(dndsout, "/cluster/projects/vannergroup/Stueckmann/out/dndscv/CPTAC_WES_allSamples_dndsout.rds")

ix_chTum = which(mutations_cptac$sampleID %in% chTum_samples)
ix_chPb = which(mutations_cptac$sampleID %in% chPb_samples)
ix_control = which(mutations_cptac$sampleID %in% control_samples)

#all CPTAC samples with CH (both PB-only and CH-TUM)
dndsout = dndscv(mutations_cptac[c(ix_chPb, ix_chTum),], refdb = "hg38", max_muts_per_gene_per_sample = Inf, outmats = T)
saveRDS(dndsout, "/cluster/projects/vannergroup/Stueckmann/out/dndscv/CPTAC_WES_allCH_dndsout.rds")

#all CPTAC samples with CH-TUM
dndsout = dndscv(mutations_cptac[ix_chTum,], refdb = "hg38", max_muts_per_gene_per_sample = Inf, outmats = T)
saveRDS(dndsout, "/cluster/projects/vannergroup/Stueckmann/out/dndscv/CPTAC_WES_CH_TUM_dndsout.rds")

#all CPTAC samples without CH
dndsout = dndscv(mutations_cptac[ix_control,], refdb = "hg38", max_muts_per_gene_per_sample = Inf, outmats = T)
saveRDS(dndsout, "/cluster/projects/vannergroup/Stueckmann/out/dndscv/CPTAC_WES_noCH_controls.rds")
