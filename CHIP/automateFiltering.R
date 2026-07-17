#automate filtering and putative driver prediction for smMIP ouput from Sagi Abelson's ML caller
library(biomaRt)
library(stringr)
library(openxlsx)
library(data.table)
library(dplyr) #required for easy reordering of sample ID into earlier column

#pull exon-specific data from biomart - run once, load from offline after
first_run = F
if (first_run) {
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  #pull tables of ASXL1/2 transcript and exons 
  exon_names = getBM(
    attributes = c("hgnc_symbol", "ensembl_exon_id", "ensembl_transcript_id", "transcript_is_canonical"),#, "exon_chrom_start", "exon_chrom_end"),
    filters = "hgnc_symbol",
    values = c("ASXL1", "ASXL2", "PPM1D"),
    mart = ensembl
  )
  ASXL1_canonical_transcript = unique(exon_names$ensembl_transcript_id[which(exon_names$hgnc_symbol == "ASXL1" & 
                                                                               exon_names$transcript_is_canonical == 1)])
  ASXL2_canonical_transcript = unique(exon_names$ensembl_transcript_id[which(exon_names$hgnc_symbol == "ASXL2" & 
                                                                               exon_names$transcript_is_canonical == 1)])
  PPM1D_canonical_transcript = unique(exon_names$ensembl_transcript_id[which(exon_names$hgnc_symbol == "PPM1D" & 
                                                                               exon_names$transcript_is_canonical == 1)])
  
  ASXL1_exons = getBM(
    attributes = c("ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"),
    filters = "ensembl_exon_id",
    values = exon_names$ensembl_exon_id[which(exon_names$ensembl_transcript_id == ASXL1_canonical_transcript)],
    mart = ensembl
  )
  #reorder ASXL1 exons by genomic coordinates and assign generic names
  ASXL1_exons = ASXL1_exons[order(ASXL1_exons$exon_chrom_start),]
  rownames(ASXL1_exons) <- paste0("exon_", 1:nrow(ASXL1_exons))
  
  ASXL2_exons = getBM(
    attributes = c("ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"),
    filters = "ensembl_exon_id",
    values = exon_names$ensembl_exon_id[which(exon_names$ensembl_transcript_id == ASXL2_canonical_transcript)],
    mart = ensembl
  )
  #reorder ASXL2 exons by genomic coordinates and assign generic names
  ASXL2_exons = ASXL2_exons[order(ASXL2_exons$exon_chrom_start),]
  rownames(ASXL2_exons) <- paste0("exon_", 1:nrow(ASXL2_exons))
  
  PPM1D_exons = getBM(
    attributes = c("ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"),
    filters = "ensembl_exon_id",
    values = exon_names$ensembl_exon_id[which(exon_names$ensembl_transcript_id == PPM1D_canonical_transcript)],
    mart = ensembl
  )
  #reorder PPM1D exons by genomic coordinates and assign generic names
  PPM1D_exons = PPM1D_exons[order(PPM1D_exons$exon_chrom_start),]
  rownames(PPM1D_exons) <- paste0("exon_", 1:nrow(PPM1D_exons))
  
  #store output for offline work
  dir.create("~/Documents/CHIP/")
  saveRDS(ASXL1_exons, "~/Documents/CHIP/ASXL1_exons.rds")
  saveRDS(ASXL2_exons, "~/Documents/CHIP/ASXL2_exons.rds")
  saveRDS(PPM1D_exons, "~/Documents/CHIP/PPM1D_exons.rds")
} else {
  #load already generated exon info
  ASXL1_exons = readRDS("~/Documents/CHIP/ASXL1_exons.rds")
  ASXL2_exons = readRDS("~/Documents/CHIP/ASXL2_exons.rds")
  PPM1D_exons = readRDS("~/Documents/CHIP/PPM1D_exons.rds")
}

#create notin operator
`%notin%` <- Negate(`%in%`)

#create function to split read support string
#returns a vector of four integers corresponding to the four SSCS read support values 
ss_cov_string <- function(X) {
  split1 = strsplit(X, "<R>", T)
  return(as.numeric(c(
    strsplit(split1[[1]][1], "<S>", T)[[1]][1],
    strsplit(split1[[1]][1], "<S>", T)[[1]][2],
    strsplit(split1[[1]][2], "<S>", T)[[1]][1],
    strsplit(split1[[1]][2], "<S>", T)[[1]][2]
  )))
}

#create function to split sample ID string; take left half (library 1 or 2) instead of both
ss_sample_string <- function(X) {
  split1 = strsplit(X, " <R> ", T)
  return(split1[[1]][1])
}

#create function to find minimum read support in cov matrix and return a pass or fail based on threshold
#returns either 'PASS' or 'FAIL' depending on if ALL biological and probe replicates for a site have >= threshold support
#default threshold is 2 (more than 1 supporting read across all replicates)
min_cov_matrix <- function(X, threshold = 2) {
  if( str_detect(X, "NaN") ) {
    return('FAIL')
  } else if ( str_detect(X, "<M>") ) {
    left_str = trimws(strsplit(X, "<M>", T)[[1]][1])
    right_str = trimws(strsplit(X, "<M>", T)[[1]][2])
    if ( min(c(
      ss_cov_string(left_str),
      ss_cov_string(right_str))) >= threshold ) {
      return("PASS")
    } else {
      return("FAIL")
    }
  } else {
    if ( min(ss_cov_string(X)) >= threshold ) {
      return("PASS")
    } else {
      return("FAIL")
    }
  }
}

#create function to convert three letter amino acid codes in Sagi's output to one letter codes in rule set
#returns a substituted string
aa_conversion = function(X) {
  mappingVect = c(A = 'Ala', R = 'Arg', N = 'Asn', D = 'Asp',
                  C = 'Cys', E = 'Glu', Q = 'Gln', G = 'Gly',
                  H = 'His', I = 'Ile', L = 'Leu', K = 'Lys',
                  M = 'Met', F = 'Phe', P = 'Pro', S = 'Ser',
                  T = 'Thr', W = 'Trp', Y = 'Tyr', V = 'Val')
  
  # count occurrences of ANY AA
  pattern <- paste(mappingVect, collapse = "|")
  n_matches <- str_count(X, pattern)
  
  if( ! is.na(X) ) { 
    if (n_matches %in% 2) {
      return(str_replace_all(X, setNames(names(mappingVect), mappingVect)))
    } else {
      return(X)
    }
  } else {
    return(X)
  }
}

#load in unfiltered variants
booklet = readxl::read_excel("~/Documents/CHIP/sagi_automated_varcalls.xlsx", sheet = 1)

#@@@Step 1: filtering out candidate mutations based on QC features@@@#

#set all unlisted gnomAD variant frequencies to be 0
booklet$gnomAD_Genome_AF[which(is.na(booklet$gnomAD_Genome_AF))] <- 0

#filter out common gnomAD variants, inconsequential variants
first_pass = booklet[which(booklet$gnomAD_Genome_AF < 0.001 & booklet$Consequence %notin% c(
  "intron_variant", "synonymous_variant", "3_prime_UTR_variant")),]

#filter out variants without sufficient read support in one or more replicates
SSCS_support = unlist(lapply(first_pass$SSCS_Cov_Matrix, min_cov_matrix))
second_pass = first_pass[which(SSCS_support == "PASS"),]

#filter based on SSCS weighted VAF > threshold (0.01)
#output is on the order of 1e2 down from 1e4 
third_pass = second_pass[which(second_pass$SSCS_Weighted_VAF > 0.01),]

#@@@Step 2: identifying candidate driver mutations based on Vlasschaert rules@@@#

#add in column to third_pass df to store Vlasschaert status, add in a column to store rationale for functional prediction
third_pass$Vlasschaert = rep(0, nrow(third_pass))
third_pass$rationale = rep(0, nrow(third_pass))

#check for matching missense variants
missense_rules = read.table("~/Documents/CHIP/CHIP_missense_vars_mmb_03072024.txt", sep = "\t", header = T)
HGVSp = gsub("p.", replacement = "", third_pass$HGVSp, fixed = T)
HGVSp_subbed = unlist(lapply(HGVSp, aa_conversion))
#regex matching to find single amino acid substitutions
#technically can omit this step and update to go directly to gene + AA matching, but this is probably faster
missense_ix = which(str_detect(HGVSp_subbed, "^[A-Z][0-9]+[A-Z]$"))
missense_genes = third_pass$Gene_Symbol[missense_ix]

#find matching AA substitutions between ruleset and proposed missense variants
match_ix <- paste(missense_genes, HGVSp_subbed) %in%
paste(missense_rules$Gene, missense_rules$AAChange)

#update dataframe with AA substitution rows
third_pass$Vlasschaert[match_ix] = 1
third_pass$rationale[match_ix] = "missense_rule"

#check for matching nonsense variants
nonsense_rules = read.table("~/Documents/CHIP/CHIP_nonsense_FS_vars_mmb_03072024.txt", sep = "\t", header = T)
#checks for matching gene names where there is a stop gained in the HGVSp column
match_ix = which(third_pass$Gene_Symbol %in% nonsense_rules$Gene & (grepl("*", 
                                                                         third_pass$HGVSp, fixed = T) | 
                                                                      grepl("Ter",
                                                                            third_pass$HGVSp, ignore.case = T)))

#update dataframe with nonsense rows
third_pass$Vlasschaert[match_ix] = 1
third_pass$rationale[match_ix] = paste0(third_pass$rationale[match_ix], "_nonsense_rule")

#check for matching frameshift variants - same list as nonsense rules (only 'NPM1' needs additional logic later)
#checks for matching gene names where there is a frameshift in the consequence column
match_ix = which(third_pass$Gene_Symbol %in% nonsense_rules$Gene & grepl("frameshift", 
                                                                         third_pass$Consequence, ignore.case = T))

#update dataframe with frameshift rows
third_pass$Vlasschaert[match_ix] = 1
third_pass$rationale[match_ix] = paste0(third_pass$rationale[match_ix], "_frameshift_rule")
                                                              
#find matching splice mutations between ruleset and proposed variants
splice_rules = read.table("~/Documents/CHIP/CHIP_splice_vars_mmb_03072024.txt", sep = "\t", header = T)
#remove ASXL1/2 from this list - exon-specific criteria, deal with this separately 
splice_rules = splice_rules[which(splice_rules$Gene %notin% c("ASXL1", "ASXL2")),]
match_ix = which(third_pass$Gene_Symbol %in% splice_rules & grepl("splic", 
                                                                         third_pass$Consequence, ignore.case = T))

#update dataframe with splicing rows
third_pass$Vlasschaert[match_ix] = 1
third_pass$rationale[match_ix] = paste0(third_pass$rationale[match_ix], "_splicing_rule")

#@@@Step 3: additional candidate driver mutations based on gene-specific Vlasschaert rules@@@#
#implement ASXL1/2 and PPM1D logic
#looking for mutations in base(s) matching exons of interest for ASXL1
#add a one base buffer to either end to ensure we are catching splice sites in our exon definitions
#should not be an issue w.r.t. including intronic mutations due to two layers of filters 
#we are ignoring intronic mutations (1) and ASXL1 missense mutations (2)
ASXL1_exon12 = which(third_pass$Gene_Symbol == "ASXL1" & 
                       as.numeric(third_pass$Pos) >= (as.numeric(ASXL1_exons["exon_12", "exon_chrom_start"]) - 1) &
                       as.numeric(third_pass$Pos) <= (as.numeric(ASXL1_exons["exon_12", "exon_chrom_end"]) + 1)) 

ASXL1_exon13 = which(third_pass$Gene_Symbol == "ASXL1" & 
                       as.numeric(third_pass$Pos) >= (as.numeric(ASXL1_exons["exon_13", "exon_chrom_start"]) - 1) &
                       as.numeric(third_pass$Pos) <= (as.numeric(ASXL1_exons["exon_13", "exon_chrom_end"]) + 1))
ASXL1_ix = c(ASXL1_exon12, ASXL1_exon13)

#check for splice site mutations, add to dataframe
match_ix = which(grepl("splic", third_pass$Consequence[ASXL1_ix], ignore.case = T))
third_pass$Vlasschaert[ASXL1_ix[match_ix]] = 1
third_pass$rationale[ASXL1_ix[match_ix]] = paste0(third_pass$rationale[match_ix], "_ASXL1_splicing_rule")

#check for nonsense mutations, add to dataframe
match_ix = which(grepl("*", third_pass$HGVSp[ASXL1_ix], fixed = T) | 
                   grepl("Ter",
                         third_pass$HGVSp[ASXL1_ix], ignore.case = T))

third_pass$Vlasschaert[ASXL1_ix[match_ix]] = 1
third_pass$rationale[ASXL1_ix[match_ix]] = paste0(third_pass$rationale[match_ix], "_ASXL1_nonsense_rule")

#check for frameshift mutations, add to dataframe
match_ix = which(grepl("frameshift", third_pass$Consequence[ASXL1_ix], ignore.case = T))
third_pass$Vlasschaert[ASXL1_ix[match_ix]] = 1
third_pass$rationale[ASXL1_ix[match_ix]] = paste0(third_pass$rationale[match_ix], "_ASXL1_frameshift_rule")

#looking for mutations in base(s) matching exons of interest for ASXL2
ASXL2_exon12 = which(third_pass$Gene_Symbol == "ASXL2" & 
                       as.numeric(third_pass$Pos) >= (as.numeric(ASXL2_exons["exon_12", "exon_chrom_start"]) - 1) &
                       as.numeric(third_pass$Pos) <= (as.numeric(ASXL2_exons["exon_12", "exon_chrom_end"]) + 1)) 

ASXL2_exon13 = which(third_pass$Gene_Symbol == "ASXL2" & 
                       as.numeric(third_pass$Pos) >= (as.numeric(ASXL2_exons["exon_13", "exon_chrom_start"]) - 1) &
                       as.numeric(third_pass$Pos) <= (as.numeric(ASXL2_exons["exon_13", "exon_chrom_end"]) + 1))
ASXL2_ix = c(ASXL2_exon12, ASXL2_exon13)

#check for splice site mutations, add to dataframe
match_ix = which(grepl("splic", third_pass$Consequence[ASXL2_ix], ignore.case = T))
third_pass$Vlasschaert[ASXL2_ix[match_ix]] = 1
third_pass$rationale[ASXL2_ix[match_ix]] = paste0(third_pass$rationale[match_ix], "_ASXL2_splicing_rule")

#check for nonsense mutations, add to dataframe
match_ix = which(grepl("*", third_pass$HGVSp[ASXL2_ix], fixed = T) | 
                   grepl("Ter",
                         third_pass$HGVSp[ASXL2_ix], ignore.case = T))
third_pass$Vlasschaert[ASXL2_ix[match_ix]] = 1
third_pass$rationale[ASXL2_ix[match_ix]] = paste0(third_pass$rationale[match_ix], "_ASXL2_nonsense_rule")

#check for frameshift mutations, add to dataframe
match_ix = which(grepl("frameshift", third_pass$Consequence[ASXL2_ix], ignore.case = T))
third_pass$Vlasschaert[ASXL2_ix[match_ix]] = 1
third_pass$rationale[ASXL2_ix[match_ix]] = paste0(third_pass$rationale[match_ix], "_ASXL2_frameshift_rule")

#looking for mutations in base(s) matching exons of interest for PPM1D
PPM1D_exon5 = which(third_pass$Gene_Symbol == "PPM1D" & 
                       as.numeric(third_pass$Pos) >= (as.numeric(PPM1D_exons["exon_5", "exon_chrom_start"]) - 1) &
                       as.numeric(third_pass$Pos) <= (as.numeric(PPM1D_exons["exon_5", "exon_chrom_end"]) + 1)) 

PPM1D_exon6 = which(third_pass$Gene_Symbol == "PPM1D" & 
                       as.numeric(third_pass$Pos) >= (as.numeric(PPM1D_exons["exon_6", "exon_chrom_start"]) - 1) &
                       as.numeric(third_pass$Pos) <= (as.numeric(PPM1D_exons["exon_6", "exon_chrom_end"]) + 1))
PPM1D_ix = c(PPM1D_exon5, PPM1D_exon6)

#check for nonsense mutations, add to dataframe
match_ix = which(grepl("*", third_pass$HGVSp[PPM1D_ix], fixed = T) | 
                   grepl("Ter",
                         third_pass$HGVSp[PPM1D_ix], ignore.case = T))
third_pass$Vlasschaert[PPM1D_ix[match_ix]] = 1
third_pass$rationale[PPM1D_ix[match_ix]] = paste0(third_pass$rationale[match_ix], "_PPM1D_nonsense_rule")

#check for frameshift mutations, add to dataframe
match_ix = which(grepl("frameshift", third_pass$Consequence[PPM1D_ix], ignore.case = T))
third_pass$Vlasschaert[PPM1D_ix[match_ix]] = 1
third_pass$rationale[PPM1D_ix[match_ix]] = paste0(third_pass$rationale[match_ix], "_PPM1D_frameshift_rule")

#implement TET2 logic
residues = as.integer(gsub("[^[:digit:]]", "", HGVSp_subbed))
match_ix = (which(third_pass$Gene_Symbol == "TET2" & grepl("^[A-Z][0-9]+[A-Z]", HGVSp_subbed) &
                    ((residues>=1104&residues<=1481)|(residues>=1843&residues<=2002))))
third_pass$Vlasschaert[match_ix] = 1
third_pass$rationale[match_ix] = paste0(third_pass$rationale[match_ix], "_TET2_missense_rule")

#implement CBL logic
residues = as.integer(gsub("[^[:digit:]]", "", HGVSp_subbed))
match_ix = (which(third_pass$Gene_Symbol == "CBL" &  grepl("^[A-Z][0-9]+[A-Z]", HGVSp_subbed) & 
                    residues>=381 & residues<=421))
third_pass$Vlasschaert[match_ix] = 1
third_pass$rationale[match_ix] = paste0(third_pass$rationale[match_ix], "_CBL_missense_rule")

#implement CBLB logic
residues = as.integer(gsub("[^[:digit:]]", "", HGVSp_subbed))
match_ix = (which(third_pass$Gene_Symbol == "CBLB"  & grepl("^[A-Z][0-9]+[A-Z]", HGVSp_subbed) 
                  & residues>=372 & residues<=412))
third_pass$Vlasschaert[match_ix] = 1
third_pass$rationale[match_ix] = paste0(third_pass$rationale[match_ix], "_CBLB_missense_rule")

#implement ZBTB33 logic
residues = as.integer(gsub("[^[:digit:]]", "", HGVSp_subbed))
match_ix = (which(third_pass$Gene_Symbol == "ZBTB33" & grepl("^[A-Z][0-9]+[A-Z]", HGVSp_subbed) &
                    ((residues>=9&residues<=126)|(residues>=332&residues<=591))))
third_pass$Vlasschaert[match_ix] = 1
third_pass$rationale[match_ix] = paste0(third_pass$rationale[match_ix], "_ZBTB_missense_rule")

#@@@Step 4: flag certain mutations for mannual review@@@#

#implement in-frame indel logic
match_ix = (which(third_pass$Gene_Symbol %in% c("ZBTB33", "CREBBP", "DNMT3A", "EP300", 
                                                "FLT3", "JAK2", "KDM6A", "KIT", "MPL") & (
                        grepl("inframe_del", third_pass$Consequence, ignore.case = T) |
                        grepl("inframe_ins", third_pass$Consequence, ignore.case = T)) 
            ))
third_pass$Vlasschaert[match_ix] = 2
third_pass$rationale[match_ix] = paste0(third_pass$rationale[match_ix], "_check_indel_rule")


#implement CSF3R	truncating c.741-791 logic
match_ix = (which(third_pass$Gene_Symbol == "CSF3R" & (
  grepl("del", third_pass$Consequence, ignore.case = T) |
    grepl("*", third_pass$HGVSp, fixed = T) |
    grepl("Ter", third_pass$HGVSp, ignore.case = T)) 
))
third_pass$Vlasschaert[match_ix] = 2
third_pass$rationale[match_ix] = paste0(third_pass$rationale[match_ix], "_CSF3R_truncating_rule")

#implement NPM1	frameshift p.W288fs logic
match_ix = (which(third_pass$Gene_Symbol == "NPM1" &
  grepl("frameshift", third_pass$Consequence, ignore.case = T)))
third_pass$Vlasschaert[match_ix] = 2
third_pass$rationale[match_ix] = paste0(third_pass$rationale[match_ix], "_CSF3R_truncating_rule")

#clean up trailing '0_' from third_pass$rationale, order by ML score
third_pass$rationale = gsub("0_", "", third_pass$rationale)
third_pass = third_pass[order(third_pass$ML_Score, decreasing = T),]

#@@@Step 5: add original sample IDs back by mapping to OICR IDs@@@#
sample_map = fread("~/Documents/CHIP/CHTC.release.02.24.2026.fastqs.map.tsv")
OICR_IDs = unlist(lapply(third_pass$Sample_ID, ss_sample_string))
#find rows of sample_map which correspond to OCIR IDs enabling samples with multiple mutations to multimap
ix_match = match(OICR_IDs, sample_map$library_id)
third_pass$Donor_ID = sample_map$donor_id[ix_match]
third_pass = third_pass %>% relocate(Donor_ID, .after = Sample_ID)


#write sheet to booklet 
wb = loadWorkbook("~/Documents/CHIP/sagi_automated_varcalls.xlsx")
addWorksheet(wb, "automatedCalling")
writeData(wb, sheet = "automatedCalling", x = third_pass, colNames = T, rowNames = F)

#remove unused sheets
#removeWorksheet(wb, "aggregated_annotated_variants_2")
#removeWorksheet(wb, "Consequence+gnomAD")

#format automated sheet based on overlap with ML method
greenStyle <- createStyle(fontColour = "#000000", fgFill = "green")
yellowStyle <- createStyle(fontColour = "#000000", fgFill = "yellow")
redStyle <- createStyle(fontColour = "#000000", fgFill = "red")
orangeStyle <- createStyle(fontColour = "#000000", fgFill = "orange")

# Apply colourStyle:
greenRows = which(third_pass$Predicted_Label == 1 & third_pass$Vlasschaert == 1)
redRows = which(third_pass$Predicted_Label == 1 & third_pass$Vlasschaert != 1)
yellowRows = which(third_pass$Predicted_Label != 1 & third_pass$Vlasschaert == 1)
orangeRows = which(third_pass$Vlasschaert == 2)

addStyle(wb, "automatedCalling", cols = ncol(third_pass) - 1, rows = greenRows +1,
         style = greenStyle, gridExpand = TRUE)
addStyle(wb, "automatedCalling", cols = ncol(third_pass) - 1, rows = yellowRows +1,
         style = yellowStyle, gridExpand = TRUE)
addStyle(wb, "automatedCalling", cols = ncol(third_pass) - 1, rows = redRows +1,
         style = redStyle, gridExpand = TRUE)
addStyle(wb, "automatedCalling", cols = ncol(third_pass) - 1, rows = orangeRows +1,
         style = orangeStyle, gridExpand = TRUE)

saveWorkbook(wb, "~/Documents/CHIP/daniel_automated_annotations.xlsx", overwrite = T)
