# Function: Import the WGS data of the F1 mice.
# Author: J. van Riet

# Libraries ----

library(dplyr)
library(ParallelLogger)

# Import functions. ----

source("R/functions.R")

# Import metadata. ----

metadata <- readxl::read_xlsx(
    path = "~/Dropbox/Work/DKFZ/Projects/Odom - LesionSegregation/Draft/Tables/SupplTable1_OverviewSequencing.xlsx",
    sheet = "WGS_Reciprocals", trim_ws = TRUE
) %>%
    # Only keep the malignant samples.
    dplyr::filter(!is.na(matched_group))

# Import somatic variants and CNV. ----

## Import somatic variants. ----
files_vcf <- list.files(path = "~/Downloads/VCF", pattern = paste(paste0(metadata$sample, "_withStrainCounts.vcf.gz$"), collapse = "|"), full.names = TRUE)
data_somaticvariants <- pbapply::pblapply(files_vcf, import_vcf, cl = 10)
names(data_somaticvariants) <- vapply(data_somaticvariants, function(x) unique(x$sample), "")

# Remove matched-normals.
data_somaticvariants <- data_somaticvariants[names(data_somaticvariants) %in% metadata$sample]

# Convert to VRangesList
data_somaticvariants <- VariantAnnotation::VRangesList(data_somaticvariants)


## Determine WGS characteristics. ----

reciprocal_results <- list()

## Flagstats. ----

files_flagstats <- list.files(path = "~/DKFZ/odomLab/LesionSegregration_F1/results/alignment/C57BL_6J_CAST_EiJ/WGS/", pattern = "*.flagstats$", full.names = TRUE)
reciprocal_results$flagstats <- tibble::as_tibble(dplyr::bind_rows(base::lapply(files_flagstats, read_flagstats)) %>% reshape2::dcast(., sample ~ variable)) %>%
    dplyr::filter(sample %in% metadata$sample)

## SNPSplit. ----

files_snpsplit <- list.files(path = "~/DKFZ/odomLab/LesionSegregration_F1/results/alignment/C57BL_6J_CAST_EiJ/WGS/", pattern = "*.SNPsplit_report.yaml$", full.names = TRUE)
reciprocal_results$snpsplit <- dplyr::bind_rows(base::lapply(files_snpsplit, read_snpsplit)) %>% dplyr::filter(sample %in% metadata$sample)


## Determine mutational burden. ----
reciprocal_results$mutationalBurden <- dplyr::bind_rows(base::lapply(data_somaticvariants, determine_mutational_burden))

## Generate the 96-context matrices. ----

reciprocal_results$mut_matrixes_96 <- generate_mutmatrices_96(data_somaticvariants)

## Perform bootstrapped COSMIC-signature analysis. ----

signatures_sbs <- MutationalPatterns::get_known_signatures(muttype = "snv", source = "COSMIC_v3.2", genome = "mm10")
signatures_dbs <- MutationalPatterns::get_known_signatures(muttype = "dbs", source = "COSMIC_v3.2", genome = "mm10")
signatures_indel <- MutationalPatterns::get_known_signatures(muttype = "indel", source = "COSMIC_v3.2", genome = "GRCh37")

reciprocal_results$signature_fit_sbs <- MutationalPatterns::fit_to_signatures_bootstrapped(reciprocal_results$mut_matrixes_96$sbs, signatures_sbs, n_boots = 50, method = "strict")
reciprocal_results$signature_fit_dbs <- MutationalPatterns::fit_to_signatures_bootstrapped(reciprocal_results$mut_matrixes_96$dbs, signatures_dbs, n_boots = 50, method = "strict")
reciprocal_results$signature_fit_indel <- MutationalPatterns::fit_to_signatures_bootstrapped(reciprocal_results$mut_matrixes_96$indel, signatures_indel, n_boots = 50, method = "strict")

## Determine no. of Ti/Tv and ratio. -----

reciprocal_results$ti_tv <- determine_ti_tv(reciprocal_results$mut_matrixes_96$sbs)

## dN/dS. ----

reciprocal_results$dNdS <- run_dnds(data_somaticvariants)


# Save the data. ----

saveRDS(data_somaticvariants, "data/reciprocal_somaticvariants.rds")
saveRDS(reciprocal_results, "data/reciprocal_results.rds")
