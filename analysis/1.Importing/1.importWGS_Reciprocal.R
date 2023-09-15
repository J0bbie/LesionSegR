# Function: Import the WGS data of the F1 mice.
# Author: J. van Riet

# Libraries ----

library(dplyr)
library(future)
library(LesionSegR)

# Parallel settings.
future::plan(future::multisession, workers = 10)

# Import metadata. ----

metadata <- readxl::read_xlsx(
    path = "~/odomLab/LesionSegregration_F1/manuscript/tables/SupplTable1_OverviewSequencing.xlsx",
    sheet = "WGS (NovaSeq)", trim_ws = TRUE
) %>%
    # Only keep the malignant samples.
    dplyr::filter(!is.na(matched_group))

# Import somatic variants and CNV. ----

## Import somatic variants. ----
files_vcf <- base::list.files(path = "~/odomLab/LesionSegregration_F1/testRun/variants/", pattern = base::paste(base::paste0(metadata$sample, ".*_haplocounted.vcf.gz$"), collapse = "|"), full.names = TRUE)
data_somaticvariants <- LesionSegR::import_vcf(files_vcf)


## Determine WGS characteristics. ----

reciprocal_results <- list()

## Flagstats. ----

files_flagstats <- list.files(path = "~/odomLab/LesionSegregration_F1/testRun/qc/", pattern = "*.flagstats$", full.names = TRUE)
reciprocal_results$flagstats <- LesionSegR::read_flagstats(files_flagstats)

## Determine mutational burden. ----
reciprocal_results$mutationalBurden <- LesionSegR::determine_mutational_burden(data_somaticvariants)


## Generate the 96-context matrices. ----

reciprocal_results$mut_matrixes_96 <- LesionSegR::generate_mutmatrices_96(data_somaticvariants)

## Perform bootstrapped COSMIC-signature analysis. ----

signatures_sbs <- MutationalPatterns::get_known_signatures(muttype = "snv", source = "COSMIC_v3.2", genome = "mm10")
signatures_dbs <- MutationalPatterns::get_known_signatures(muttype = "dbs", source = "COSMIC_v3.2", genome = "mm10")
signatures_indel <- MutationalPatterns::get_known_signatures(muttype = "indel", source = "COSMIC_v3.2", genome = "GRCh37")

reciprocal_results$signature_fit_sbs <- MutationalPatterns::fit_to_signatures_bootstrapped(reciprocal_results$mut_matrixes_96$sbs, signatures_sbs, n_boots = 100, method = "strict")
reciprocal_results$signature_fit_dbs <- MutationalPatterns::fit_to_signatures_bootstrapped(reciprocal_results$mut_matrixes_96$dbs, signatures_dbs, n_boots = 100, method = "strict")
reciprocal_results$signature_fit_indel <- MutationalPatterns::fit_to_signatures_bootstrapped(reciprocal_results$mut_matrixes_96$indel, signatures_indel, n_boots = 100, method = "strict")

## Determine no. of Ti/Tv and ratio. -----

reciprocal_results$ti_tv <- LesionSegR::determine_ti_tv(reciprocal_results$mut_matrixes_96$sbs)

## dN/dS. ----

reciprocal_results$dNdS <- LesionSegR::run_dnds(data_somaticvariants, path_db = "~/Downloads/refCDS_ENSEMBLv109_GRCm39.rda")


# Save the data. ----

saveRDS(data_somaticvariants, "data/reciprocal_somaticvariants.rds")
saveRDS(reciprocal_results, "data/reciprocal_results.rds")
