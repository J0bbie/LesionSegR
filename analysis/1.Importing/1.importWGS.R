# Function: Import the WGS data of the F1 mice.
# Author: J. van Riet

# Libraries ----

library(dplyr)
library(future)
library(LesionSegR)
library(VariantAnnotation)

# Parallel settings.
future::plan(future::multicore, workers = 10)

# Import metadata. ----

metadata <- readxl::read_xlsx(
    path = "~/odomLab/LesionSegregration_F1/manuscript/tables/SupplTable1_OverviewSequencing.xlsx",
    sheet = "WGS (NovaSeq)", trim_ws = TRUE
) %>%
    # Only keep the malignant samples.
    dplyr::filter(!is.na(matched_group)) %>% 
    dplyr::mutate(sample_strain = paste(sample, strain1, strain2, sep = '_'))

# Import somatic variants and CNV. ----

## Import somatic variants. ----
files_vcf <- base::list.files(path = "~/odomLab/LesionSegregration_F1/data/workflow/variants/", pattern = base::paste(base::paste0(metadata$sample, ".*_haplocounted.vcf.gz$"), collapse = "|"), full.names = TRUE)
data_somaticvariants <- LesionSegR::import_vcf(files_vcf)


## Determine WGS characteristics. ----

data_results <- list()

## Flagstats. ----

files_flagstats <- list.files(path = "~/odomLab/LesionSegregration_F1/data/workflow/alignment/WGS/", pattern = "*.flagstat$", full.names = TRUE)
data_results$flagstats <- LesionSegR::read_flagstats(files_flagstats)

## Haplotag ----

files_haplotag<- list.files(path = "~/odomLab/LesionSegregration_F1/data/workflow/logs/haplotyping/", pattern = "^haplotag.*log", full.names = TRUE)
data_results$haplotag <- LesionSegR::read_haplotag_log(files_haplotag) %>% 
    dplyr::filter(sample %in% metadata$sample_strain)


## Determine mutational burden. ----
data_results$mutationalBurden <- LesionSegR::determine_mutational_burden(data_somaticvariants)

## Generate the 96-context matrices. ----

data_results$mut_matrixes_96 <- LesionSegR::generate_mutmatrices_96(data_somaticvariants)

## Perform bootstrapped COSMIC-signature analysis. ----

signatures_sbs <- MutationalPatterns::get_known_signatures(muttype = "snv", source = "COSMIC_v3.2", genome = "mm10")
signatures_indel <- MutationalPatterns::get_known_signatures(muttype = "indel", source = "COSMIC_v3.2", genome = "GRCh37")

data_results$signature_fit_sbs <- MutationalPatterns::fit_to_signatures_bootstrapped(data_results$mut_matrixes_96$sbs, signatures_sbs, n_boots = 100, method = "strict")
data_results$signature_fit_indel <- MutationalPatterns::fit_to_signatures_bootstrapped(data_results$mut_matrixes_96$indel, signatures_indel, n_boots = 100, method = "strict")

## Determine no. of Ti/Tv and ratio. -----

data_results$ti_tv <- LesionSegR::determine_ti_tv(data_results$mut_matrixes_96$sbs)

## dN/dS. ----

data_results$dNdS <- LesionSegR::run_dnds(data_somaticvariants, path_db = "/omics/groups/OE0538/internal/projects/sharedData/GRCm39/annotation/refCDS_ENSEMBLv110_GRCm39.rda")
data_results$dNdS <- data_results$dNdS$finalOutput

# Save the data. ----

saveRDS(data_somaticvariants, "~/odomLab/LesionSegregration_F1/data/rdata/data_somaticvariants.rds")
saveRDS(data_results, "~/odomLab/LesionSegregration_F1/data/rdata/data_results.rds")