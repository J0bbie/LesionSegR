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
files_vcf <- list.files(path = "~/Downloads/VCF", pattern = "*withStrainCounts.vcf.gz$", full.names = TRUE)
data_somaticvariants <- pbapply::pblapply(files_vcf, import_vcf, cl = 4)
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

## Determine no. of Ti/Tv and ratio. -----

reciprocal_results$TiTv <- NULL

## dN/dS. ----

reciprocal_results$dNdS <- NULL


# Save the data. ----

saveRDS(data_somaticvariants, "data/reciprocal_somaticvariants.rds")
saveRDS(reciprocal_results, "data/reciprocal_results.rds")
