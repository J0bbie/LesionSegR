library(dplyr)
library(future)
library(LesionSegR)
library(VariantAnnotation)

# Parallel settings.
future::plan(future::multisession, workers = 10)

# Import metadata of Ultima/Illumuna samples. ----

metadata <- readxl::read_xlsx(
    path = "~/odomLab/LesionSegregration_F1/manuscript/tables/SupplTable1_OverviewSequencing.xlsx",
    sheet = "WGS (NovaSeq)", trim_ws = TRUE
) %>%
    # Only keep the malignant samples.
    dplyr::filter(!is.na(matched_group)) %>% 
    dplyr::filter(sample %in% c('AS-949292', 'AS-949296', 'AS-949300'))

# Import somatic variants and CNV. ----

## Import somatic variants (Illumina). ----
files_vcf <- base::list.files(path = "~/odomLab/LesionSegregration_F1/data/workflow/variants/", pattern = base::paste(base::paste0(metadata$sample, ".*_haplocounted.vcf.gz$"), collapse = "|"), full.names = TRUE)
somatics_illumina <- LesionSegR::import_vcf(files_vcf)

## Import somatic variants (Ultima). ----
files_vcf <- base::list.files(path = "~/odomLab/LesionSegregration_F1/UltimaGenomics/variantCalling/Mutect2_haplotyped/", pattern = ".*_haplocounted.vcf.gz$", full.names = TRUE)
somatics_ultima <- LesionSegR::import_vcf_mutect(files_vcf)

## Convert sample_names to similar names.
IRanges::subsetByOverlaps(somatics_illumina$`AS-949292_B6_CAST_EiJ`, somatics_ultima$CABL1t1)
somatics_illumina$`AS-949296_B6_CAST_EiJ` somatics_ultima$CABL2t1
