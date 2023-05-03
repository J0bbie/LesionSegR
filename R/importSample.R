# Function: Generate overview of samples.

# Libraries ----

library(dplyr)

# Import list of sequencing samples ----

samples_sequencing <- list()
samples_sequencing$WGS <- readxl::read_xlsx("/omics/groups/OE0538/internal/projects/LesionSegregration_F1/metadata/SupplTable1_OverviewSequencing.xlsx", sheet = "WGS", trim_ws = T)
samples_sequencing$WTS <- readxl::read_xlsx("/omics/groups/OE0538/internal/projects/LesionSegregration_F1/metadata/SupplTable1_OverviewSequencing.xlsx", sheet = "WTS", trim_ws = T)

# Determine matching WTS.
samples_sequencing$WGS$hasMatchingWTS <- gsub(".*-", "", samples_sequencing$WGS$description) %in% gsub(".*-", "", samples_sequencing$WTS$description)
samples_sequencing$WTS$hasMatchingWGS <- gsub(".*-", "", samples_sequencing$WTS$description) %in% gsub(".*-", "", samples_sequencing$WGS$description)


# Read in VCF. ----

pVCF <- "/omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/WGS/variantCalling/B6_CAST-EiJ/AS-949280_withStrainCounts.vcf.gz"

z <- VariantAnnotation::readVcf(pVCF, genome = "GRCm39")
z <- z[z@fixed$FILTER == "PASS", ]

# Expand multi-allelic hits to separate rows.
z <- VariantAnnotation::expand(z)

# Generate overview of UA / G1 / G2 per haplotype.
x <- tibble::as_tibble(data.frame(
    seqname = GenomicRanges::seqnames(z),
    pos = GenomicRanges::start(z),
    ref = VariantAnnotation::ref(z),
    alt = VariantAnnotation::alt(z),
    UA = VariantAnnotation::geno(z)$UA[, , 2],
    G1 = VariantAnnotation::geno(z)$G1[, , 2],
    G2 = VariantAnnotation::geno(z)$G2[, , 2]
))

x <- x %>% dplyr::mutate(strainAlt = ifelse(G1 >= 3 | G2 >= 3, ifelse(G1 / G2 >= 3, "G1", ifelse(G2 / G1 >= 3, "G2", "AU")), "AU"))
