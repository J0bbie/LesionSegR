library(dplyr)

f <- "/omics/groups/OE0538/internal/projects/sharedData/GRCm39/genome/SNPSplit/B6_CAST-EiJ/all_SNPs_CAST_EiJ_GRCm39.txt.gz"

snps <- readr::read_tsv(f, col_names = c("ID", "chr", "start", "width", "allele"), show_col_types = FALSE) %>%
    dplyr::mutate(ref = stringr::str_sub(allele, 1, 1), alt = stringr::str_sub(allele, 3, 3)) %>%
    dplyr::select(chr, start = start, end = start, ref, alt)

snps <- GenomicRanges::makeGRangesFromDataFrame(snps, keep.extra.columns = TRUE) %>%
    GenomicRanges::sort() %>% 
    VariantAnnotation::makeVRangesFromGRanges()

VariantAnnotation::sampleNames(snps) <- "MGP"

VariantAnnotation::writeVcf(snps, "/omics/groups/OE0538/internal/projects/sharedData/GRCm39/genome/SNPSplit/B6_CAST-EiJ/all_SNPs_CAST_EiJ_GRCm39.vcf.gz", index = TRUE)    
