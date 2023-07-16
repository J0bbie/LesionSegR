
# Check status of heterogeneity of CAST-SNPs in all mice. ----

f <- "/omics/groups/OE0538/internal/projects/sharedData/GRCm39/genome/SNPSplit/B6_CAST-EiJ/all_SNPs_CAST_EiJ_GRCm39.txt.gz"
snps <- readr::read_tsv(f, col_names = c("ID", "chr", "start", "width", "allele"), show_col_types = FALSE) %>%
    dplyr::mutate(ref = stringr::str_sub(allele, 1, 1), alt = stringr::str_sub(allele, 3, 3)) %>%
    dplyr::mutate(haplotype = paste(pmin(ref, alt), pmax(ref, alt), sep = "/")) %>%
    dplyr::select(chr, start = start, end = start, ref, alt, haplotype) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE) %>%
    GenomicRanges::sort()

# Determine possibility of SNP
z = IRanges::subsetByOverlaps(unlist(reciprocal_data), snps)


# Retrieve the no. of nucleotide at each SNP position.

get_pileup <- function(f, snps) {
    Rsamtools::pileup(
        Rsamtools::BamFile(f),
        scanBamParam = Rsamtools::ScanBamParam(which = snps),
        pileupParam = Rsamtools::PileupParam(distinguish_strands = FALSE, include_deletions = FALSE, max_depth = 100, include_insertions = FALSE)
    ) %>%
        dplyr::group_by(seqnames, pos) %>% 
        dplyr::filter(count >= 2) %>% 
        dplyr::summarise(
            nNucleotides = dplyr::n_distinct(nucleotide),
            haplotype = paste(pmin(nucleotide), collapse = "/")
            ) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(sample = gsub("_.*", "", basename(f)))
}

pileups <- list()
pileups$normal1 <- get_pileup("/omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/WGS/alignment/B6_CAST-EiJ/AS-949304_sortedByCoord_markDup_BQSR_alleleFlagged.bam", snps)
pileups$normal2 <- get_pileup("/omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/WGS/alignment/B6_CAST-EiJ/AS-949306_sortedByCoord_markDup_BQSR_alleleFlagged.bam", snps)

# saveRDS(pileups, "~/odomLab/LesionSegregration_F1/pileups_snps.rds")


pileups <- base::readRDS("~/odomLab/LesionSegregration_F1/pileups_snps.rds")
z <- dplyr::full_join(pileups$normal1, pileups$normal2, by = c("seqnames", "pos"))
zz = GenomicRanges::as.data.frame(snps) %>% dplyr::select(seqnames, pos = start, B6 = ref, CAST = alt, haplotype_MGP = haplotype)

zzz <- z %>% dplyr::right_join(zz, by = c("seqnames", "pos"))

zzz <- zzz %>% 
    dplyr::mutate(
        matchAny = dplyr::if_else(haplotype.x == haplotype_MGP | haplotype.y == haplotype_MGP, T, F),
        matchBoth = dplyr::if_else(haplotype.x == haplotype_MGP & haplotype.y == haplotype_MGP, T, F),
        matchB6 = dplyr::if_else(haplotype.x == B6 & haplotype.y == B6, T, F),
        matchCAST = dplyr::if_else(haplotype.x == CAST & haplotype.y == CAST, T, F)
    )

table(zzz$matchAny)
table(zzz$matchBoth)


x <- zzz[!zzz$matchBoth,]
table(x$haplotype.x == x$B6)
table(x$haplotype.x == x$CAST)