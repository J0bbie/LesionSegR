# Function: Generates a figure detailing the haplotype assignment of somatic variants (G1 vs. G2) per sample.
# Author: J. van Riet

# Load libraries. ----

library(dplyr)
library(patchwork)
library(extrafont)

# Import functions. ----

source("R/functions.R")
source("R/themes.R")


# Import data. ----

reciprocal_results <- base::readRDS("data/reciprocal_results.rds")
reciprocal_data <- base::readRDS("data/reciprocal_somaticvariants.rds")


# Import metadata. ----

metadata <- readxl::read_xlsx(
    path = "~/Dropbox/Work/DKFZ/Projects/Odom - LesionSegregation/Draft/Tables/SupplTable1_OverviewSequencing.xlsx",
    sheet = "WGS", trim_ws = TRUE
) %>%
    dplyr::filter(ASID %in% reciprocal_results$mutationalBurden$sample)


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

## QC - SNPsplit ----

reciprocal_results$snpsplit %>%
    dplyr::filter(!grepl("f1-liver-normal", description)) %>%
    dplyr::left_join(reciprocal_results$mutationalBurden) %>%
    dplyr::left_join(reciprocal_results$flagstats) %>%
    dplyr::mutate(
        PercReadsWithN = formattable::percent((.$N_containing_reads / .$total_reads), digits = 2),
        percent_g1 = formattable::percent(percent_g1 / 100, digits = 2),
        percent_g2 = formattable::percent(percent_g2 / 100, digits = 2),
        percent_unassignable = formattable::percent(percent_unassignable / 100, digits = 2),
        TMB = formattable::color_tile("grey90", max.color = "red")(round(TMB, 2)),
        total_reads = formattable::color_tile("lightpink", max.color = "hotpink")(base::formatC(.$total_reads, drop0trailing = TRUE, big.mark = ",")),
        properPairPerc = formattable::percent(`properly paired` / `paired in sequencing`, digits = 2)
    ) %>%
    dplyr::select(
        ASID = sample,
        Descr. = description,
        "Total reads" = total_reads,
        "Mapped" = mapped,
        "Paired" = `paired in sequencing`,
        "Properly Paired (%)" = properPairPerc,
        "N-reads" = PercReadsWithN,
        "B6" = percent_g1,
        "CAST-EiJ" = percent_g2,
        "UA" = percent_unassignable,
        TMB,
        "Total" = totalMutations,
        "B6 " = totalG1,
        "CAST-EiJ " = totalG2,
        "UA " = totalUA,
    ) %>%
    dplyr::arrange(-Total) %>% 
    knitr::kable(escape = FALSE, align = "llcccccccccc", format.args = list(decimal.mark = '.', big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15, full_width = TRUE, html_font = "Lexend", bootstrap_options = c("striped")) %>%
    kableExtra::add_header_above(line_sep = 10, c(" " = 3, "Flagstats" = 3, "SNPSplit (% reads)" = 4, "Somatic variants (no.)" = 5))


# QC - Visualize the distribution of strain-assigned somatic variants. ----

dplyr::bind_rows(
    base::lapply(reciprocal_data, function(x) tibble::as_tibble(S4Vectors::mcols(x)[c("sample", "G1_Alt", "G2_Alt", "origin_mutant")]))
) %>%
    dplyr::inner_join(metadata, by = c("sample" = "ASID")) %>% 
    dplyr::mutate(
        delta = G1_Alt - G2_Alt,
        assignment = dplyr::case_when(
            origin_mutant == "G1" ~ "B6",
            origin_mutant == "G2" ~ "CAST",
            origin_mutant == "UA" ~ "Unassigned"
        )
    ) %>%
    # Remove 0 deltas
    dplyr::filter(delta != 0) %>%
    # Plot distribution.
    ggplot2::ggplot(., ggplot2::aes(x = delta, fill = assignment)) +
    ggplot2::geom_histogram(binwidth = .5, na.rm = TRUE, color = "grey10", lwd = ggplot2::rel(.33)) +
    ggplot2::scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
    ggplot2::scale_y_continuous(expand = c(0, 0.3), trans = scales::sqrt_trans(), breaks = c(0, 10, 50, 100, 250, 500, 1000, 2500, 5000, 7500), limits = c(0, 7500)) +
    ggplot2::scale_fill_manual(values = color_scheme, name = NULL) +
    ggplot2::labs(
        x = "B6 - CAST SNP-assigned reads",
        y = "Frequency<br><sup>(No. of somatic variants)</sup>"
    ) +
    ggplot2::facet_wrap(~ description, ncol = 3) +
    theme_job
