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
    sheet = "WGS", trim_ws = TRUE
)

# Read in somatic variants. ----

## Import somatic variants. ----
files_vcf <- list.files(path = "~/Downloads/VCF/", pattern = "*withStrainCounts.vcf.gz$", full.names = TRUE)
data_reciprocal <- pbapply::pblapply(files_vcf, import_vcf, cl = 4)

## Determine mutational burden. ----
mutational_info <- dplyr::bind_rows(lapply(data_reciprocal, calculate_tmb)) %>%
    dplyr::inner_join(metadata, by = c("sample" = "ASID"))

## Mutational signature analysis. ----
data_mutsigs <- lapply(data_reciprocal, determine_knownsignatures)
MutationalPatterns::plot_96_profile(data_mutsigs[[12]]$sbs$mutMatrix, condensed = TRUE, ymax = 0.2)
MutationalPatterns::plot_contribution(data_mutsigs[[12]]$sbs$fit_res$contribution, data_mutsigs$sbs$signature, mode = "relative")


## QC - SNPsplit ----

# Read yaml files.
read_snpsplit <- function(x) {
    tibble::as_tibble(yaml::read_yaml(x)$Tagging) %>%
        dplyr::mutate(
            sample = basename(x) %>% stringr::str_replace_all(pattern = "_.*", replacement = "")
        )
}

files_snpsplit <- list.files(path = "~/DKFZ/odomLab/LesionSegregration_F1/data/WGS/alignment/B6_CAST-EiJ/", pattern = "*.SNPsplit_report.yaml$", full.names = TRUE)
data_snpsplit <- dplyr::bind_rows(lapply(files_snpsplit, read_snpsplit))
data_snpsplit <- data_snpsplit %>% dplyr::inner_join(metadata, by = c("sample" = "ASID"))

data_snpsplit %>%
    dplyr::filter(!grepl("f1-liver-normal", description)) %>%
    dplyr::left_join(mutational_info) %>%
    dplyr::mutate(
        PercReadsWithN = formattable::percent((.$N_containing_reads / .$total_reads), digits = 2),
        percent_g1 = formattable::percent(percent_g1 / 100, digits = 2),
        percent_g2 = formattable::percent(percent_g2 / 100, digits = 2),
        percent_unassignable = formattable::percent(percent_unassignable / 100, digits = 2),
        TMB = formattable::color_tile("grey90", max.color = "red")(round(TMB, 2)),
        total_reads = formattable::color_tile("lightpink", max.color = "hotpink")(total_reads)

    ) %>%
    dplyr::select(
        ASID = sample,
        Descr. = description,
        "Total reads" = total_reads,
        "N-reads" = PercReadsWithN,
        "B6" = percent_g1,
        "CAST-EiJ" = percent_g2,
        "UA" = percent_unassignable,
        TMB,
        "Total" = totalMutations,
        "B6 " = totalG1,
        "CAST-EiJ " = totalG2,
        "AU " = totalAU,
    ) %>%
    knitr::kable(escape = FALSE, align = "llcccccccccc") %>%
    kableExtra::kable_styling(font_size = 15, full_width = TRUE, html_font = "Lexend", bootstrap_options = c("striped")) %>%
    kableExtra::add_header_above(line_sep = 10, c(" " = 3, "SNPSplit (%)" = 4, "Somatic variants (no.)" = 5))


# Detect recurrent mutations. ---

recurrent_mutations <- dplyr::bind_rows(lapply(data_reciprocal, function(x) {
    tibble::as_tibble(data.frame(chr = seqnames(x), x@ranges, DP = x@totalDepth, AD = x@altDepth, ref = x@ref, alt = x@alt, origin_mutant = x$origin_mutant, sample = x$sample))
})) %>%
    dplyr::left_join(metadata, by = c("sample" = "ASID")) %>%
    dplyr::group_by(chr, start, end, ref, alt) %>%
    dplyr::summarise(
        countSamples = dplyr::n_distinct(sample),
        totalCasB6 = sum(Type == "Tumor (Liver; Reciprocal CAS/B6)"),
        totalB6CAS = sum(Type == "Tumor (Liver; Reciprocal B6/CAS)"),
        totalG1 = sum(origin_mutant == "G1"),
        totalG2 = sum(origin_mutant == "G2"),
        totalUA = sum(origin_mutant == "UA"),
        minDP = max(DP),
        minAD = max(AD)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(countSamples > 1)


# Determine reciprocal-specific recurrencies.
IRanges::subsetByOverlaps(data_reciprocal[[8]], GenomicRanges::makeGRangesFromDataFrame(subset(recurrent_mutations, countSamples == 6)))



# CNV -------------------------------------------------------------------------

# Import the CNV profiles.
sprintf("\tImporting CNV profiles") %>% ParallelLogger::logTrace()

files_cnv <- list.files(path = "~/DKFZ/odomLab/LesionSegregration_F1/data/WGS/copynumbers/B6_CAST-EiJ/", pattern = "*modelFinal.seg", full.names = TRUE)

files_h5 <- list.files(path = "~/DKFZ/odomLab/LesionSegregration_F1/data/WGS/copynumbers/B6_CAST-EiJ/", pattern = "*hdf5", full.names = TRUE)

# Load read-counts per interval.
file <- files_h5[1]

x <- rhdf5::H5Fopen(file, flags = "H5F_ACC_RDONLY")
intervals <- data.frame(x$intervals$transposed_index_start_end)
intervals[, 1] <- x$intervals$indexed_contig_names[intervals[, 1] + 1]
targetCoverage <- GenomicRanges::GRanges(intervals[, 1], IRanges::IRanges(intervals[, 2], intervals[, 3]))
targetCoverage$counts <- x$counts$values[, 1]

z <- rhdf5::h5dump(path_cnv)

path_cnv <- files_cnv[1]
sample_cnv <- readr::read_delim(path_cnv, comment = "@") %>%
    dplyr::mutate(
        CONTIG = factor(CONTIG, levels = gtools::mixedsort(unique(CONTIG)))
    )

# Plot the CNV profile.
sample_cnv %>%
    ggplot2::ggplot(ggplot2::aes(x = START, xend = END, yend = LOG2_COPY_RATIO_POSTERIOR_10, y = LOG2_COPY_RATIO_POSTERIOR_10)) +
    ggplot2::geom_segment() +
    ggplot2::theme_bw() +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        legend.position = "none"
    ) +
    ggplot2::labs(
        x = "Position",
        y = "Copy number",
        title = "CNV profile"
    ) +
    ggplot2::facet_grid(~CONTIG, space = "free") +
    ggplot2::scale_x_continuous(labels = scales::comma)
