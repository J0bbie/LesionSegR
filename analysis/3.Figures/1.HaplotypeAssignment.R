# Function: Generates a figure detailing the haplotype assignment of somatic variants (G1 vs. G2) per sample.
# Author: J. van Riet

# Load libraries. ----

library(dplyr)
library(patchwork)
library(extrafont)

# Import functions. ----

source("analysis/themes.R")


# Import data. ----

data_somaticvariants <- base::readRDS("~/odomLab/LesionSegregration_F1/data/rdata/data_somaticvariants.rds")
data_results <- base::readRDS("~/odomLab/LesionSegregration_F1/data/rdata/data_results.rds")


# Import metadata. ----

metadata <- readxl::read_xlsx(
    path = "~/odomLab/LesionSegregration_F1/manuscript/tables/SupplTable1_OverviewSequencing.xlsx",
    sheet = "WGS (NovaSeq)", trim_ws = TRUE
) %>%
    # Only keep the malignant samples.
    dplyr::filter(!is.na(matched_group)) %>% 
    dplyr::mutate(
        ASID = sample, 
        sample = sprintf('%s_%s_%s', sample, strain1, strain2),
        cross = trimws(base::gsub("Reciprocal ",'', base::gsub("\\)",'', base::gsub(".*; ",'', group))))
    )


## QC - Overview ----

metadata %>%
    dplyr::inner_join(data_results$mutationalBurden) %>%
    dplyr::inner_join(data_results$flagstats) %>%
    dplyr::inner_join(
        data_results$haplotag %>% 
            dplyr::group_by(sample) %>% 
            dplyr::summarise(
                total_tagged_reads = sum(total_tagged_reads),
                total_tagged_variants = sum(total_tagged_variants)
            )
    ) %>%
    dplyr::mutate(
        sample_name = base::gsub('f1-liver-tumor-','', sample_name),
        mapped_reads = formattable::percent(Mapped / `Total reads`,2),
        paired_reads = formattable::percent(`Paired in sequencing` / `Total reads`,2),
        haplotaggable_reads = formattable::percent(`Total haplotaggable reads` / `Total reads`,2),
        percent_total_tagged_reads = formattable::percent(total_tagged_reads / `Total haplotaggable reads`,2),
        percent_h1 = formattable::percent(totalH1 / totalMutations, digits = 2),
        percent_h2 = formattable::percent(totalH2 / totalMutations, digits = 2),
        percent_unassignable = formattable::percent(totalUA / totalMutations, digits = 2),
        TMB = formattable::color_tile("grey90", max.color = "red")(round(TMB, 2)),
        total_reads = formattable::color_tile("lightpink", max.color = "hotpink")(base::formatC(`Total reads`, drop0trailing = TRUE, big.mark = ",", format = 'd'))
    ) %>%
    dplyr::select(
        ASID,
        "Cross" = cross,
        "Strain (H1)" = strain1,
        "Strain (H2)" = strain2,
        Descr. = sample_name,
        "Total reads" = total_reads,
        "Mapped" = mapped_reads,
        "Paired" = paired_reads,
        "Primary" = haplotaggable_reads,
        "Haplo-reads" = percent_total_tagged_reads,
        "Σ(SNPs)" = total_tagged_variants,
        TMB,
        "Total" = totalMutations,
        "H1" = percent_h1,
        "H2" = percent_h2,
        "UA" = percent_unassignable,
    ) %>%
    dplyr::arrange(-Total) %>%
    knitr::kable(escape = FALSE, align = "llcclccccccccccc", format.args = list(decimal.mark = ".", big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15, full_width = TRUE, html_font = "Lexend", bootstrap_options = c("striped")) %>%
    kableExtra::add_header_above(line_sep = 3, c(" " = 5, "Flagstats" = 4, "Haplotag" = 2, "Somatic variants (no.)" = 5))


# QC - Visualize the distribution of strain-assigned somatic variants. ----

data_distributions <- dplyr::bind_rows(
    base::lapply(data_somaticvariants, function(x){
        tibble::as_tibble(data.frame(sample = Biobase::sampleNames(x), chrom = seqnames(x), S4Vectors::mcols(x)[c("H1_Alt", "H2_Alt", "origin_mutant")]))
    } )
) %>%
    dplyr::inner_join(metadata) %>%
    dplyr::mutate(
        delta = H1_Alt - H2_Alt,
        origin_mutant2 = 'B6',
        origin_mutant2 = dplyr::if_else(cross == 'CAST/C3H' & origin_mutant == 'H1', 'CAST', origin_mutant2),
        origin_mutant2 = dplyr::if_else(cross == 'CAST/C3H' & origin_mutant == 'H2', 'C3H', origin_mutant2),
        origin_mutant2 = dplyr::if_else(cross == 'B6/CAST' & origin_mutant == 'H2', 'CAST', origin_mutant2),
        origin_mutant2 = dplyr::if_else(cross == 'CAST/B6' & origin_mutant == 'H2', 'CAST', origin_mutant2)
    ) %>%
    dplyr::filter(delta != 0, origin_mutant != 'UA')

# Plot distribution.
ggplot2::ggplot(data_distributions %>% dplyr::filter(!grepl('chrX', chrom)), ggplot2::aes(x = delta, fill = origin_mutant2)) +
    ggplot2::geom_histogram(binwidth = .5, na.rm = TRUE, color = "grey10", lwd = ggplot2::rel(.33)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 3500)) +
    ggplot2::scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
    ggplot2::scale_fill_manual(values = color_scheme, name = NULL) +
    ggplot2::labs(
        x = "H1 - H2 assigned reads",
        y = "Frequency<br><sup>(No. of somatic variants)</sup>"
    ) +
    ggplot2::facet_wrap(cross ~ sample_name, ncol = 5, as.table = T) +
    theme_job

# Plot distribution.
ggplot2::ggplot(data_distributions %>% dplyr::filter(sample == 'AS-614477_CAST_EiJ_C3H_HeJ'), ggplot2::aes(x = delta, fill = origin_mutant2)) +
    ggplot2::geom_histogram(binwidth = .5, na.rm = TRUE, color = "grey10", lwd = ggplot2::rel(.33)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 250)) +
    ggplot2::scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
    ggplot2::scale_fill_manual(values = color_scheme, name = NULL) +
    ggplot2::labs(
        x = "H1 - H2 assigned reads",
        y = "Frequency<br><sup>(No. of somatic variants)</sup>"
    ) +
    ggplot2::facet_wrap(~ chrom, ncol = 4) +
    theme_job
