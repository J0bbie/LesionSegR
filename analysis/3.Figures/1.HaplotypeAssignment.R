# Function: Generates a figure detailing the haplotype assignment of somatic variants (H1 vs. H2) per sample.
# Author: J. van Riet

# Load libraries. ----

library(dplyr)
library(patchwork)
library(extrafont)

# Import functions. ----

source("analysis/themes.R")


# Import data. ----

data_combined <- base::readRDS("~/odomLab/LesionSegregration_F1/data/rdata/data_combined.rds")

## QC - Overview ----

data_combined$metadata %>%
    dplyr::left_join(data_combined$tumorburden) %>%
    dplyr::inner_join(data_combined$flagstats) %>%
    dplyr::inner_join(
        data_combined$haplotyping %>% 
            dplyr::group_by(sample) %>% 
            dplyr::summarise(
                total_tagged_reads = sum(total_tagged_reads),
                total_tagged_variants = sum(total_tagged_variants)
            )
    ) %>%
    dplyr::mutate(
        mapped_reads = formattable::percent(Mapped / `Total reads`,2),
        paired_reads = formattable::percent(`Paired in sequencing` / `Total reads`,2),
        haplotaggable_reads = formattable::percent(`Total haplotaggable reads` / `Total reads`,2),
        percent_total_tagged_reads = formattable::percent(total_tagged_reads / `Total haplotaggable reads`,2),
        percent_h1 = formattable::percent(totalH1 / totalMutations, digits = 2),
        percent_h2 = formattable::percent(totalH2 / totalMutations, digits = 2),
        percent_unassignable = formattable::percent(totalUA / totalMutations, digits = 2),
        TMB = round(TMB, 2),
        total_reads = base::formatC(`Total reads`, drop0trailing = TRUE, big.mark = ",", format = 'd')
    ) %>%
    dplyr::select(
        "Sample" = sample_name,
        "Sequencing" = sequencing_type,
        "Group" = group,
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
    knitr::kable(escape = FALSE, align = "lllcclccccccccc", format.args = list(decimal.mark = ".", big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15, full_width = TRUE, html_font = "Lexend", bootstrap_options = c("striped")) %>%
    kableExtra::add_header_above(line_sep = 3, c(" " = 3, "Flagstats" = 4, "Haplotag" = 2, "Somatic variants (no.)" = 5))


# QC - Visualize the distribution of strain-assigned somatic variants. ----

data_distributions <- data_combined$somaticvariants %>% 
    dplyr::inner_join(metadata) %>%
    dplyr::mutate(
        delta = H1_Alt - H2_Alt,
        origin_mutant2 = 'B6',
        origin_mutant2 = dplyr::if_else(grepl('CAST/C3H', group) & origin_mutant == 'H1', 'CAST', origin_mutant2),
        origin_mutant2 = dplyr::if_else(grepl('CAST/C3H', group) & origin_mutant == 'H2', 'C3H', origin_mutant2),
        origin_mutant2 = dplyr::if_else(grepl('B6/CAST', group) & origin_mutant == 'H2', 'CAST', origin_mutant2),
        origin_mutant2 = dplyr::if_else(grepl('CAST/B6', group) & origin_mutant == 'H2', 'CAST', origin_mutant2),
    ) %>%
    dplyr::filter(delta != 0, origin_mutant != 'UA')

# Plot distribution.
ggplot2::ggplot(data_distributions %>% dplyr::filter(!grepl('chrX|chrY', seqnames)), ggplot2::aes(x = delta, fill = origin_mutant2)) +
    ggplot2::geom_histogram(binwidth = .5, na.rm = TRUE, color = "grey10", lwd = ggplot2::rel(.33)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 4000)) +
    ggplot2::scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
    ggplot2::scale_fill_manual(values = color_scheme, name = NULL) +
    ggplot2::labs(
        x = "H1 - H2 assigned reads",
        y = "Frequency<br><sup>(No. of somatic variants)</sup>"
    ) +
    ggplot2::facet_wrap(group ~ sample_name, ncol = 5, as.table = T) +
    theme_job

# Plot distribution - Per chromosome.
ggplot2::ggplot(data_distributions %>% dplyr::filter(sample == 'AS-1132164'), ggplot2::aes(x = delta, fill = origin_mutant2)) +
    ggplot2::geom_histogram(binwidth = .5, na.rm = TRUE, color = "grey10", lwd = ggplot2::rel(.33)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 400)) +
    ggplot2::scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
    ggplot2::scale_fill_manual(values = color_scheme, name = NULL) +
    ggplot2::labs(
        x = "H1 - H2 assigned reads",
        y = "Frequency<br><sup>(No. of somatic variants)</sup>"
    ) +
    ggplot2::facet_wrap(~ seqnames, ncol = 4) +
    theme_job
