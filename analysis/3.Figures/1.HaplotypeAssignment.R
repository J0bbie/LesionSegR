# Function: Generates a figure detailing the haplotype assignment of somatic variants (G1 vs. G2) per sample.
# Author: J. van Riet

# Load libraries. ----

library(dplyr)
library(patchwork)
library(extrafont)

# Import functions. ----

source("analysis/functions.R")
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
        sample = sprintf('%s_%s_%s', sample, strain1, strain2)
    )


## QC - SNPsplit ----

metadata %>%
    dplyr::inner_join(data_results$mutationalBurden) %>%
    dplyr::inner_join(data_results$flagstats, by = c("ASID" = "sample")) %>%
    dplyr::mutate(
        mapped_reads = formattable::percent(Mapped / `Total reads`,2),
        paired_reads = formattable::percent(`Paired in sequencing` / `Total reads`,2),
        percent_h1 = formattable::percent(totalH1 / totalMutations, digits = 2),
        percent_h2 = formattable::percent(totalH2 / totalMutations, digits = 2),
        percent_unassignable = formattable::percent(totalUA / totalMutations, digits = 2),
        TMB = formattable::color_tile("grey90", max.color = "red")(round(TMB, 2)),
        total_reads = formattable::color_tile("lightpink", max.color = "hotpink")(base::formatC(.$`Total reads`, drop0trailing = TRUE, big.mark = ","))
    ) %>%
    dplyr::select(
        ASID,
        "Strain (H1)" = strain1,
        "Strain (H2)" = strain2,
        Descr. = sample_name,
        "Total reads" = total_reads,
        "Mapped" = mapped_reads,
        "Paired" = paired_reads,
        TMB,
        "Total" = totalMutations,
        "H1" = percent_h1,
        "H2" = percent_h2,
        "UA" = percent_unassignable,
    ) %>%
    dplyr::arrange(-Total) %>%
    knitr::kable(escape = FALSE, align = "lcclcccccccc", format.args = list(decimal.mark = ".", big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15, full_width = TRUE, html_font = "Lexend", bootstrap_options = c("striped")) %>%
    kableExtra::add_header_above(line_sep = 10, c(" " = 4, "Flagstats" = 3, "Somatic variants (no.)" = 5))


# QC - Visualize the distribution of strain-assigned somatic variants. ----

dplyr::bind_rows(
    base::lapply(data_somaticvariants, function(x) tibble::as_tibble(data.frame(sample = Biobase::sampleNames(x), S4Vectors::mcols(x)[c("H1_Alt", "H2_Alt", "origin_mutant")])))
) %>%
    dplyr::inner_join(metadata) %>%
    dplyr::mutate(
        delta = H1_Alt - H2_Alt,
    ) %>%
    dplyr::filter(delta != 0) %>%
    # Plot distribution.
    ggplot2::ggplot(., ggplot2::aes(x = delta, fill = origin_mutant)) +
    ggplot2::geom_histogram(binwidth = .5, na.rm = TRUE, color = "grey10", lwd = ggplot2::rel(.33)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 4100)) +
    ggplot2::scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
    ggplot2::scale_fill_manual(values = color_scheme, name = NULL) +
    ggplot2::labs(
        x = "H1 - H2 assigned reads",
        y = "Frequency<br><sup>(No. of somatic variants)</sup>"
    ) +
    ggplot2::facet_wrap(~ sample_name, ncol = 3) +
    theme_job
