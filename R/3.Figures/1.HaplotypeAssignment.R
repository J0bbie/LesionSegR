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
    sheet = "WGS_Reciprocals", trim_ws = TRUE
)


## QC - SNPsplit ----

reciprocal_results$snpsplit %>%
    dplyr::inner_join(metadata) %>%
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
        Descr. = sample_name,
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
    knitr::kable(escape = FALSE, align = "llcccccccccc", format.args = list(decimal.mark = ".", big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15, full_width = TRUE, html_font = "Lexend", bootstrap_options = c("striped")) %>%
    kableExtra::add_header_above(line_sep = 10, c(" " = 3, "Flagstats" = 3, "SNPSplit (% reads)" = 4, "Somatic variants (no.)" = 5))


# QC - Visualize the distribution of strain-assigned somatic variants. ----

dplyr::bind_rows(
    base::lapply(reciprocal_data, function(x) tibble::as_tibble(S4Vectors::mcols(x)[c("sample", "G1_Alt", "G2_Alt", "origin_mutant")]))
) %>%
    dplyr::inner_join(metadata) %>%
    dplyr::mutate(
        delta = G1_Alt - G2_Alt,
        assignment = dplyr::case_when(
            origin_mutant == "G1" ~ "B6",
            origin_mutant == "G2" ~ "CAST",
            origin_mutant == "UA" ~ "Unassigned"
        )
    ) %>%
    dplyr::filter(delta != 0) %>%
    # Plot distribution.
    ggplot2::ggplot(., ggplot2::aes(x = delta, fill = assignment)) +
    ggplot2::geom_histogram(binwidth = .5, na.rm = TRUE, color = "grey10", lwd = ggplot2::rel(.33)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 4100)) +
    ggplot2::scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
    ggplot2::scale_fill_manual(values = color_scheme, name = NULL) +
    ggplot2::labs(
        x = "B6 - CAST SNP-assigned reads",
        y = "Frequency<br><sup>(No. of somatic variants)</sup>"
    ) +
    ggplot2::facet_wrap(~ sample_name, ncol = 3) +
    theme_job
