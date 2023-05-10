# Libraries ----

library(dplyr)
library(ParallelLogger)

# Import functions. ----

source("R/functions.R")
source("R/themes.R")


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


# QC - Visualize the distribution of strain-assigned somatic variants. ----

dplyr::bind_rows(
    lapply(data_reciprocal, function(x) tibble::as_tibble(S4Vectors::mcols(x)[c("sample", "G1_Alt", "G2_Alt", "origin_mutant")]))
) %>%
    dplyr::inner_join(metadata, by = c("sample" = "ASID")) %>% 
    dplyr::mutate(
        delta = G1_Alt - G2_Alt,
        assignment = dplyr::case_when(
            origin_mutant == "G1" ~ "C57BL/6",
            origin_mutant == "G2" ~ "CAST (EiJ)",
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
        x = "C57BL/6 - CAST (EiJ) SNP-assigned reads",
        y = "Frequency<br><sup>(No. of somatic variants)</sup>"
    ) +
    ggplot2::facet_wrap(~ description, ncol = 3) +
    theme_job()

