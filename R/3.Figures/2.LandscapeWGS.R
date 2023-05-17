# Function: Generates a landscape figure of the WGS samples.
# Author: J. van Riet

# Load libraries. ----

library(dplyr)
library(ParallelLogger)
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

# Sorting of samples. ----

order_samples <- reciprocal_results$mutationalBurden %>%
    dplyr::inner_join(metadata, by = c("sample" = "ASID")) %>%
    dplyr::arrange(Type, -TMB) %>%
    dplyr::pull(sample)


# Generation of tracks. ----

tracks <- list()

## Track - Strain. ----

metadata %>%
    dplyr::distinct(ASID, Type) %>%
    dplyr::mutate(sampleId = factor(ASID, levels = order_samples)) %>%
    ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = "Tissue", fill = Type)) +
    ggplot2::geom_tile(width = .9, colour = "grey25", lwd = .5, na.rm = TRUE) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = color_scheme, guide = ggplot2::guide_legend(title = NULL, title.position = "top", title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    theme_anno_job

## Track - Genome-wide TMB. ---

reciprocal_results$mutationalBurden %>%
    dplyr::mutate(sampleId = factor(sample, levels = order_samples)) %>%
    ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = TMB)) +
    ggplot2::geom_segment(ggplot2::aes(xend = sampleId, yend = 0), lwd = 1, lty = "dotted", color = "grey25") +
    ggplot2::geom_point(shape = 21, fill = "black", color = "black", size = 4) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30, 5), limits = c(0, 31)) +
    ggplot2::geom_hline(yintercept = 10, color = "black", lty = "dotted", lwd = .5) +
    ggplot2::labs(y = "Tumor Mutational Burden") +
    theme_anno_job

## Track - Haplotype assignment. ----

reciprocal_results$mutationalBurden %>%
    reshape2::melt(., id.vars = c("sample", "totalMutations"), measure.vars = c("totalG1", "totalG2", "totalUA")) %>%
    dplyr::mutate(
        sampleId = factor(sample, levels = order_samples),
        value_perc = value / totalMutations,
        variable = factor(variable, levels = c("totalG1", "totalG2", "totalUA"), labels = c("C57BL/6", "CAST (EiJ)", "Unassigned"))
    ) %>%
    ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = value_perc, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = .9, colour = "grey25", lwd = .5) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_fill_manual(values = color_scheme, name = NULL) +
    ggplot2::labs(y = "Haplotype assignment") +
    theme_anno_job

## Track - Genome-wide ploidy. ----


## Track - Mutational patterns. ----

### Track - Mutational signatures (SBS). ----

### Track - Mutational signatures (DBS). ----

### Track - Mutational signatures (ID). ----
