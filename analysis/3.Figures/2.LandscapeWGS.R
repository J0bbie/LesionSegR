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
    sheet = "WGS_Reciprocals", trim_ws = TRUE
) %>%
    dplyr::filter(!is.na(matched_group))

# Sorting of samples. ----

order_samples <- reciprocal_results$mutationalBurden %>%
    dplyr::inner_join(metadata) %>%
    dplyr::arrange(group, -TMB) %>%
    dplyr::pull(sample)


# Generation of tracks. ----

tracks <- list()


## Track - Genome-wide TMB. ---

tracks$tmb <- reciprocal_results$mutationalBurden %>%
    dplyr::mutate(sampleId = factor(sample, levels = order_samples)) %>%
    ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = TMB)) +
    ggplot2::geom_segment(ggplot2::aes(xend = sampleId, yend = 0), lwd = 1, color = "grey25") +
    ggplot2::geom_point(shape = 21, fill = "black", color = "black", size = 3) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30, 5), limits = c(0, 31)) +
    ggplot2::geom_hline(yintercept = 10, color = "black", lty = "dotted", lwd = .5) +
    ggplot2::labs(y = "TMB") +
    theme_anno_job

## Track - Haplotype assignment. ----

tracks$haplotype <- reciprocal_results$mutationalBurden %>%
    reshape2::melt(., id.vars = c("sample", "totalMutations"), measure.vars = c("totalG1", "totalG2", "totalUA")) %>%
    dplyr::mutate(
        sampleId = factor(sample, levels = order_samples),
        value_perc = value / totalMutations,
        variable = factor(variable, levels = c("totalG1", "totalG2", "totalUA"), labels = c("Unassigned", "B6", "CAST"))
    ) %>%
    ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = value_perc, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = .9, colour = "grey25", lwd = .5) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_fill_manual(values = color_scheme, name = NULL) +
    ggplot2::labs(y = "Haplotype<br><sub>(Somatic variants)</sub>") +
    theme_anno_job

## Track - Mice ID. ----

tracks$mice <- metadata %>%
    dplyr::mutate(sampleId = factor(sample, levels = order_samples)) %>%
    ggplot2::ggplot(., ggplot2::aes(x = sample, y = "Mice", fill = mice_id)) +
    ggplot2::geom_tile(width = .9, colour = "grey25", lwd = .5, na.rm = TRUE) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(name = 'Dark2', n = 8), guide = ggplot2::guide_legend(title = NULL, title.position = "top", title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    theme_anno_job

## Track - Strain. ----

tracks$strain <- metadata %>%
    dplyr::inner_join(reciprocal_results$mutationalBurden) %>%
    dplyr::distinct(sample, sample_name, group) %>%
    dplyr::mutate(sampleId = factor(sample, levels = order_samples)) %>%
    ggplot2::ggplot(., ggplot2::aes(x = sample_name, y = "Tissue", fill = group)) +
    ggplot2::geom_tile(width = .9, colour = "grey25", lwd = .5, na.rm = TRUE) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = color_scheme, guide = ggplot2::guide_legend(title = NULL, title.position = "top", title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    theme_anno_job +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, vjust = .5, size = 10, face = "bold"))

### Combine tracks.


tracks$tmb +
tracks$haplotype +
tracks$mice +
tracks$strain +
patchwork::plot_layout(ncol = 1, guides = "collect", heights = c(1, 1, .2, .2)) +
patchwork::plot_annotation(tag_levels = "a")
