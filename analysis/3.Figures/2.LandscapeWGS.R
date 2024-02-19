# Function: Generates a landscape figure of the WGS samples.
# Author: J. van Riet

# Load libraries. ----

library(dplyr)
library(patchwork)
library(extrafont)

# Import functions. ----

source("analysis/themes.R")


# Import data. ----

data_combined <- base::readRDS("~/odomLab/LesionSegregration_F1/data/rdata/data_combined.rds")

# Sorting of samples. ----

order_samples <- data_combined$tumorburden %>%
    dplyr::inner_join(metadata) %>%
    dplyr::arrange(group, -TMB) %>%
    dplyr::pull(sample)


# Generation of tracks. ----

tracks <- list()


## Track - Genome-wide TMB. ---

tracks$tmb <- data_combined$tumorburden %>%
    dplyr::mutate(sampleId = factor(sample, levels = order_samples)) %>%
    ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = TMB)) +
    ggplot2::geom_segment(ggplot2::aes(xend = sampleId, yend = 0), lwd = 1, color = "grey25") +
    ggplot2::geom_point(shape = 21, fill = "black", color = "black", size = 3) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30, 5), limits = c(0, 31)) +
    ggplot2::geom_hline(yintercept = 10, color = "black", lty = "dotted", lwd = .5) +
    ggplot2::labs(y = "TMB") +
    theme_anno_job

## Track - Haplotype assignment. ----

tracks$haplotype <- data_combined$tumorburden %>%
    dplyr::inner_join(metadata) %>% 
    reshape2::melt(., id.vars = c("sample", "totalMutations", 'group'), measure.vars = c("totalH1", "totalH2", "totalUA")) %>%
    dplyr::mutate(
        variable = dplyr::if_else(variable == 'totalH1' & grepl('CAST/C3H', group), 'CAST', variable),
        variable = dplyr::if_else(variable == 'totalH2' & grepl('CAST/C3H', group), 'C3H', variable),
        variable = dplyr::if_else(variable == 'totalH1' & !grepl('CAST/C3H', group), 'B6', variable),
        variable = dplyr::if_else(variable == 'totalH2' & !grepl('CAST/C3H', group), 'CAST', variable),
        variable = dplyr::if_else(variable == 'totalUA', 'Unassigned', variable)
    ) %>% 
    dplyr::mutate(
        sampleId = factor(sample, levels = order_samples),
        value_perc = value / totalMutations,
        variable = factor(variable, levels = c("B6", "CAST", "C3H", "Unassigned"))
    ) %>%
    ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = value_perc, fill = variable)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = .9, colour = "grey25", lwd = .5) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_fill_manual(values = color_scheme, name = NULL) +
    ggplot2::labs(y = "Haplotype<br><sub>(Somatic variants)</sub>") +
    theme_anno_job


## Track - Strain. ----

tracks$strain <- data_combined$metadata %>%
    dplyr::distinct(sample, sample_name, group) %>%
    dplyr::filter(sample %in% order_samples) %>% 
    dplyr::mutate(sampleId = factor(sample, levels = order_samples)) %>%
    ggplot2::ggplot(., ggplot2::aes(x = sample_name, y = "Tissue", fill = group)) +
    ggplot2::geom_tile(width = .9, colour = "grey25", lwd = .5, na.rm = TRUE) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = color_scheme, guide = ggplot2::guide_legend(title = NULL, title.position = "top", title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    theme_anno_job +
    ggplot2::theme(
        legend.position = 'bottom',
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = .5, size = 10, face = "bold")
    )

### Combine tracks.

tracks$tmb +
    tracks$haplotype +
    tracks$strain +
    patchwork::plot_layout(ncol = 1, guides = "keep", heights = c(.33, 1, .1)) +
    patchwork::plot_annotation(tag_levels = "a")
