# Function: Import the WGS data of the F1 mice.
# Author: J. van Riet

# Libraries ----

library(dplyr)
library(patchwork)

# Import functions. ----

source("analysis/functions.R")
source("analysis/themes.R")

# Import somatic variants. ----

data_somaticvariants <- base::readRDS("~/odomLab/LesionSegregration_F1/data/rdata/data_somaticvariants.rds")


## Report potential cross-specific artifacts. ----

recurrent_mutations <- dplyr::bind_rows(base::lapply(data_somaticvariants, function(x) {
    tibble::as_tibble(
        data.frame(
            chr = GenomeInfoDb::seqnames(x),
            start = GenomicRanges::start(x),
            end = GenomicRanges::end(x),
            ref = VariantAnnotation::ref(x),
            alt = VariantAnnotation::alt(x),
            sample = Biobase::sampleNames(x)
        )
    )
})) %>%
    dplyr::inner_join(metadata) %>%
    dplyr::group_by(chr, start, end, ref, alt) %>%
    dplyr::summarise(
        countSamples = dplyr::n_distinct(sample),
        totalCASTC3H = sum(group == 'Tumor (Liver; CAST/C3H)'),
        totalCasB6 = sum(group == "Tumor (Liver; Reciprocal CAST/B6)"),
        totalB6Cas = sum(group == "Tumor (Liver; Reciprocal B6/CAST)")
    ) %>%
    dplyr::ungroup()


# Plot no. of shared cross-specific mutations. ----

recurrent_mutations %>%
    dplyr::group_by(totalCasB6, totalB6Cas, totalCASTC3H) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(totalCasB6, totalB6Cas, fill = list(count = 0)) %>%
    ggplot2::ggplot(ggplot2::aes(x = totalCasB6, y = totalB6Cas, fill = count, label = count)) +
    ggplot2::geom_tile(color = "black", lwd = 0.5) +
    ggplot2::geom_text(size = 5) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_gradient(low = "white", high = "#F70289") +
    ggplot2::labs(
        x = "No. of shared somatic mutations<br>CAS/B6",
        y = "No. of shared somatic mutations<br>B6/CAS",
    ) +
    ggplot2::theme(
        legend.position = "none"
    )

# # Remove litter-specific artifacts.

# b6cas_mice <- metadata %>%
#     dplyr::filter(group == "Tumor (Liver; Reciprocal B6/CAS)") %>%
#     dplyr::pull(sample)

# data_somaticvariants <- VariantAnnotation::VRangesList(lapply(data_somaticvariants, function(x) {
#     if (unique(x$sample) %in% b6cas_mice) {
#         remove <- recurrent_mutations %>%
#             dplyr::filter(totalB6Cas == 6)
#     } else {
#         remove <- recurrent_mutations %>%
#             dplyr::filter(totalCasB6 == 6)
#     }

#     IRanges::subsetByOverlaps(x, GenomicRanges::makeGRangesFromDataFrame(remove), invert = TRUE)
# }))

# z = data_somaticvariants$`AS-949280`
# z = IRanges::subsetByOverlaps(z, GenomicRanges::makeGRangesFromDataFrame(recurrent_mutations %>% dplyr::filter(totalB6Cas == 6)), invert = FALSE)


# Sample-wise comparison of shared mutations. ----

count_overlap <- function(set1, set2) {
    x <- paste0(GenomeInfoDb::seqnames(set1), set1@ranges, set1@ref, set1@alt)
    y <- paste0(GenomeInfoDb::seqnames(set2), set2@ranges, set2@ref, set2@alt)
    return(base::length(intersect(x, y)))
}

samplewise_intersect <- dplyr::bind_rows(base::lapply(seq_along(data_somaticvariants), function(i) {
    dplyr::bind_rows(base::lapply(seq_along(data_somaticvariants), function(j) {
        is <- metadata %>%
            dplyr::filter(sample == unique(data_somaticvariants[[i]]$sample)) %>%
            dplyr::pull(sample_name)
        js <- metadata %>%
            dplyr::filter(sample == unique(data_somaticvariants[[j]]$sample)) %>%
            dplyr::pull(sample_name)

        if (i != j) {
            return(tibble::as_tibble(data.frame(i = is, j = js, shared = count_overlap(data_somaticvariants[[i]], data_somaticvariants[[j]]))))
        } else {
            return(tibble::as_tibble(data.frame(i = is, j = js, shared = NA)))
        }
    }))
}))

# Plot sample-wise comparison.
reshape2::acast(samplewise_intersect, i ~ j, value.var = "shared") %>%
    pheatmap::pheatmap(.,
        cluster_rows = TRUE, cluster_cols = TRUE, treeheight_row = 50, treeheight_col = 50,
        fontsize = 10, cellwidth = 25, cellheight = 25,
        border_color = "grey25", scale = "none", display_numbers = TRUE, fontsize_number = 7, number_format = "%.0f", na_col = "white", na_rm = TRUE,
        color = grDevices::colorRampPalette(c("grey99", "hotpink"))(50)
    )
