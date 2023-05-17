# Function: Import the WGS data of the F1 mice.
# Author: J. van Riet

# Libraries ----

library(dplyr)
library(ParallelLogger)
library(patchwork)

# Import functions. ----

source("R/functions.R")
source("R/themes.R")

# Import somatic variants and CNV. ----

## Import somatic variants. ----
files_vcf <- list.files(path = "~/Downloads/VCF/", pattern = "*withStrainCounts.vcf.gz$", full.names = TRUE)
data_somaticvariants <- pbapply::pblapply(files_vcf, import_vcf, cl = 4)


## Report and remove litter-specific artifacts. ----

recurrent_mutations <- dplyr::bind_rows(base::lapply(data_somaticvariants, function(x) {
    tibble::as_tibble(
        data.frame(
            chr = GenomeInfoDb::seqnames(x),
            start = GenomicRanges::start(x),
            end = GenomicRanges::end(x),
            ref = VariantAnnotation::ref(x),
            alt = VariantAnnotation::alt(x),
            sample = unique(x$sample)
        )
    )
})) %>%
    dplyr::left_join(metadata, by = c("sample" = "ASID")) %>%
    dplyr::group_by(chr, start, end, ref, alt) %>%
    dplyr::summarise(
        countSamples = dplyr::n_distinct(sample),
        totalCasB6 = sum(Type == "Tumor (Liver; Reciprocal CAS/B6)"),
        totalB6Cas = sum(Type == "Tumor (Liver; Reciprocal B6/CAS)")
    ) %>%
    dplyr::ungroup()


# Plot no. of shared cross-specific mutations. ----

reshape2::acast(
    recurrent_mutations, totalCasB6 ~ totalB6Cas,
    value.var = "countSamples",
    fun.aggregate = length
) %>%
    pheatmap::pheatmap(.,
        cluster_rows = FALSE, cluster_cols = FALSE,
        fontsize = 10, cellwidth = 30, cellheight = 30,
        border_color = "grey25", scale = "none", display_numbers = TRUE, fontsize_number = 8, number_format = "%.0f", na_col = "white", na_rm = TRUE,
        color = "white", legend = FALSE
    )


## Remove litter-specific artifacts. ----

data_somaticvariants <- lapply(data_somaticvariants, function(x) {
    if (unique(x$sample) %in% metadata$ASID[metadata$Type == "Tumor (Liver; Reciprocal CAS/B6)"]) {
        remove <- subset(recurrent_mutations, totalCasB6 == 6)
    } else {
        remove <- subset(recurrent_mutations, totalB6Cas == 6)
    }

    remove <- paste0(remove$chr, remove$start, remove$ref, remove$alt)
    x$id <- paste0(as.character(GenomeInfoDb::seqnames(x)), GenomicRanges::start(x), x@ref, x@alt)

    # Remove variants.
    x <- x[!x$id %in% remove, ]
    x$id <- NULL

    return(x)
})


# Sample-wise comparison of shared mutations. ----

count_overlap <- function(set1, set2) {
    x <- paste0(GenomeInfoDb::seqnames(set1), set1@ranges, set1@ref, set1@alt)
    y <- paste0(GenomeInfoDb::seqnames(set2), set2@ranges, set2@ref, set2@alt)
    return(base::length(intersect(x, y)))
}

samplewise_intersect <- dplyr::bind_rows(base::lapply(seq_along(data_somaticvariants), function(i) {
    dplyr::bind_rows(base::lapply(seq_along(data_somaticvariants), function(j) {
        is <- metadata %>%
            dplyr::filter(ASID == unique(data_somaticvariants[[i]]$sample)) %>%
            dplyr::pull(description)
        js <- metadata %>%
            dplyr::filter(ASID == unique(data_somaticvariants[[j]]$sample)) %>%
            dplyr::pull(description)

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


# Convert to VRangesList. ----

data_somaticvariants <- VariantAnnotation::VRangesList(data_somaticvariants)
names(data_somaticvariants) <- vapply(data_somaticvariants, function(x) unique(x$sample), "")


## Determine WGS characteristics. ----

reciprocal_results <- list()

## Flagstats. ----

read_flagstats <- function(x) {
    readr::read_csv(x, col_names = FALSE, show_col_types = FALSE) %>%
        dplyr::mutate(
            X1 = base::trimws(X1),
            X2 = as.integer(base::trimws(gsub(":.*", "", X2))),
            sample = basename(x) %>% stringr::str_replace_all(pattern = "_.*", replacement = "")
            ) %>%
        dplyr::select(variable = X1, value = X2, sample)
}

files_flagstats <- list.files(path = "~/DKFZ/odomLab/LesionSegregration_F1/data/WGS/alignment/B6_CAST-EiJ/", pattern = "*.flagstats$", full.names = TRUE)
data_flagstats <- tibble::as_tibble(dplyr::bind_rows(base::lapply(files_flagstats, read_flagstats)) %>% reshape2::dcast(., sample ~ variable))

reciprocal_results$flagstats <- data_flagstats

## SNPSplit. ----

# Read yaml files.
read_snpsplit <- function(x) {
    tibble::as_tibble(yaml::read_yaml(x)$Tagging) %>%
        dplyr::mutate(
            sample = basename(x) %>% stringr::str_replace_all(pattern = "_.*", replacement = "")
        )
}

files_snpsplit <- list.files(path = "~/DKFZ/odomLab/LesionSegregration_F1/data/WGS/alignment/B6_CAST-EiJ/", pattern = "*.SNPsplit_report.yaml$", full.names = TRUE)
data_snpsplit <- dplyr::bind_rows(base::lapply(files_snpsdataplit, read_snpsplit))
reciprocal_results$snpsplit <- data_snpsplit %>% dplyr::inner_join(metadata, by = c("sample" = "ASID"))


## Determine mutational burden. ----
reciprocal_results$mutationalBurden <- dplyr::bind_rows(base::lapply(data_somaticvariants, determine_mutational_burden))

## Mutational signature analysis. ----
reciprocal_results$mutSigs <- determine_knownsignatures(base::unlist(data_somaticvariants))

## dN/dS. ----

reciprocal_results$dNdS <- NULL

## Determine no. of Ti/Tv and ratio. -----

reciprocal_results$TiTv <- NULL


# Save the data. ----

saveRDS(data_somaticvariants, "data/reciprocal_somaticvariants.rds")
saveRDS(reciprocal_results, "data/reciprocal_results.rds")


# Close logger ------------------------------------------------------------

ParallelLogger::unregisterLogger()


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
