# Function: Generate an overview of the WGS landscape of the F1-mice.
# Author: J. van Riet

# Libraries ----

library(dplyr)
library(ParallelLogger)

# Import functions. ----

source("R/functions.R")

# Import metadata. ----

metadata <- readxl::read_xlsx(
    path = "~/DKFZ/odomLab/LesionSegregration_F1/metadata/SupplTable1_OverviewSequencing.xlsx",
    sheet = "WGS", trim_ws = TRUE
)

# Retrieve the VCF files.

files_vcf <- list.files(path = "~/Downloads/VCF/", pattern = "*withStrainCounts.vcf.gz$", full.names = TRUE)


# Read in somatic variants. ----

data_reciprocal <- lapply(files_vcf[1:2], import_vcf)

# Determine tumor mutational burden.
mutational_info <- dplyr::bind_rows(lapply(data_reciprocal, calculate_tmb))

# Perform mutational signature analysis.
data_mutsigs <- determine_knownsignatures(data_reciprocal[[2]])
MutationalPatterns::plot_96_profile(data_mutsigs$sbs$mutMatrix, condensed = TRUE, ymax = 0.05)
MutationalPatterns::plot_contribution(data_mutsigs$sbs$fit_res$contribution, data_mutsigs$sbs$signature, mode = "relative")



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
    ggplot2::facet_grid(~ CONTIG, space = "free") +
    ggplot2::scale_x_continuous(labels = scales::comma)
