# Libraries ----

library(dplyr)


# Import functions. ----

source("R/functions.R")

# Read in VCF. ----

path_vcf <- "~/Downloads/VCF/AS-949280_withStrainCounts.vcf.gz"
x <- import_vcf(path_vcf)

# Determine allelic origin. ----

# Compare methods.
test <- tibble::as_tibble(x) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        ratio = determine_origin_ratio(G1, G2),
        chisq = determine_origin_chisq(G1, G2)
    )

table(test$ratio)
table(test$chisq)

knitr::kable(
    table(Ratio = test$ratio, Chisq = test$chisq),
    align = "l", caption = "Comparison of methods to determine allelic origin.",
    booktabs = TRUE
)

test[test$ratio != test$chisq, ]

# Plot distribution of G1.
test %>%
    dplyr::filter(ratio != "AU") %>%
    ggplot2::ggplot(ggplot2::aes(x = G1, fill = ratio)) +
    ggplot2::geom_histogram(binwidth = 1) +
    ggplot2::scale_fill_manual(values = c("G1" = "red", "G2" = "blue")) +
    ggplot2::theme_bw()
