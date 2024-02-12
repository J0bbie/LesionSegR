library(LesionSegR)
library(dplyr)

# Load metadata. ----

metadata <- readr::read_tsv("~/jvanriet/git/snakemake-lesionsegregation/workflow/examples/example_samplesheet.tsv", show_col_types = FALSE) %>%
    dplyr::mutate(
        sample_strain = paste(sample_name, strain1, strain2, sep = '_'),
        seqname_strain = paste(sequencing_name, strain1, strain2, sep = '_'),
    )

# Import data. ----
data <- base::readRDS("~/odomLab/LesionSegregration_F1/data/rdata/data_combined.rds")

# Generate overview of known driver genes. ----
known_drivers <- subset(data_somaticvariants, SYMBOL %in% c("Kras", "Hras", "Egfr", "Braf", "Stag2", "Pilra", "Slc43a2", "Brix1", "Rps18-ps6") & ConsequenceAll == "missense_variant")
known_drivers <- tibble::as_tibble(known_drivers) %>% dplyr::inner_join(metadata)

table(known_drivers$group, known_drivers$SYMBOL)