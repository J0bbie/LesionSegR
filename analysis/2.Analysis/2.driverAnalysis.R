library(LesionSegR)
library(dplyr)

# Load metadata. ----

metadata <- readxl::read_xlsx(
    path = "~/Dropbox/Work/DKFZ/Projects/Odom - LesionSegregation/Draft/Tables/SupplTable1_OverviewSequencing.xlsx",
    sheet = "WGS_Reciprocals", trim_ws = TRUE
) %>%
    # Only keep the malignant samples.
    dplyr::filter(!is.na(matched_group))

# Generate overview of known driver genes. ----
data_somaticvariants <- base::readRDS("data/reciprocal_somaticvariants.rds")
data_somaticvariants <- base::unlist(data_somaticvariants)


known_drivers <- subset(data_somaticvariants, SYMBOL %in% c("Kras", "Hras", "Egfr", "Braf", "Stag2", "Pilra", "Slc43a2", "Brix1", "Rps18-ps6") & ConsequenceAll == "missense_variant")

known_drivers <- tibble::as_tibble(known_drivers) %>% dplyr::inner_join(metadata)


table(known_drivers$group, known_drivers$SYMBOL)

results <- base::readRDS("data/reciprocal_results.rds")
results$mutationalBurden
