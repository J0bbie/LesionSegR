# Function: Generate an overview of the sequenced mice.
# Author: J. van Riet

# Libraries ----

library(dplyr)
library(waffle)

# Import metadata. ----

metadata_wgs <- readxl::read_xlsx(
    path = "~/Dropbox/Work/DKFZ/Projects/Odom - LesionSegregation/Draft/Tables/SupplTable1_OverviewSequencing.xlsx",
    sheet = "WGS", trim_ws = TRUE
) %>% 
# Only take along the latest sequencing run on the NovaSeq.
dplyr::filter(sequencingPlatform == "Illumina NovaSeq 6000") %>% 
dplyr::mutate(miceID = gsub(" \\(.*", "", miceID))

metadata_wts <- readxl::read_xlsx(
    path = "~/Dropbox/Work/DKFZ/Projects/Odom - LesionSegregation/Draft/Tables/SupplTable1_OverviewSequencing.xlsx",
    sheet = "WTS", trim_ws = TRUE
)

# Number of samples with matching WGS / WTS.
metadata_wgs <- metadata_wgs %>% dplyr::mutate(hasWTS = description %in% metadata_wts$description)
metadata_wts <- metadata_wts %>% dplyr::mutate(hasWGS = description %in% metadata_wgs$description)
