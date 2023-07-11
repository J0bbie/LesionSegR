# Function: Import the WGS data of the F1 mice.
# Author: J. van Riet

# Libraries ----

library(dplyr)
library(ParallelLogger)
library(patchwork)

# Import functions. ----

source("R/functions.R")
source("R/themes.R")

# Import somatic variants. ----

reciprocal_data <- base::readRDS("data/reciprocal_somaticvariants.rds")
reciprocal_results <- base::readRDS("data/reciprocal_results.rds")


# Generate the 96-context matrices. ----

mut_matrixes_96 <- generate_mutmatrices_96(reciprocal_data)

# Visualize the 96-context matrices. ----
sigminer::show_catalogue(reciprocal_results$mut_matrixes_96$sbs, mode = "SBS", style = "cosmic", samples = c("AS-949282", "AS-949280"))

# Perform custom signature analysis (SBS). ----

## Run bootstrapped NMF analysis. ----
analysis_sbs <- sigminer::bp_extract_signatures(t(mut_matrixes_96$sbs), n_bootstrap = 5, n_nmf_run = 5, range = 2:8, cores = 10)

## Determine the optimal number of signatures. ----
sigminer::bp_show_survey2(reciprocal_results$mutsigs_sbs, highlight = 5)

## Plot the signatures. ----
analysis_sbs_sigs <- sigminer::bp_get_sig_obj(analysis_sbs, signum = 5)
sigminer::show_sig_profile(analysis_sbs_sigs, mode = "SBS", style = "cosmic")

sigminer::show_sig_exposure(analysis_sbs_sigs, style = "cosmic", hide_samps = FALSE)

sim <- sigminer::get_sig_similarity(analysis_sbs_sigs, sig_db = "SBS_mm10")
pheatmap::pheatmap(sim$similarity)
