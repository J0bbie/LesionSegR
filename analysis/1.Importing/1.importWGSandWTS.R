# Function: Import the WGS data of the F1 mice.
# Author: J. van Riet

# Libraries ----

library(dplyr)
library(future)
library(LesionSegR)
library(VariantAnnotation)

# Parallel settings.
future::plan(future::multisession, workers = 8)

# Import metadata. ----

metadata <- readr::read_tsv("/omics/groups/OE0538/internal/users/e480l/projects/DEN_tumors/snakemake/TEST_files_Novaseq/DNA/TEST_38558_Novaseq_DNA_samplesheet.tsv", show_col_types = FALSE) %>%
    dplyr::bind_rows(readr::read_tsv("/omics/groups/OE0538/internal/users/e480l/projects/DEN_tumors/snakemake/TEST_files_Novaseq/RNA/TEST_38415_Novaseq_RNA_sample_sheet.tsv", show_col_types = FALSE)) %>% 
    dplyr::mutate(
        sample = sequencing_name,
        sample_name = tolower(sample_name),
        sample_strain = paste(sample_name, strain1, strain2, sep = '_'),
        seqname_strain = paste(sequencing_name, strain1, strain2, sep = '_')
    )

workflow_dir <- "/omics/odcf/analysis/OE0538_projects/DO-0006/f1_b6_mcas/e480l/projects/DEN_tumors/TEST_same_folders/"

# Subset GTF on genes to analyze. ----

gtf <- rtracklayer::import("/omics/groups/OE0538/internal/projects/sharedData/GRCm39/annotation/gencode.vM34.basic.annotation.gtf", format = "gff") %>%
    tibble::as_tibble() %>%
    dplyr::filter(
        type == "gene",
        !base::is.na(gene_id)
    ) %>% 
    # Subset on classes-of-interest.
    dplyr::filter(
        base::grepl("protein_coding|lncRNA|IG_.*_gene", gene_type)
    ) %>% 
    # Additional filtering.
    dplyr::filter(
        !base::grepl("^Gm[0-9]", gene_name) | gene_name == 'Gm2a'
    ) %>%
    # Remove RIKEN lncRNA genes.
    dplyr::filter(
        ! (base::grepl("Rik$", gene_name) & gene_type == 'lncRNA')
    ) %>%
    # Remove genes without a gene-name.
    dplyr::filter(!grepl("ENSMUS", gene_name)) %>% 
    dplyr::distinct(
        seqnames, start, end, strand, gene_id, gene_type, gene_name
    ) %>%
    dplyr::filter(!duplicated(gene_name))

## Import WGS and WTS data. ----

data_combined <- import_samples(metadata, workflow_dir, gtf = gtf)

## Mutational signatures. ----

data_combined$signatures_sbs <- sigminer::sig_fit(data_combined$mutmatrices_sbs, sig_index = 'ALL', sig_db = 'SBS_mm10', auto_reduce = T, type = "relative")
sigminer::show_catalogue(data_combined$mutmatrices_sbs, mode = "SBS", style = "cosmic", x_label_angle = 90, samples = colnames(data_combined$mutmatrices_sbs))

## Driver analysis: dN/dS. ----

data_combined$dNdS <- LesionSegR::run_dnds(data_combined$somaticvariants, path_db = "/omics/groups/OE0538/internal/projects/sharedData/GRCm39/annotation/refCDS_ENSEMBLv110_GRCm39.rda")
data_combined$dNdS <- data_combined$dNdS$finalOutput

# Save the data. ----

data_combined$metadata <- metadata
saveRDS(data_combined, "~/odomLab/LesionSegregration_F1/data/rdata/data_combined.rds")
