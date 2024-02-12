# Function: Import the WGS data of the F1 mice.
# Author: J. van Riet

# Libraries ----

library(dplyr)
library(future)
library(LesionSegR)
library(VariantAnnotation)

# Parallel settings.
future::plan(future::multisession, workers = 5)

# Import metadata. ----

metadata <- readr::read_tsv("~/test/test.tsv", show_col_types = FALSE) %>%
    dplyr::mutate(
        sample_strain = paste(sample_name, strain1, strain2, sep = '_'),
        seqname_strain = paste(sequencing_name, strain1, strain2, sep = '_'),
    )

workflow_dir <- "/omics/groups/OE0538/internal/users/e480l/projects/DEN_tumors/snakemake/DNA/Nextseq/38571/003_SM_output/"

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

## Determine WGS characteristics. ----

data_combined <- import_samples(metadata, workflow_dir, gtf = gtf)

## dN/dS. ----

data_combined$dNdS <- LesionSegR::run_dnds(data_combined$somaticvariants, path_db = "/omics/groups/OE0538/internal/projects/sharedData/GRCm39/annotation/refCDS_ENSEMBLv110_GRCm39.rda")
data_combined$dNdS <- data_results$dNdS$finalOutput

# Save the data. ----

saveRDS(data_combined, "~/odomLab/LesionSegregration_F1/data/rdata/data_combined.rds")