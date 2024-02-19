
# Load metadata. ----

metadata <- readr::read_tsv("/omics/groups/OE0538/internal/users/e480l/projects/DEN_tumors/snakemake/all_Novaseq_samples/all_Novaseq_samples_BlxCast.tsv") %>% 
    dplyr::mutate(
        tissue = factor(dplyr::case_when(
            grepl('Normal', group) ~ 'Normal',
            grepl('Near', group) ~ 'Near_adjacent', 
            .default = 'Tumor'
        )),
        strain = factor(dplyr::if_else(grepl('CAST/B6', group), 'CAST_B6', 'B6_CAST')),
        tissue_strain = factor(paste(tissue, strain, sep = '_')),
        mice_id = gsub('CasBl|BlCas', '', sample_name)
    )

# Subset GTF ----

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


# Import counts. ----

files_featurecounts <- list.files(
    "/omics/odcf/analysis/OE0538_projects/DO-0006/f1_b6_mcas/e480l/projects/DEN_tumors/RNA/Samples_1-50_52-87-89-90_92/Novaseq/38415/SM_output_new_VCF/38415_Novaseq_RNA/counting/WTS/", 
    pattern = '_counts.txt$', 
    full.names = T
)

data_counts <- dplyr::bind_rows(pbapply::pblapply(files_featurecounts, function(x){
    data <- data.table::fread(x, col.names = c('gene_id', 'chr', 'start', 'end', 'strand', 'length', 'totalCounts', 'H1', 'H2', 'UA')) %>% 
        dplyr::distinct(gene_id, totalCounts) %>% 
        dplyr::inner_join(gtf, by = 'gene_id') %>% 
        dplyr::mutate(
            gene_name = factor(gene_name),
            sample = factor(gsub('_.*', '', basename(x)))
        ) %>% 
        dplyr::select(gene_name, totalCounts, sample)

    return(data)
}, cl = 10)) %>% 
    tidyr::pivot_wider(., names_from = c(sample), values_from = totalCounts)

# Generate DESeq2 object. ----

m <- as.matrix(data_counts[2:ncol(data_counts)])
rownames(m) <- data_counts$gene_name

colData <- metadata %>% dplyr::filter(sequencing_name %in% colnames(m)) %>%  dplyr::arrange(ordered(sequencing_name, colnames(m)))
dds <- DESeq2::DESeqDataSetFromMatrix(m, colData = colData, design = ~tissue)
dds <- DESeq2::DESeq(dds, test = 'Wald', parallel = TRUE)

# Perform PCA on vst-counts. ----
dds_counts <- DESeq2::vst(dds, blind = T) %>% 
    SummarizedExperiment::assay(.)

# Select discriminating genes. ----

genes_de <- DESeq2::results(dds, contrast = c('tissue', 'Tumor', 'Normal'), tidy = T) %>% dplyr::filter(padj < 0.01, abs(log2FoldChange) >= .5, lfcSE < 1)
dds_counts_de <- dds_counts[rownames(dds_counts) %in% genes_de$row,]
dds_counts_var <- dds_counts[rownames(dds_counts) %in% names(sort(-rowVars(dds_counts))[1:1000]),]

dds_pca <- broom::tidy(prcomp(t(dds_counts_de))) %>% 
    tidyr::pivot_wider(names_from = PC, names_prefix = 'PC') %>% 
    dplyr::mutate(sequencing_name = row) %>% 
    dplyr::inner_join(metadata)

# Plot ----

ggplot2::ggplot(dds_pca, ggplot2::aes(x = PC1, y = PC2, fill = strain, label = tissue)) +
    ggplot2::geom_point(shape = 21, color = 'black') +
    ggforce::geom_mark_ellipse(mapping = aes(fill = tissue), expand = 0.01, alpha = 0.1) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, label.position = 'top', nrow = 1)) +
    ggplot2::scale_fill_manual(values = c('Normal' = 'skyblue', 'Near_adjacent' = 'darkblue', 'Tumor' = 'orange', 'B6_CAST' = 'hotpink', 'CAST_B6' = 'darkred')) +
    ggplot2::scale_y_continuous(limits = c(-50, 50)) +
    ggplot2::scale_x_continuous(limits = c(-100, 100)) +
    scir::theme_ggplot()
