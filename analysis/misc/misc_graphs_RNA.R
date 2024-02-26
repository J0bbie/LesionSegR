library(magrittr)
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

dds_counts_var <- dds_counts[rownames(dds_counts) %in% names(sort(-matrixStats::rowVars(dds_counts))[1:10]),]

gene_list <- c(names(sort(-matrixStats::rowVars(dds_counts))[1:20]), c("Egfr", "Kras", "Hras", "Braf"))
dds_counts_selected <- dds_counts[rownames(dds_counts) %in% gene_list,]

dds_pca <- broom::tidy(prcomp(t(dds_counts_var))) %>% 
    tidyr::pivot_wider(names_from = PC, names_prefix = 'PC') %>% 
    dplyr::mutate(sequencing_name = row) %>% 
    dplyr::inner_join(metadata)

# Plot ----

ggplot2::ggplot(dds_pca, ggplot2::aes(x = PC1, y = PC2, fill = strain, label = tissue)) +
    ggplot2::geom_point(shape = 21, color = 'black') +
    ggforce::geom_mark_ellipse(mapping = ggplot2::aes(fill = tissue), expand = ggplot2::unit(1, 'mm'), alpha = 0.1) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, label.position = 'top', nrow = 1)) +
    ggplot2::scale_fill_manual(values = c('Normal' = 'skyblue', 'Near_adjacent' = 'darkblue', 'Tumor' = 'orange', 'B6_CAST' = 'hotpink', 'CAST_B6' = 'darkred')) +
    ggplot2::scale_y_continuous(limits = c(-8, 8)) +
    ggplot2::scale_x_continuous(limits = c(-20, 20)) +
    ggplot2::theme(legend.position="bottom") #+
    #scir::theme_ggplot()

DESeq2::plotCounts(dds, gene = 'H19', intgroup = 'tissue', normalized = T)






# no magic anymore ------------------------

# PCA on most variable genes --------------
pca <- stats::prcomp(t(dds_counts_de))
pca_data <- tibble::as_tibble(pca$x, rownames = 'sample') %>%
    #add metadata 
    dplyr::inner_join(metadata, by = c(sample = "sequencing_name"))


# Plot PCA ----
ggplot2::ggplot(dds_pca, ggplot2::aes(x = PC1, y = PC2, fill = strain, label = tissue)) +
    ggplot2::geom_point(shape = 21, color = 'black') +
    ggforce::geom_mark_ellipse(mapping = ggplot2::aes(fill = tissue), expand = 0.01, alpha = 0.1) +
    ggplot2::scale_fill_manual(values = c('Normal' = 'skyblue', 'Near_adjacent' = 'darkblue', 'Tumor' = 'orange', 'B6_CAST' = 'hotpink', 'CAST_B6' = 'darkred')) 
   

# heatmap -----------------
dds_counts_heatmap <-  tibble::as_tibble(dds_counts, rownames = NA) %>% 
    tibble::rownames_to_column(var = "gene") %>%
    tidyr::pivot_longer(cols = colnames(.)[-1]) %>% 
    dplyr::inner_join(metadata, by = c(name = "sequencing_name"))

# add log expression of genes 
dds_counts_heatmap$log.expression <- log(dds_counts_heatmap$value)

# DEG
dds_counts_heatmap_de <-  tibble::as_tibble(dds_counts_de, rownames = NA) %>% 
    tibble::rownames_to_column(var = "gene") %>%
    tidyr::pivot_longer(cols = colnames(.)[-1]) %>% 
    dplyr::inner_join(metadata, by = c(name = "sequencing_name"))

# add log expression of genes 
dds_counts_heatmap_de$log.expression <- log(dds_counts_heatmap_de$value)


# mots variable genes
dds_counts_heatmap_var <-  tibble::as_tibble(dds_counts_var, rownames = NA) %>% 
    tibble::rownames_to_column(var = "gene") %>%
    tidyr::pivot_longer(cols = colnames(.)[-1]) %>% 
    dplyr::inner_join(metadata, by = c(name = "sequencing_name"))

# add log expression of genes 
dds_counts_heatmap_var$log.expression <- log(dds_counts_heatmap_var$value)

#heatmap 

# exp.heatmap <- ggplot2::ggplot(data = dds_counts_heatmap, mapping = ggplot2::aes(x = sample_name,
#                                                      y = gene,
#                                                      fill = log.expression)) +
#     ggplot2::geom_tile() +
#     ggplot2::xlab(label = "Sample") + # Add a nicer x-axis title
#     ggplot2::theme(axis.title.y = ggplot2::element_blank(), # Remove the y-axis title
#           axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5)) # Rotate the x-axis labels
# 
# exp.heatmap



exp.heatmap_per_tissue <- ggplot2::ggplot(data = dds_counts_heatmap_var, mapping = ggplot2::aes(x = sample_name,
                                                     y = gene,
                                                     fill = log.expression)) +
    ggplot2::geom_tile() +
    ggplot2::xlab(label = "Sample") +
    # facet_grid makes two panels, one for control, one for flu:
    ggplot2::facet_grid(~ tissue, switch = "x", scales = "free_x", space = "free_x") + 
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.5))

exp.heatmap_per_tissue


# H19 graphs --------------
H19_counts <- setNames(tibble::as_tibble(dds_counts[rownames(dds_counts) == "H19",], rownames = NA), "H19") %>% 
    tibble::rownames_to_column(var = "sequencing_name") %>%
    tidyr::pivot_longer(cols = H19) %>% 
    dplyr::inner_join(metadata, by = "sequencing_name")

H19_counts %>% ggplot2::ggplot(ggplot2::aes(x= tissue, y = value, fill = strain, color = strain, label = sample_name)) + 
    ggplot2::geom_jitter(shape = 21, color = 'black') 
    


ggplot2::ggplot(dds_pca, ggplot2::aes(x = PC1, y = PC2, fill = strain, label = tissue)) +
    ggplot2::geom_point(shape = 21, color = 'black') +
    ggforce::geom_mark_ellipse(mapping = ggplot2::aes(fill = tissue), expand = ggplot2::unit(1, 'mm'), alpha = 0.1) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, label.position = 'top', nrow = 1)) +
    ggplot2::scale_fill_manual(values = c('Normal' = 'skyblue', 'Near_adjacent' = 'darkblue', 'Tumor' = 'orange', 'B6_CAST' = 'hotpink', 'CAST_B6' = 'darkred')) +
    ggplot2::scale_y_continuous(limits = c(-8, 8)) +
    ggplot2::scale_x_continuous(limits = c(-20, 20)) +
    ggplot2::theme(legend.position="bottom") 






### trying to extract H1 and H2 counts 
data_counts_H1 <- dplyr::bind_rows(pbapply::pblapply(files_featurecounts, function(x){
    data <- data.table::fread(x, col.names = c('gene_id', 'chr', 'start', 'end', 'strand', 'length', 'totalCounts', 'H1', 'H2', 'UA')) %>% 
        dplyr::distinct(gene_id, H1) %>% 
        dplyr::inner_join(gtf, by = 'gene_id') %>% 
        dplyr::mutate(
            gene_name = factor(gene_name),
            sample = factor(gsub('_.*', '', basename(x)))
        ) %>% 
        dplyr::select(gene_name, H1, sample)
    
    return(data)
}, cl = 10)) %>% 
    tidyr::pivot_wider(., names_from = c(sample), values_from = H1)


data_counts_H2 <- dplyr::bind_rows(pbapply::pblapply(files_featurecounts, function(x){
    data <- data.table::fread(x, col.names = c('gene_id', 'chr', 'start', 'end', 'strand', 'length', 'totalCounts', 'H1', 'H2', 'UA')) %>% 
        dplyr::distinct(gene_id, H2) %>% 
        dplyr::inner_join(gtf, by = 'gene_id') %>% 
        dplyr::mutate(
            gene_name = factor(gene_name),
            sample = factor(gsub('_.*', '', basename(x)))
        ) %>% 
        dplyr::select(gene_name, H2, sample)
    
    return(data)
}, cl = 10)) %>% 
    tidyr::pivot_wider(., names_from = c(sample), values_from = H2)






# Generate DESeq2 object. ----

m <- as.matrix(data_counts_H1[2:ncol(data_counts_H1)])
rownames(m) <- data_counts_H1$gene_name

colData <- metadata %>% dplyr::filter(sequencing_name %in% colnames(m)) %>%  dplyr::arrange(ordered(sequencing_name, colnames(m)))
dds <- DESeq2::DESeqDataSetFromMatrix(m, colData = colData, design = ~tissue)
dds <- DESeq2::DESeq(dds, test = 'Wald', parallel = TRUE)

# Perform PCA on vst-counts. ----
dds_counts_H1 <- DESeq2::vst(dds, blind = T) %>% 
    SummarizedExperiment::assay(.)


m <- as.matrix(data_counts_H2[2:ncol(data_counts_H2)])
rownames(m) <- data_counts_H2$gene_name

colData <- metadata %>% dplyr::filter(sequencing_name %in% colnames(m)) %>%  dplyr::arrange(ordered(sequencing_name, colnames(m)))
dds <- DESeq2::DESeqDataSetFromMatrix(m, colData = colData, design = ~tissue)
dds <- DESeq2::DESeq(dds, test = 'Wald', parallel = TRUE)

# Perform PCA on vst-counts. ----
dds_counts_H2 <- DESeq2::vst(dds, blind = T) %>% 
    SummarizedExperiment::assay(.)



H19_counts_H1 <- setNames(tibble::as_tibble(dds_counts_H1[rownames(dds_counts_H1) == "H19",], rownames = NA), "H19") %>% 
    tibble::rownames_to_column(var = "sequencing_name") %>%
    tidyr::pivot_longer(cols = H19) %>% 
    dplyr::inner_join(metadata, by = "sequencing_name")

H19_counts_H2 <- setNames(tibble::as_tibble(dds_counts_H2[rownames(dds_counts_H2) == "H19",], rownames = NA), "H19") %>% 
    tibble::rownames_to_column(var = "sequencing_name") %>%
    tidyr::pivot_longer(cols = H19) %>% 
    dplyr::inner_join(metadata, by = "sequencing_name")

H19_counts_H1 %>% ggplot2::ggplot(ggplot2::aes(x= tissue, y = value, fill = strain, color = strain, label = sample_name)) + 
    ggplot2::geom_jitter(shape = 21, color = 'black') +
    ggplot2::ggtitle("H19, Bl6 haplotype") +
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_fill_manual(values = c('B6_CAST' = 'hotpink', 'CAST_B6' = 'skyblue')) +
    ggplot2::theme(legend.position="bottom") 


H19_counts_H2 %>% ggplot2::ggplot(ggplot2::aes(x= tissue, y = value, fill = strain, color = strain, label = sample_name)) + 
    ggplot2::geom_jitter(shape = 21, color = 'black') +
    ggplot2::ggtitle("H19, CAST_EiJ haplotype") +
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_fill_manual(values = c('B6_CAST' = 'hotpink', 'CAST_B6' = 'skyblue')) +
    ggplot2::theme(legend.position="bottom") 


TLR7_counts_H1 <- setNames(tibble::as_tibble(dds_counts_H1[rownames(dds_counts_H1) == "Tlr7",], rownames = NA), "Tlr7") %>% 
    tibble::rownames_to_column(var = "sequencing_name") %>%
    tidyr::pivot_longer(cols = Tlr7) %>% 
    dplyr::inner_join(metadata, by = "sequencing_name")

TLR7_counts_H1 %>% ggplot2::ggplot(ggplot2::aes(x= tissue, y = value, fill = strain, color = strain, label = sample_name)) + 
    ggplot2::geom_jitter(shape = 21, color = 'black') +
    ggplot2::ggtitle(label = "TLR7, Bl6 haplotype") +
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_fill_manual(values = c('B6_CAST' = 'hotpink', 'CAST_B6' = 'skyblue')) +
    ggplot2::theme(legend.position="bottom") 



TLR7_counts_H2 <- setNames(tibble::as_tibble(dds_counts_H2[rownames(dds_counts_H2) == "Tlr7",], rownames = NA), "Tlr7") %>% 
    tibble::rownames_to_column(var = "sequencing_name") %>%
    tidyr::pivot_longer(cols = Tlr7) %>% 
    dplyr::inner_join(metadata, by = "sequencing_name")

TLR7_counts_H2 %>% ggplot2::ggplot(ggplot2::aes(x= tissue, y = value, fill = strain, color = strain, label = sample_name)) + 
    ggplot2::geom_jitter(shape = 21, color = 'black') +
    ggplot2::ggtitle(label = "TLR7, CAST_EiJ haplotype") +
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_fill_manual(values = c('B6_CAST' = 'hotpink', 'CAST_B6' = 'skyblue')) +
    ggplot2::theme(legend.position="bottom") 



# Gpc3

GPC3_counts <- setNames(tibble::as_tibble(dds_counts[rownames(dds_counts) == "Gpc3",], rownames = NA), "Gpc3") %>% 
    tibble::rownames_to_column(var = "sequencing_name") %>%
    tidyr::pivot_longer(cols = Gpc3) %>% 
    dplyr::inner_join(metadata, by = "sequencing_name")

GPC3_counts %>% ggplot2::ggplot(ggplot2::aes(x= tissue, y = value, fill = strain, color = strain, label = sample_name)) + 
    ggplot2::geom_jitter(shape = 21, color = 'black') +
    ggplot2::ggtitle(label = "GPC3, total counts") +
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_fill_manual(values = c('B6_CAST' = 'hotpink', 'CAST_B6' = 'skyblue')) +
    ggplot2::theme(legend.position="bottom") 


Cyp2f2_counts <- setNames(tibble::as_tibble(dds_counts[rownames(dds_counts_H2) == "Gpc3",], rownames = NA), "Gpc3") %>% 
    tibble::rownames_to_column(var = "sequencing_name") %>%
    tidyr::pivot_longer(cols = Gpc3) %>% 
    dplyr::inner_join(metadata, by = "sequencing_name")

Cyp2f2_counts %>% ggplot2::ggplot(ggplot2::aes(x= tissue, y = value, fill = strain, color = strain, label = sample_name)) + 
    ggplot2::geom_jitter(shape = 21, color = 'black') +
    ggplot2::ggtitle(label = "GPC3, total counts") +
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::scale_fill_manual(values = c('B6_CAST' = 'hotpink', 'CAST_B6' = 'skyblue')) +
    ggplot2::theme(legend.position="bottom") 

