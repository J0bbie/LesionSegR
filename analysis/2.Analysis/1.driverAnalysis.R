library(LesionSegR)
library(dplyr)

# Load metadata. ----

metadata <- readr::read_tsv("/omics/groups/OE0538/internal/users/e480l/projects/DEN_tumors/snakemake/all_Novaseq_samples/all_Novaseq_samples_BlxCast.tsv", show_col_types = FALSE) %>%
    dplyr::mutate(
        sample_strain = paste(sample_name, strain1, strain2, sep = '_'),
        seqname_strain = paste(sequencing_name, strain1, strain2, sep = '_'),
    )

# Import data. ----
data_combined <- base::readRDS("/omics/groups/OE0538/internal/users/e480l/projects/DEN_tumors/snakemake/all_Novaseq_samples/all_data_LesionSegR/data_combined.rds")

mean(data_combined$tumorburden$totalMutations)
data_combined$tumorburden$totalMutations
# Generate overview of known driver genes. ----
known_drivers <- subset(data_combined$somaticvariants, SYMBOL %in% c("Kras", "Hras", "Egfr", "Braf") & ConsequenceAll == "missense_variant")
known_drivers <- tibble::as_tibble(known_drivers) %>% dplyr::inner_join(metadata, by = c('sample' = 'sequencing_name'))

table(known_drivers$sample, known_drivers$SYMBOL)
table(known_drivers$SYMBOL, known_drivers$matched_group)

sigminer::show_catalogue(data_combined$mutmatrices_sbs, mode = 'SBS', style = 'cosmic', samples = colnames(data_combined$mutmatrices_sbs)[sample(80, 10)])
sigminer::show_catalogue(data_combined$mutmatrices_sbs, 
                         mode = 'SBS', 
                         style = 'cosmic', 
                         x_lab = "Sequence context")  

x <- known_drivers %>% 
    dplyr::group_by(sample) %>% 
    dplyr::summarise(
        mut_type = paste(sort(unique(SYMBOL)), collapse = ', ')
    )

x <- x %>% dplyr::inner_join(metadata %>% dplyr::select(sequencing_name, sample_name), by = c('sample' = 'sequencing_name'))
dds_pca <- dds_pca %>% dplyr::left_join(x, by = c('sample_name' = 'sample_name'))

dds_pca %>% 
   # dplyr::filter(!grepl(', ', mut_type)) %>% 
    dplyr::filter(tissue == 'Tumor') %>% 
    ggplot2::ggplot(., ggplot2::aes(x = PC1, y = PC2, fill = mut_type)) +
    ggplot2::geom_point(shape = 21, color = 'black', size = 5) +
    ggplot2::scale_fill_brewer(palette = 'Dark2') +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size=20))+
    ggplot2::theme(legend.position="bottom") 
