library(dplyr)
library(future)
library(LesionSegR)
library(VariantAnnotation)

source("analysis/themes.R")

# Import metadata of Ultima/Illumina samples. ----

metadata <- readxl::read_xlsx(
    path = "~/odomLab/LesionSegregration_F1/manuscript/tables/SupplTable1_OverviewSequencing.xlsx",
    sheet = "WGS (NovaSeq)", trim_ws = TRUE
) %>%
    # Only keep the malignant samples.
    dplyr::filter(!is.na(matched_group)) %>% 
    dplyr::filter(sample %in% c('AS-949292', 'AS-949296', 'AS-949300'))

# Flagstat overview. ----

data_ill <-  LesionSegR::read_flagstats(list.files('~/test/', pattern = 'flagstat', full.names = T)) %>% 
    dplyr::mutate(sample = gsub('_.*', '', sample))
data_ug <- LesionSegR::read_flagstats(list.files('~/odomLab/LesionSegregration_F1/UltimaGenomics/alignment/', pattern = 'flagstat', full.names = T)) %>% 
    dplyr::mutate(sample = gsub('-.*', '', sample))

haplo_ill <- LesionSegR::read_haplotag_log(list.files('~/odomLab/LesionSegregration_F1/data/workflow/logs/haplotyping/', pattern = 'log', full.names = T)) %>% 
    dplyr::mutate(sample = gsub('_.*', '', sample)) %>% 
    dplyr::filter(chrom == 'chr10', sample %in% c('AS-949292', 'AS-949296', 'AS-949300'))

haplo_ug <- LesionSegR::read_haplotag_log(list.files('~/odomLab/LesionSegregration_F1/UltimaGenomics/logs/haplotyping/', pattern = 'log', full.names = T)) %>%
    dplyr::filter(chrom == 'chr10') %>% 
    dplyr::mutate(sample = gsub('haplotag_', '', gsub('-.*', '', sample)))

dplyr::bind_rows(data_ill, data_ug) %>% 
    dplyr::left_join(dplyr::bind_rows(haplo_ug, haplo_ill)) %>% 
    dplyr::mutate(
        sample = ifelse(grepl('AS-', sample), sprintf("(IL) %s", sample), sprintf("(UG) %s", sample)),
        sample = ifelse(grepl('292', sample), "(IL) CABL1t1", sample),
        sample = ifelse(grepl('296', sample), "(IL) CABL2t1", sample),
        sample = ifelse(grepl('300', sample), "(IL) CABL3t1", sample),
) %>% 
    dplyr::mutate(
        "Strain" = "B6xCAST_EiJ",
        mapped_reads = formattable::percent(Mapped / `Total reads`,2),
        paired_reads = formattable::percent(`Paired in sequencing` / `Total reads`,2),
        total_reads = formattable::color_tile("lightpink", max.color = "hotpink")(base::formatC(`Total reads`, drop0trailing = TRUE, big.mark = ",", format = 'd')),
        Duplicates = formattable::percent(Duplicates / `Total reads`,2),
        total_tagged_reads = formattable::color_tile("lightpink", max.color = "skyblue")(formattable::percent(total_tagged_reads / `Total reads`,2)),
    ) %>%
    dplyr::select(
        Sample = sample,
        `Strain`,
        "Total reads" = total_reads,
        "Mapped" = mapped_reads,
        "Paired" = paired_reads,
        Duplicates,
        "Haplo-reads" = total_tagged_reads,
        "Σ(SNPs)" = total_tagged_variants,
    ) %>%
    knitr::kable(escape = FALSE, align = "llccccccc", format.args = list(decimal.mark = ".", big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15, full_width = TRUE, html_font = "Roboto", bootstrap_options = c("striped")) %>%
    kableExtra::add_header_above(line_sep = 3, c(" " = 3, "Flagstats" = 3, "Haplotag" = 2))


# Import somatic variants. ----

## Import somatic variants (Illumina). ----
somatics_illumina <- base::readRDS("~/odomLab/LesionSegregration_F1/data/rdata/data_somaticvariants.rds")

# Subset on chr 10.
somatics_illumina <- VariantAnnotation::VRangesList(lapply(somatics_illumina, function(x){
    x[GenomicRanges::seqnames(x) == 'chr10',]
}))


## Import somatic variants (Ultima). ----
files_vcf <- base::list.files(path = "~/test/", pattern = ".*haplocounted.removedblacklist.vcf.gz$", full.names = TRUE)
somatics_ultima <- LesionSegR::import_vcf_mutect(files_vcf)

# Generate Venn-diagram. ----

## Determine overlapping variants.
IRanges::subsetByOverlaps(somatics_ultima$CABL1t1, paul292, invert = F)
IRanges::subsetByOverlaps(somatics_ultima$CABL2t1, paul292, invert = F)
IRanges::subsetByOverlaps(somatics_ultima$CABL3t1, paul292, invert = F)

## Function.
generateVenn <- function(x, y){
    x_snv <- x[x$mutType == 'SNV']
    y_snv <- y[y$mutType == 'SNV']
    set1 = paste0(seqnames(x_snv), start(x_snv), end(x_snv), ref(x_snv), alt(x_snv))
    set2 = paste0(seqnames(y_snv), start(y_snv), end(y_snv), ref(y_snv), alt(y_snv))
    
    x_id <- x[x$mutType != 'SNV']
    y_id <- y[y$mutType != 'SNV']
    set3 = paste0(seqnames(x_id), start(x_id), end(x_id), ref(x_id), alt(x_id))
    set4 = paste0(seqnames(y_id), start(y_id), end(y_id), ref(y_id), alt(y_id))
    
    VennDiagram::venn.diagram(
        
        x = list(set1, set2, set3, set4), 
        filename = NULL,
        category.names = c('(UG) SNV', '(IL) SNV', '(UG) InDel', '(IL) InDel'),
        output = TRUE ,
        lwd = .5,
        imagetype="png" ,
        height = 350, 
        width = 350,
        resolution = 200,
        disable.logging = T,
        # Colors.
        col=c("#CA534F", '#54A585', '#4FAEDA', '#FCD494'),
        fill = c(ggplot2::alpha("#CA534F",0.1), ggplot2::alpha('#54A585',0.1), ggplot2::alpha('#4FAEDA',0.1), ggplot2::alpha('#FCD494',0.1)),
        
        # Numbers
        cex = .6, 
        fontface = "bold",
        fontfamily = "Roboto",
        sub.fontfamily = "Roboto",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer"
    )
    
}

## Combine Venns. ----
a = generateVenn(x = somatics_ultima$CABL1t1, y = somatics_illumina$`AS-949292_B6_CAST_EiJ`)
b = generateVenn(x = somatics_ultima$CABL2t1, y = somatics_illumina$`AS-949296_B6_CAST_EiJ`)
c = generateVenn(x = somatics_ultima$CABL3t1, y = somatics_illumina$`AS-949300_B6_CAST_EiJ`)

cowplot::plot_grid(grid::grobTree(a), grid::grobTree(b), grid::grobTree(c), label_size = 10, label_fontfamily = 'Roboto', nrow = 1)


# Determine no. of haplotyped somatic variants.

data_distributions_ill <- dplyr::bind_rows(
    base::lapply(somatics_illumina, function(x){
        tibble::as_tibble(data.frame(sample = Biobase::sampleNames(x), chrom = seqnames(x), S4Vectors::mcols(x)[c("H1_Alt", "H2_Alt", "origin_mutant")]))
    } )
) %>%
    dplyr::mutate(
        sample = as.character(sample),
        sample = ifelse(grepl(292, sample), 'CABL1t1', sample),
        sample = ifelse(grepl(296, sample), 'CABL2t1', sample),
        sample = ifelse(grepl(300, sample), 'CABL3t1', sample),
        delta = H1_Alt - H2_Alt,
        origin_mutant2 = ifelse(delta < 0, 'CAST', 'B6'),
        source = 'Illumina',
        origin_mutant2 = ifelse(origin_mutant == 'UA', 'UA', origin_mutant2)
    ) %>%
    dplyr::filter(delta != 0) %>% 
    dplyr::filter(grepl('chr10', chrom)) %>% 
    dplyr::filter(grepl('CABL', sample))


data_distributions_ug <- dplyr::bind_rows(
    base::lapply(somatics_ultima, function(x){
        tibble::as_tibble(data.frame(sample = Biobase::sampleNames(x), chrom = seqnames(x), S4Vectors::mcols(x)[c("H1_Alt", "H2_Alt", "origin_mutant")]))
    } )
) %>%
    dplyr::mutate(
        delta = H1_Alt - H2_Alt,
        origin_mutant2 = ifelse(delta < 0, 'CAST', 'B6'),
        origin_mutant2 = ifelse(origin_mutant == 'UA', 'UA', origin_mutant2),
        
        source = 'Ultima Genomics'
    ) %>%
    dplyr::filter(delta != 0) %>% 
    dplyr::filter(grepl('chr10', chrom))

dplyr::bind_rows(data_distributions_ug, data_distributions_ill) %>% 
    # Plot distribution.
    ggplot2::ggplot(., ggplot2::aes(x = delta, fill = origin_mutant2)) +
    ggplot2::geom_histogram(binwidth = .5, na.rm = TRUE, color = "grey10", lwd = ggplot2::rel(.33)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) +
    ggplot2::scale_x_continuous(limits = c(-25, 25), breaks = seq(-25, 25, 5)) +
    ggplot2::scale_fill_manual(values = color_scheme, name = NULL) +
    ggplot2::labs(
        x = "H1 - H2 assigned reads",
        y = "Frequency<br><sup>(No. of somatic variants)</sup>"
    ) +
    ggplot2::facet_wrap(source~sample, ncol = 3, as.table = T) +
    theme_job

## Mutational signature analysis. ----

mutmatrix_ultima <- LesionSegR::generate_mutmatrices_96(somatics_ultima)
mutmatrix_illumina <- LesionSegR::generate_mutmatrices_96(somatics_illumina)

subset_samples <- function(x){
    colnames(x) <- gsub('_B6_CAST_EiJ', '', colnames(x))
    x <- x[,colnames(x) %in% c('AS-949292', 'AS-949296', 'AS-949300')]
    
    colnames(x) <- c('(IL) CABL1t1', '(IL) CABL2t1', '(IL) CABL3t1')
    return(x)
}

sigminer::show_catalogue(
    subset_samples(mutmatrix_illumina$sbs), 
    mode = "SBS", 
    style = "cosmic",
    samples = c('(IL) CABL1t1', 'CABL1t1', '(IL) CABL2t1', 'CABL2t1', '(IL) CABL3t1', 'CABL3t1'), method = 'W'
)
