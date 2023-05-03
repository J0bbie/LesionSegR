# Function: Make a prelim. overview of possible sequencing samples.

# Load libraries. ----
library(dplyr)


# Retrieve project samples and clean-up descriptions. ----
sequencingRuns <- list()

## B6 / CAS ----

sequencingRuns$f1_b6_mcas <- tibble::tibble(f = list.files(sprintf('/omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/sequencingRuns/f1_b6_mcas/sequencing/%s/view-by-pid/', c('whole_genome_sequencing', 'rna_sequencing')), recursive = T, all.files = F, pattern = 'fastq.gz$', full.names = T)) %>% 
    dplyr::mutate(
        ASID = gsub('-LR.*', '', basename(f)),
        sequencingType = ifelse(grepl('whole_genome_sequencing', f), 'WGS', 'Other'),
        sequencingType = ifelse(grepl('rna_sequencing', f), 'WTS', sequencingType),
        description = stringr::str_split(f, pattern = '/',simplify = T)[, 16],
        sequencingRun = 'f1_b6_mcas'
    ) %>% 
    dplyr::distinct(ASID, sequencingType, description, sequencingRun)

## C3H / CAS ----
sequencingRuns$f1_mcas_c3h <- tibble::tibble(f = list.files(sprintf('/omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/sequencingRuns/f1_mcas_c3h/sequencing/%s/view-by-pid/', c('whole_genome_sequencing', 'rna_sequencing')), recursive = T, all.files = F, pattern = 'fastq.gz$', full.names = T)) %>% 
    dplyr::mutate(
        ASID = gsub('-LR.*', '', basename(f)),
        sequencingType = ifelse(grepl('whole_genome_sequencing', f), 'WGS', 'Other'),
        sequencingType = ifelse(grepl('rna_sequencing', f), 'WTS', sequencingType),
        description = stringr::str_split(f, pattern = '/',simplify = T)[, 16],
        sequencingRun = 'f1_mcas_c3h'
    ) %>% 
    dplyr::distinct(ASID, sequencingType, description, sequencingRun)

## Normal mmus ----
sequencingRuns$mmus <- tibble::tibble(f = list.files('/omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/sequencingRuns/mmus/sequencing/whole_genome_sequencing/view-by-pid/', recursive = T, all.files = F, pattern = 'fastq.gz$', full.names = T)) %>% 
    dplyr::mutate(
        ASID = gsub('-LR.*', '', basename(f)),
        sequencingType = 'WGS',
        description = stringr::str_split(f, pattern = '/',simplify = T)[, 16],
        sequencingRun = 'mmus',
    ) %>% 
    dplyr::distinct(ASID, sequencingType, description, sequencingRun)

## Reciprocal crosses. ----
sequencingRuns$run230224 <- readr::read_tsv('/omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/sequencingRuns/230224/33210_meta.tsv') %>% 
    dplyr::filter(!is.na(SEQUENCING_TYPE)) %>% 
    dplyr::mutate(
        ASID = gsub('-LR.*', '', basename(FASTQ_FILE)),
        sequencingType = 'WGS',
        sequencingPlatform = 'Illumina NovaSeq 6000',
        description = SAMPLE_NAME,
        sequencingRun = '230224',
    ) %>% 
    dplyr::distinct(ASID, sequencingType, sequencingPlatform, description, sequencingRun)


## WTS. ----
sequencingRuns$run230222 <- readr::read_tsv('/omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/sequencingRuns/230222_VH00211_234_AAC77JYM5/230222_VH00211_234_AAC77JYM5_meta.tsv') %>% 
    dplyr::filter(!is.na(SEQUENCING_TYPE)) %>% 
    dplyr::mutate(
        ASID = gsub('-LR.*', '', basename(FASTQ_FILE)),
        sequencingType = 'WTS',
        sequencingPlatform = 'Illumina NextSeq 2000',
        description = SAMPLE_NAME,
        sequencingRun = '230222',
    ) %>% 
    dplyr::distinct(ASID, sequencingType, sequencingPlatform, description, sequencingRun)

# Output table.
write.table(x = dplyr::bind_rows(sequencingRuns), row.names = F, append = F, quote = F, sep = '\t', file = 'asd.txt')
