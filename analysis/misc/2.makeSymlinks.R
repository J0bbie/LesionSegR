# Function: Generate symlinks of included samples (lane fastq.fz files) to one main folder.

# Load libraries. ----
library(dplyr)

# Create symlinks. ----

## WGS ----
x <- readxl::read_xlsx('/omics/groups/OE0538/internal/projects/LesionSegregration_F1/manuscript/tables/SupplTable1_OverviewSequencing.xlsx', sheet = 'WGS (NovaSeq)')
files <- tibble::tibble(f = list.files('/omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/sequencingRuns/', recursive = T, pattern = 'fastq.gz$', full.names = T)) %>% 
    dplyr::mutate(ASID = gsub('-LR.*', '', basename(f))) %>% 
    dplyr::inner_join(x, by = c('ASID' = 'sample')) %>% 
    dplyr::mutate(CMD = sprintf('ln -s %s /omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/rawreads/WGS/', f))

write.table(files$CMD, file = '~/test/symlinkWGS.txt', row.names = F, sep = '\t', quote = F, col.names = F)

## WTS ----

x <- readxl::read_xlsx('/omics/groups/OE0538/internal/projects/LesionSegregration_F1/manuscript/tables/SupplTable1_OverviewSequencing.xlsx', sheet = 'WTS (All)') %>% 
    dplyr::filter(sequencing_platform == "Illumina NextSeq 2000")

files <- tibble::tibble(f = list.files('/omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/sequencingRuns/', recursive = T, pattern = 'fastq.gz$', full.names = T)) %>% 
    dplyr::mutate(ASID = gsub('-LR.*', '', basename(f))) %>% 
    dplyr::inner_join(x, by = c('ASID' = 'sample')) %>% 
    dplyr::mutate(CMD = sprintf('ln -s %s /omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/rawreads/WTS/', f))

write.table(files$CMD, file = '~/test/symlinkWTS.txt', row.names = F, sep = '\t', quote = F, col.names = F)


## Ultima ----

# Convert CRAM to .fastq
x <- readxl::read_xlsx('/omics/groups/OE0538/internal/projects/LesionSegregration_F1/metadata/SupplTable1_OverviewSequencing.xlsx', sheet = 'Ultima') %>% dplyr::filter(sequencingType == 'WGS')
files <- tibble::tibble(f = list.files('/omics/gpcf/midterm/034400/data/UltimaGenomics/wgs/cram_files/', pattern = '*cabl.*cram$', full.names = T, ignore.case = T), recursive = T, pattern = 'cram$', full.names = T) %>% 
    dplyr::mutate(
        UltimaID = gsub('.*-', '', gsub('-UGA.*', '', basename(f))),
        runID = gsub('-.*', '', basename(f))
        ) %>% 
    dplyr::inner_join(x) %>% 
    dplyr::mutate(CMD = sprintf('samtools fastq --reference /omics/gpcf/midterm/034400/data/UltimaGenomics/Reference/GRCm38.primary_assembly.genome.fa %s > /omics/groups/OE0538/internal/projects/LesionSegregration_F1/data/UltimaGenomics/reads/untrimmed/%s_run%s.fastq', f, ASID, runID))

write.table(files$CMD, file = '~/test/UltimaWGS.txt', row.names = F, sep = '\t', quote = F, col.names = F)
