#' @title Import and read haplotag log.
#' @description Determine the no. of tagged reads and used variants during haplotagging.
#' @param x (character): Paths to one or multiple flagstats file.
#' @return (tibble): Tibble with columns variable, value, and sample.
#' @importFrom dplyr %>%
#' @export
read_haplotag_log <- function(x){
    # Check input.
    checkmate::assertFile(x)
    
    futile.logger::flog.info(glue::glue("Importing {length(x)} haplotag log files."))
    
    data_haplotag <- tibble::as_tibble(dplyr::bind_rows(future.apply::future_lapply(x, function(sample_path) {
        data <- readr::read_tsv(sample_path, col_names = 'row', col_types = 'c') %>% 
            dplyr::filter(grepl('Found.*reads', row)) %>% 
            dplyr::mutate(
                total_tagged_reads = as.integer(base::gsub(' reads.*', '', base::gsub('Found ', '', row))),
                total_tagged_variants = as.integer(base::gsub(' variants', '',base::gsub('.*covering ', '', row))),
                sample = factor(base::gsub('\\.log', '', base::gsub('haplotag_WGS_', '', base::basename(sample_path))))
            )
        
        data$chrom <- paste0('chr', c(1:19, 'X', 'Y', 'M'))
        
        data <- data %>% 
            dplyr::distinct(chrom, total_tagged_reads, total_tagged_variants, sample)
        
        return(data)
        
    }))) %>% 
        dplyr::mutate(chrom = factor(chrom, levels = paste0('chr', c(1:19, 'X', 'Y', 'M'))))
    
    return(data_haplotag)
}