#' @title Import and read flagstats
#' @description Import and read flagstats
#' @param x (character): Paths to one or multiple flagstats file.
#' @return (tibble): Tibble with columns variable, value, and sample.
#' @importFrom dplyr %>%
#' @export
read_flagstats <- function(x) {
    # Check input.
    checkmate::assertFile(x)
    
    futile.logger::flog.info(glue::glue("Importing {length(x)} flagstat files."))
    
    # Read flagstats.
    data_flagstats <- tibble::as_tibble(dplyr::bind_rows(future.apply::future_lapply(x, function(sample_path) {
        # Import and clean flagstats.
        data <- readr::read_csv(sample_path, col_names = FALSE, show_col_types = FALSE) %>%
            dplyr::mutate(
                value = as.integer(base::trimws(gsub("\\+.*", "", X1))),
                sample = basename(sample_path) %>% stringr::str_replace_all(pattern = "_sortedByCoord_markDup_haplotagged.bam.flagstat", replacement = "")
            ) %>%
            dplyr::select(value, sample)
        
        # Import and clean no. of primary-aligned reads per chromosome.
        if(file.exists(paste0(sample_path, '_primaries_per_chr.txt'))){
            data_count <- readr::read_delim(paste0(sample_path, '_primaries_per_chr.txt'), delim = ' ', col_names =c('count', 'chrom'), show_col_types = FALSE) %>%
                dplyr::filter(!grepl('chrY|chrM', chrom)) %>% 
                dplyr::summarise(
                    value = sum(count),
                    variable = "Total haplotaggable reads",
                    sample = basename(sample_path) %>% stringr::str_replace_all(pattern = "_sortedByCoord_markDup_haplotagged.bam.flagstat", replacement = "")
                )
        }else{
            data_count <- tibble::tibble(
                value = NA, 
                variable = "Total haplotaggable reads",
                sample = basename(sample_path) %>% stringr::str_replace_all(pattern = "_sortedByCoord_markDup_haplotagged.bam.flagstat", replacement = "")
            )
        }
        
        data$variable <- c(
            "Total reads",
            "Secondary reads",
            "Supplementary reads",
            "Duplicates",
            "Mapped",
            "Paired in sequencing",
            "Read 1",
            "Read 2",
            "Properly paired",
            "With itself and mate mapped",
            "Singletons",
            "Mate mapped to diff. chr",
            "Mate mapped to diff. chr (mapQ>=5")
        
        data <- dplyr::bind_rows(data, data_count)
        
        return(data)
    })))
    
    # Spread wide.
    data_flagstats <- data_flagstats %>%
        tidyr::spread(key = variable, value = value)
    
    # Return.
    return(data_flagstats)
}
