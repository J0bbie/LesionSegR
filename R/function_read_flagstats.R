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
                sample = gsub("_.*", "", basename(sample_path))
            ) %>%
            dplyr::select(value, sample)
        
        # Import and clean no. of primary-aligned reads per chromosome.
        if(file.exists(paste0(sample_path, '_primaries_per_chr.txt'))){
            data_count <- data.table::fread(paste0(sample_path, '_primaries_per_chr.txt'), strip.white = T, col.names = c('count', 'chrom')) %>%
                dplyr::filter(!grepl('chrY|chrM', chrom)) %>% 
                dplyr::summarise(
                    value = sum(count),
                    variable = "Total haplotaggable reads",
                    sample = unique(data$sample)
                )
        }else{
            data_count <- tibble::tibble(
                value = NA, 
                variable = "Total haplotaggable reads",
                sample = unique(data$sample)
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
