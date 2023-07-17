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
        # Import and clean.
        readr::read_csv(sample_path, col_names = FALSE, show_col_types = FALSE) %>%
            dplyr::mutate(
                X1 = base::trimws(X1),
                X2 = as.integer(base::trimws(gsub(":.*", "", X2))),
                sample = basename(sample_path) %>% stringr::str_replace_all(pattern = "_.*", replacement = "")
            ) %>%
            dplyr::select(variable = X1, value = X2, sample)
    }))) %>%
        tidyr::spread(key = variable, value = value)

    # Return.
    return(data_flagstats)
}
