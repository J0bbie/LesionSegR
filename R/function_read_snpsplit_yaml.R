#' @title Import and clean SNPsplit reporting yamls
#' @description Import and clean SNPsplit reporting yamls
#' @param x (character): Path to flagstats file.
#' @return (tibble): Tibble containing the SNPsplit reporting yaml.
#' @importFrom dplyr %>%
#' @export
read_snpsplit_yaml <- function(x) {

    # Check input.
    checkmate::assertFile(x)

    futile.logger::flog.info(glue::glue("Importing {length(x)} SNPSplit yaml. files."))

    # Read flagstats.
    data_snpsplit <- dplyr::bind_rows(future.apply::future_lapply(x, function(sample_path){

        tibble::as_tibble(yaml::read_yaml(sample_path)$Tagging) %>%
            dplyr::mutate(
                sample = basename(sample_path) %>% stringr::str_replace_all(pattern = "_.*", replacement = "")
            )
    }))
}
